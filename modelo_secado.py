"""
Modelo de secado de leña: difusión de humedad 2D (r, z) transiente.

Ecuación:   dM/dt = D_ef * ( d2M/dr2 + (1/r) dM/dr + d2M/dz2 )
Dominio:    0 <= r <= R ;  0 <= z <= L/2   (simetría en z = L/2)
Método:     volúmenes finitos centrados en celda + Euler implícito.

Todas las unidades en SI (m, s, kg/kg). M está en base seca.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq
from scipy.special import j0, j1


# ----------------------------------------------------------------------
# 1. Humedad de equilibrio (isoterma de sorción Hailwood-Horrobin, Simpson 1973)
# ----------------------------------------------------------------------
def humedad_equilibrio(T_C, HR):
    """Contenido de humedad de equilibrio M_eq [kg/kg, base seca].

    T_C : temperatura del aire [°C]
    HR  : humedad relativa [fracción 0-1]
    Coeficientes de Simpson (1973) con T en °F.
    """
    T = T_C * 9.0 / 5.0 + 32.0
    W = 330.0 + 0.452 * T + 0.00415 * T**2
    K = 0.791 + 4.65e-4 * T - 1.34e-6 * T**2
    K1 = 6.34 + 7.75e-4 * T - 9.3e-5 * T**2
    K2 = 1.09 + 2.84e-2 * T - 9.04e-5 * T**2
    h = HR
    term1 = K * h / (1.0 - K * h)
    term2 = (K1 * K * h + 2.0 * K1 * K2 * K**2 * h**2) / (
        1.0 + K1 * K * h + K1 * K2 * K**2 * h**2)
    return (1800.0 / W) * (term1 + term2) / 100.0


# ----------------------------------------------------------------------
# 2. Coeficiente de transferencia de masa efectivo (base M)
# ----------------------------------------------------------------------
def h_masa_efectivo(T_C, HR, u, d, rho_s, D_AB=2.5e-5, nu=1.5e-5):
    """Estima h_m [m/s] en base M para una correlación de Churchill-Bernstein.

    Se calcula primero h_aire [m/s] (base densidad de vapor) con Sherwood,
    y se convierte a base M: h_m = h_aire * rho_v_sat / (rho_s * dM/dHR).
    u: velocidad del viento [m/s]; d: diámetro del trozo [m];
    rho_s: densidad de la madera seca [kg/m3].
    """
    Re = u * d / nu
    Sc = nu / D_AB
    Sh = 0.3 + (0.62 * Re**0.5 * Sc**(1 / 3) / (1 + (0.4 / Sc)**(2 / 3))**0.25
                * (1 + (Re / 282000.0)**(5 / 8))**(4 / 5))
    h_aire = Sh * D_AB / d
    # densidad de vapor saturado [kg/m3] (Magnus + gas ideal)
    p_sat = 610.94 * np.exp(17.625 * T_C / (T_C + 243.04))
    rho_v_sat = p_sat * 0.018015 / (8.314 * (T_C + 273.15))
    dM_dHR = (humedad_equilibrio(T_C, min(HR + 0.05, 0.99))
              - humedad_equilibrio(T_C, max(HR - 0.05, 0.01))) / 0.1
    return h_aire * rho_v_sat / (rho_s * dM_dHR)


# ----------------------------------------------------------------------
# 3. Solver numérico
# ----------------------------------------------------------------------
def resolver(R, L, D, h_m, M0, M_eq, Nr, Nz, dt, t_max, M_obj=0.25,
             guardar_cada=None):
    """Resuelve el problema 2D con volúmenes finitos + Euler implícito.

    Parámetros (unidades): R, L [m]; D [m2/s]; h_m [m/s]; M0, M_eq [kg/kg];
    Nr, Nz: nº de celdas radiales y axiales (mitad del largo);
    dt [s]; t_max [s]; M_obj: humedad objetivo [kg/kg].

    Retorna dict con: t [s], M_prom(t), M_max(t) (centro), t_prom, t_max_obj
    (tiempos en que se alcanza M_obj), campo final M(r,z) y malla (rc, zc),
    y opcionalmente snapshots de M para graficar.
    """
    dr, dz = R / Nr, (L / 2.0) / Nz
    rc = (np.arange(Nr) + 0.5) * dr          # centros de celda en r
    zc = (np.arange(Nz) + 0.5) * dz          # centros de celda en z
    idx = lambda i, j: i * Nz + j            # numeración de la celda (i, j)

    # Conductancia de borde (difusión interna en medio-celda + resistencia del aire)
    U = 1.0 / (dr / (2.0 * D) + 1.0 / h_m)   # lateral  [m/s]
    Uz = 1.0 / (dz / (2.0 * D) + 1.0 / h_m)  # extremo  [m/s]

    N = Nr * Nz
    A = sp.lil_matrix((N, N))
    b = np.zeros(N)
    for i in range(Nr):
        re, rw = (i + 1) * dr, i * dr
        for j in range(Nz):
            k = idx(i, j)
            diag = 0.0
            # cara este (r + dr/2)
            if i < Nr - 1:
                c = D * re / (rc[i] * dr**2)
                A[k, idx(i + 1, j)] += c
                diag -= c
            else:  # r = R : Robin
                c = U * re / (rc[i] * dr)
                diag -= c
                b[k] += c * M_eq
            # cara oeste (r - dr/2); en r=0 el área es 0 (simetría)
            if i > 0:
                c = D * rw / (rc[i] * dr**2)
                A[k, idx(i - 1, j)] += c
                diag -= c
            # cara superior (z + dz/2): en z = L/2 simetría -> sin flujo
            if j < Nz - 1:
                c = D / dz**2
                A[k, idx(i, j + 1)] += c
                diag -= c
            # cara inferior (z - dz/2)
            if j > 0:
                c = D / dz**2
                A[k, idx(i, j - 1)] += c
                diag -= c
            else:  # z = 0 : Robin
                c = Uz / dz
                diag -= c
                b[k] += c * M_eq
            A[k, k] += diag
    A = A.tocsc()

    # Euler implícito: (I - dt A) M^{n+1} = M^n + dt b
    lu = spla.splu(sp.identity(N, format="csc") - dt * A)

    w = np.repeat(rc, Nz)                    # peso de volumen ~ r
    w = w / w.sum()
    M = np.full(N, M0)
    n_pasos = int(round(t_max / dt))
    t = np.zeros(n_pasos + 1)
    Mp = np.zeros(n_pasos + 1)
    Mx = np.zeros(n_pasos + 1)
    Mp[0], Mx[0] = M0, M0
    snaps = {}
    for n in range(1, n_pasos + 1):
        M = lu.solve(M + dt * b)
        t[n] = n * dt
        Mp[n] = w @ M
        Mx[n] = M[idx(0, Nz - 1)]            # centro del trozo (más lento)
        if guardar_cada and n % guardar_cada == 0:
            snaps[t[n]] = M.reshape(Nr, Nz).copy()

    def cruce(serie):
        """Tiempo en que 'serie' baja de M_obj (interpolación lineal)."""
        k = np.argmax(serie <= M_obj)
        if serie[k] > M_obj:
            return np.nan
        return t[k - 1] + (serie[k - 1] - M_obj) / (serie[k - 1] - serie[k]) * dt

    return dict(t=t, M_prom=Mp, M_centro=Mx, t_prom=cruce(Mp),
                t_centro=cruce(Mx), campo=M.reshape(Nr, Nz), rc=rc, zc=zc,
                snaps=snaps)


# ----------------------------------------------------------------------
# 4. Solución analítica (validación): producto cilindro infinito x placa
# ----------------------------------------------------------------------
def _raices(f, n, beta_max=500.0, paso=2e-3):
    """Primeras n raíces positivas de f(beta)=0 (barrido + brentq)."""
    x = np.arange(1e-6, beta_max, paso)
    y = f(x)
    cambios = np.where(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0][:n]
    return np.array([brentq(f, x[c], x[c + 1]) for c in cambios])


def analitica_M_prom(t, R, L, D, h_m, M0, M_eq, n_raices=150):
    """Humedad promedio analítica de un cilindro finito con convección.

    Para D constante, M0 uniforme y borde Robin, la solución es el producto
    de la solución del cilindro infinito (radio R) y la de una placa de
    semiespesor l = L/2. Raíces: cilindro  b*J1(b) = Bi*J0(b);
    placa  b*tan(b) = Bi.
    """
    t = np.atleast_1d(t).astype(float)
    Bi_r, Bi_z = h_m * R / D, h_m * (L / 2) / D
    br = _raices(lambda b: b * j1(b) - Bi_r * j0(b), n_raices)
    bz = _raices(lambda b: b * np.sin(b) - Bi_z * np.cos(b), n_raices)
    Fo_r, Fo_z = D * t / R**2, D * t / (L / 2) ** 2
    cil = (4 * Bi_r**2 / (br**2 * (br**2 + Bi_r**2)))[:, None] \
        * np.exp(-np.outer(br**2, Fo_r))
    pla = (2 * Bi_z**2 / (bz**2 * (bz**2 + Bi_z**2 + Bi_z)))[:, None] \
        * np.exp(-np.outer(bz**2, Fo_z))
    theta = cil.sum(axis=0) * pla.sum(axis=0)
    return M_eq + (M0 - M_eq) * theta
