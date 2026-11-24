
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq
from scipy.special import j0, j1


# ----------------------------------------------------------------------
# 1. Humedad de equilibrio (Realizada con isoterma de sorción Hailwood-Horrobin)
# ----------------------------------------------------------------------
#T_C : temperatura del aire
#HR  : humedad relativa
#Se pasa el T_C de °C a F

def humedad_equilibrio(T_C, HR):
    
    T = T_C * 9.0 / 5.0 + 32.0
    W = 330.0 + 0.452 * T + 0.00415 * T**2
    K = 0.791 + 4.65e-4 * T - 1.34e-6 * T**2
    K1 = 6.34 + 7.75e-4 * T - 9.3e-5 * T**2
    K2 = 1.09 + 2.84e-2 * T - 9.04e-5 * T**2
    h = HR
    termino1 = K * h / (1.0 - K * h)
    termino2 = (K1 * K * h + 2.0 * K1 * K2 * K**2 * h**2) / (
        1.0 + K1 * K * h + K1 * K2 * K**2 * h**2)
    return (1800.0 / W) * (termino1 + termino2) / 100.0


# ----------------------------------------------------------------------
# 2. Coeficiente de transferencia de masa efectivo (base M)
# ----------------------------------------------------------------------
#Se calcula primero h_aire con Sherwood
#y se convierte a h_m = h_aire * rho_v_sat / (rho_s * dM/dHR)
#u: velocidad del viento 
#d: diámetro del trozo
#rho_s: densidad de la madera seca
#T_C : temperatura del aire
#HR: Humedad relativa
    
def h_masa_efectivo(T_C, HR, u, d, rho_s, D_AB, nu):
    
    Nu = d * 
    Re = u * d / Nu
    Sc = Nu / D_AB
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
# Retorna dict con: t [s], M_prom(t), M_max(t) (centro), t_prom, 
# t_max_obj (tiempos en que se alcanza M_obj), campo final M(r,z) y malla (rc, zc)
# y opcionalmente snapshots de M para graficar.
def resolver(R, L, D, h_m, M0, M_eq, Nr, Nz, dt, t_max, M_obj=0.25,
             guardar_cada=None):
                 
    dr = R / Nr
    dz = (L / 2.0) / Nz
    rc = (np.arange(Nr) + 0.5) * dr    # Centros de celda en r
    zc = (np.arange(Nz) + 0.5) * dz    # Centros de celda en z
    idx = lambda i, j: i * Nz + j    # Numeración de la celda (i, j)

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
        k = np.argmax(serie <= M_obj)
        if serie[k] > M_obj:
            return np.nan
        return t[k - 1] + (serie[k - 1] - M_obj) / (serie[k - 1] - serie[k]) * dt

    return dict(t=t, M_prom=Mp, M_centro=Mx, t_prom=cruce(Mp),
                t_centro=cruce(Mx), campo=M.reshape(Nr, Nz), rc=rc, zc=zc,
                snaps=snaps)
