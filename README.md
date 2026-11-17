# Secado de leña de eucalipto en Talca: modelo 2D transiente

Proyecto semestral IIQ2003 – Fenómenos de Transporte (PUC Chile, 2026-2).

## 1. Descripción y relevancia sustentable

En Talca y Maule, la combustión residencial de leña es la principal fuente de MP2,5. La leña húmeda quema peor y emite más: según el Ministerio del Medio Ambiente, un calefactor de doble cámara pasa de 5 g MP2,5/kg con leña a 25 % de humedad a 10 g/kg con leña a 40 %. La NCh 2907 define como leña seca la que tiene ≤ 25 % de humedad (base seca).

Este repositorio modela cuánto tarda un trozo de leña de *Eucalyptus globulus* en llegar a 25 % según las condiciones de verano e invierno de Talca, y cómo influyen el diámetro y las propiedades de transporte.

## 2. Sistema físico y modelo

Trozo cilíndrico de radio `R` y largo `L`, con humedad inicial uniforme `M0`, expuesto a aire con temperatura y humedad relativa constantes. El agua difunde hacia la superficie (lateral y extremos) y se pierde al aire.

```
        z
        ^   plano medio (z = L/2): simetría
        |-------------------------+
        |                         |
        |    M(r, z, t)           |  r = R: pérdida al aire (Robin)
        |                         |
        +-------------------------+---> r
      z = 0: pérdida al aire     r = 0: simetría
```

Se simula un cuarto del corte axial (0 ≤ r ≤ R, 0 ≤ z ≤ L/2).

**Ecuación** (derivada con un balance de masa sobre un anillo; detalle en el informe):

∂M/∂t = D_ef · (∂²M/∂r² + (1/r) ∂M/∂r + ∂²M/∂z²)

**Condiciones:**

| Posición | Condición |
|---|---|
| t = 0 | M = M₀ |
| r = 0 | ∂M/∂r = 0 |
| z = L/2 | ∂M/∂z = 0 |
| r = R | −D_ef ∂M/∂r = h_m (M_sup − M_eq) |
| z = 0 | D_ef ∂M/∂z = h_m (M_sup − M_eq) |

`M` es contenido de humedad en base seca [kg agua / kg madera seca]. `M_eq` sale de la isoterma de Hailwood–Horrobin (coeficientes de Simpson, 1973) a la T y HR del escenario. `h_m` se estima con una correlación de Sherwood (Churchill–Bernstein) y se convierte a base `M`.

Supuestos: cilindro homogéneo sin corteza, sin contracción ni agrietamiento, difusión isótropa con `D_ef` constante, T y HR del aire constantes, sin acoplamiento calor-masa, sin lluvia ni radiación solar.

## 3. Método numérico

Se sigue la metodología del curso (Módulo 4, PVIC parabólicos):

- **Esquema FTCS** (forward time, centred space) sobre una malla de **nodos** (r_i = iΔr, z_j = jΔz). La evolución se escribe `M^{n+1} = M^n + Δt · A M^n` con `A` una matriz dispersa para los nodos interiores.
- **Eje (r = 0):** el término (1/r)∂M/∂r se reemplaza por ∂²M/∂r² (límite por simetría), así que no hay singularidad.
- **Bordes:** después de cada paso se corrigen los nodos de borde con diferencias finitas de segundo orden. En r = R y z = 0 se resuelve la condición Robin para el nodo de borde; en z = L/2, ∂M/∂z = 0 se impone como `M_N = (4 M_{N-1} − M_{N-2}) / 3`.
- **Estabilidad:** criterio de Fourier. El paso de tiempo se toma como 90 % del máximo `Δt ≤ 1 / (4D/Δr² + 2D/Δz²)`. Como `D_ef ~ 10⁻¹⁰ m²/s`, este límite es de varias horas (3472 s con la malla elegida), así que FTCS es viable aquí.
- **Verificación cruzada:** el mismo problema se resuelve con volúmenes finitos y Euler implícito (`resolver` en `modelo_secado.py`). Es incondicionalmente estable.
- **Tiempo de secado `t*`:** primer instante en que la humedad promedio (integral ponderada por volumen) o la del centro baja de 25 %, con interpolación lineal.
- **Cantidades de interés ingenieril:** masa de agua en el trozo, caudal de agua perdida por la superficie (flujo con derivada de segundo orden) y balance macroscópico de agua.

El método SOR del Módulo 5 aplica a problemas estacionarios (elípticos), por eso no se usa en este modelo transiente.

## 4. Validación, independencia de grilla y balance

- **Validación analítica:** con `D_ef` constante y `M₀` uniforme el problema tiene solución exacta: producto de la solución del cilindro infinito (raíces de `β J₁(β) = Bi J₀(β)`) y de una placa de semiespesor L/2 (raíces de `β tan β = Bi`). El `t*` de FTCS coincide con el analítico en 0,3 % (16,92 vs. 16,97 días, verano); el esquema implícito también (17,03 días).
- **Independencia de grilla (FTCS):** con celdas cuadradas, `t*` difiere de la malla más fina en 5,3 % (Nr = 5), 1,3 % (Nr = 10) y 0,25 % (Nr = 20). Se usa Nr = 20, Nz = 60.
- **Balance de agua:** el agua perdida por acumulación se compara con la integral del caudal por la superficie. La diferencia baja con el refinamiento (18 %, 8,5 %, 4,0 %, 1,9 % para Nr = 5, 10, 20, 40). Con la malla elegida es de 4 %; el error viene sobre todo del gradiente muy fuerte de los primeros instantes.

Los gráficos están en `figures/` y se generan al ejecutar el notebook.

## 5. Resultados (valores de referencia, ver advertencias)

| Escenario | T [°C] | HR [%] | M_eq [%] | t* promedio [días] | t* centro [días] |
|---|---|---|---|---|---|
| Verano (enero) | 21,5 | 55 | 10,0 | 17 | 34 |
| Invierno (junio) | 9,2 | 84,4 | 18,1 | 28 | 46 |

Sensibilidad (verano): `t*` promedio va de 34 días con `D_ef` = 1,5×10⁻¹⁰ m²/s a 11 días con 4,5×10⁻¹⁰ m²/s, y de 6,5 días (diámetro 6 cm) a 31 días (diámetro 14 cm). Cada trozo del caso base (1,30 kg de madera seca) debe perder 0,71 kg de agua para llegar a 25 %.

**Advertencias:** estos números dependen de supuestos que deben verificarse con una fuente antes de la entrega final: humedad inicial `M₀`, densidad `ρ_s`, HR de verano, velocidad del viento, geometría del trozo, y `D_ef` (medido en *E. nitens*, no en *E. globulus*). Con T y HR constantes el modelo no representa los ciclos diarios.

## 6. Cómo ejecutar

```bash
pip install -r requirements.txt
jupyter notebook secado_lena.ipynb
```

Ejecutar las celdas en orden. Los parámetros están en la sección 1 del notebook, con nombre, unidad y fuente.

## Estructura

```
├── README.md
├── requirements.txt
├── modelo_secado.py      # isoterma, h_m, FTCS (método principal), FV implícito (verificación), solución analítica
├── secado_lena.ipynb     # parámetros, validación, grilla, escenarios, sensibilidad
└── figures/              # gráficos generados
```

## Referencias

Ananías, R. A. & Vergara, F. P. (1991). Un modelo de secado natural de eucalipto (*Eucalyptus globulus* L.). *Ciencia & Investigación Forestal*, 5(2), 229–236.
Simpson, W. T. (1973). Predicting equilibrium moisture content of wood by mathematical models. *Wood and Fiber*, 5(1), 41–49.
Material del curso: Transferencia de momentum, calor y masa computacional (Módulos 4 y 5), PUC.
Demás fuentes de contexto: ver lista de referencias del informe (Entrega 1).
