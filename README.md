# Secado de leña de eucalipto en Talca: modelo 2D transiente

Proyecto semestral IIQ2003 – Fenómenos de Transporte (PUC Chile, 2026-2).

## 1. Descripción y relevancia sustentable

En Talca, la combustión residencial de leña es la principal fuente de genreación de calor que se ocupan en las casas. La leña húmeda quema peor y emite más: según el Ministerio del Medio Ambiente, un calefactor de doble cámara pasa de 5 g MP2,5/kg con leña a 25 % de humedad a 10 g/kg con leña a 40 %. La NCh 2907 define como leña seca la que tiene ≤ 25 % de humedad (base seca).

Este repositorio modela cuánto tarda un trozo de leña de *Eucalyptus globulus* en llegar a 25 % según las condiciones de verano e invierno de Talca, y cómo influyen el diámetro y las propiedades de transporte.

## 2. Sistema físico y modelo

Se modelará un trozo de laño como un cilindro de radio `R` y largo `L`, con humedad inicial uniforme `M0`, expuesto a aire con temperatura y humedad relativa constantes. El agua difunde hacia la superficie (lateral y extremos) y se pierde al aire.

Se simula un cuarto del corte axial (0 ≤ r ≤ R, 0 ≤ z ≤ L/2).

**Ecuación**:

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
