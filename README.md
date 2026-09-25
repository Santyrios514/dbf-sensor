# sensor-lastre · Barracuda DBF 2026-27 (UPB)

Volumen utilizable y lastre del sensor remolcado: para cada configuración geométrica
(L × D) y cada material de relleno calcula el volumen útil para lastre, la masa, el CG,
el margen estático, la ventana de longitudes de lastre que cumple el SM y el trim
respecto al punto de remolque.

Especificación de referencia: `SPEC_lastre_sensor.md` v1.0.

| Fase | Contenido | Estado |
|---|---|---|
| A | configuración, perfiles analíticos, erosión multicapa, tests analíticos (8.1) | ✅ |
| B | script 2 en modo `--sin-orlab`, Barrowman interno, tests del modelo (8.2) | ✅ |
| C | `or_bridge.py` y script 1 (OpenRocket vía orlab/JPype) | pendiente: falta `modelos/sensor_base.ork` |
| D | integración 1 → 2, comparación OpenRocket vs. interno (8.3) | pendiente |

## Instalación

```bash
pip install -e ".[dev]"            # numpy, scipy ≥ 1.12, pandas, pyyaml, matplotlib, pytest
pip install -e ".[openrocket]"     # solo para el script 1: orlab + JPype
```

Para el script 1 (fase C) se necesita JDK 17 o 21 y el jar de OpenRocket: `python -m orlab fetch`,
o bien la ruta en la variable `ORLAB_JAR` o en `openrocket.jar` del YAML.

## Uso

```bash
# sin OpenRocket: perfiles analíticos + Barrowman interno
python scripts/02_volumen_lastre.py --config config/config.yaml --sin-orlab

# con la salida del script 1
python scripts/02_volumen_lastre.py --config config/config.yaml \
    --or-resumen data/or_resumen.csv --or-perfiles data/or_perfiles.csv

# opciones: --rellenos plomo_macizo,W90_macizo  --solo L270_D60,L350_D65  --sin-figuras
pytest
```

Todos los parámetros están en `config/config.yaml`. El barrido `L_mm × D_mm` genera ids
`L{L}_D{D}`, y `configuraciones` agrega casos explícitos con *deep merge* sobre `geometria_base`
(también acepta overrides de `pared`, `electronica`, `lastre`, etc.).

## Salidas (`data/`)

| Archivo | Contenido |
|---|---|
| `resultados.csv` | una fila por configuración × relleno. Las columnas `_geo` corresponden a ℓ_geo (límite geométrico); las demás, a ℓ_SM,max (**volumen utilizable real**) |
| `presupuestos.csv` | para cada masa de lastre objetivo: ℓ necesaria, `no_cabe`, CG, SM y trim. **Es la salida de diseño más útil** |
| `masas_capas.csv` | masa y CG de cada capa de pared por estación |
| `curvas/{id}__{relleno}.csv` | barrido ℓ → V_b, m_b, m_total, x_CG, SM, α_trim |
| `ventanas.csv` | (complementario) todos los intervalos de ℓ que cumplen el SM |

Figuras (`figs/`): `cg_sm__*` (x_CG(ℓ) y SM(ℓ) con la ventana sombreada), `perfil__*`
(perfil, capas y región de lastre) y `mapa__{relleno}` (mapas de calor L × D).

### Banderas (`banderas`, separadas por `;`)

| Bandera | Significado |
|---|---|
| `inviable_geo` | ℓ_geo ≤ 0: no cabe lastre entre x_b0 y la primera restricción |
| `SM_inalcanzable` | ninguna ℓ ∈ [0, ℓ_geo] cumple SM_min (ni SM_max, si está definido) |
| `no_unimodal` | x_CG(ℓ) tiene más de un valle, o la ventana está partida en varios intervalos (ver `ventanas.csv`) |
| `inestable_respecto_remolque` | x_CP ≤ x_T: no hay efecto veleta respecto al amarre |
| `dif_CP_alta` | \|x_CP,OR − x_CP,interno\| > `verificacion.dif_CP_max_frac_L`·L |

`cabe_largo` y `cabe_diametro` son banderas de envolvente (sección 3.12): se reportan, no filtran.

## Modelo (resumen)

- **Pared multicapa:** erosión morfológica exacta del meridiano con un disco de radio T_k
  (espesor normal constante), por estación y con el mínimo en una ventana ±T_max en las uniones.
  Maneja la punta de la nariz, las esquinas convexas y la popa cerrada.
- **Integrales:** Simpson compuesto sobre una malla uniforme (`dx_mm`); integrales acumuladas
  para ∀_b(ℓ) y Φ_1(ℓ).
- **CP interno:** Barrowman en la forma general que usa OpenRocket para cuerpos de revolución,
  C_Nα = 2(A_a − A_f)/A_ref y x = x_0 + (ℓA_a − ∀)/(A_a − A_f), más la fórmula de aletas de la
  sección 3.10.
- **Ventana:** barrido denso + `minimize_scalar` (ℓ*) + `brentq` en cada cambio de signo de
  max(x_CG − x_max, x_min − x_CG).

### Desviaciones respecto a la especificación v1.0 (y por qué)

1. **k_n de la nariz de potencia.** La spec da n/(n+1), que es incorrecto: para n = 1 (cono)
   da 0.5 en vez de 2/3. Lo correcto es 2n/(2n+1). En lugar de una tabla de k_n se usa la forma
   general 1 − ∀_n/(A L_n), que reproduce 2/3, 1/3, 2n/(2n+1) y ≈ 0.466 (ogiva esbelta).
2. **Posición del CP de la cola cónica.** En la fórmula de la spec la razón de diámetros está
   invertida: Barrowman usa d_F/d_R (delantero/trasero), no d_a/D. Con d_a/D = 0.6 el error es
   0.083·L_t (5 mm con L_t = 60 mm). Se usa la forma general, equivalente a la fórmula correcta.
3. **Factor de llenado con mamparos:** ∀_b = f_ll·(∫A_i − ∀_mamparos), y no f_ll∫A_i − ∀_mamparos.
   Solo difiere si f_ll < 1 y hay mamparos en la región de lastre.
4. **`no_unimodal`.** Con la electrónica detrás del lastre, x_CG(ℓ) sube unas décimas de mm en
   los primeros milímetros (en la punta, ρ_b A_i (x_b0 − x_CG) + m_e > 0) y luego baja. Esa
   joroba es física y aparece en todas las configuraciones, pero no parte la ventana, así que
   no levanta la bandera.
5. **ℓ_SM,min y ℓ_SM,max** se definen como el ínfimo y el supremo del conjunto que cumple el SM
   (coinciden con la spec en el caso unimodal). Con SM_max o ventanas partidas, los intervalos
   completos quedan en `ventanas.csv`.
6. **`mamparos` vacío por defecto:** el ejemplo de la spec, con `espesor_mm: 0`, contradice la
   regla "espesores > 0".

## Hipótesis y limitaciones

- El CP de Barrowman/OpenRocket vale para flujo subsónico, ángulos pequeños y flujo libre. No
  incluye la estela del avión, la interferencia del cable ni la curvatura de la catenaria.
- El trim de la sección 3.11 es estático y lineal. La estabilidad dinámica del cuerpo remolcado
  (modo péndulo, acoplamiento con el cable) queda fuera del alcance. **Con los supuestos iniciales
  los α_trim en ℓ_SM,max salen de ~14° a más de 200°, muy fuera del rango lineal:** esos valores
  indican que el amarre está mal ubicado, no un ángulo real.
- La electrónica se modela como masa uniforme que ocupa todo el diámetro.
- En la interferencia aleta–cuerpo se usa R del cuerpo, aunque con `borde_salida_desde_popa_mm: 0`
  la raíz de las aletas queda casi toda sobre el cono de cola (de radio menor).

## Notas de interpretación y fabricación

- **El volumen casi nunca es la restricción activa.** Llenar el volumen geométrico da entre
  ~2 y ~8 kg de plomo. La restricción real suele ser el presupuesto de masa: MTOW, tensión del
  cable, torque del winch y carga sobre la compuerta. Por eso `presupuestos.csv` es la salida
  de diseño más útil.
- **Plomo fundido (327 °C) dentro de un casco de PLA o PETG no es viable:** la transición vítrea
  de esos plásticos está en ~60–80 °C. El lastre macizo se funde o mecaniza aparte (molde metálico
  o de grafito) y luego se inserta y se pega. Alternativas: aleación de tungsteno mecanizada (más
  masa en menos longitud, por tanto CG más adelantado) o perdigón con epoxy vertido en capas
  delgadas para controlar la exotermia. Manipular plomo con guantes y no lijarlo.
- **Anclaje del lastre:** el lastre debe quedar retenido estructuralmente, no solo pegado, frente
  a las cargas de despliegue y retracción y al tirón del cable.
- **El punto de remolque domina el trim.** El peso (mg ≈ 25 N con 2.6 kg) supera ampliamente la
  rigidez aerodinámica q·S_ref·C_Nα ≈ 9 N/rad (30 m/s, D = 60 mm). Para α_trim ≤ 5° el amarre debe
  quedar a pocos milímetros del CG:
  |x_CG − x_T| ≲ α_max·q·S_ref·C_Nα·(x_CP − x_T)/(m g). Revisar `x_T_trim_cero_mm`.
