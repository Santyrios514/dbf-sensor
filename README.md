# sensor-lastre · Barracuda DBF 2026-27 (UPB)

Volumen utilizable y lastre del sensor remolcado. Para cada configuración (geometría del
`.ork` más un barrido de parámetros) y cada material de relleno calcula:

- el volumen útil para lastre;
- la masa, el CG y el margen estático;
- la ventana de longitudes de lastre que cumple el SM;
- el techo de estabilidad por geometría, $SM_\infty$;
- el trim respecto al punto de remolque.

Especificaciones: `SPEC_lastre_sensor.md` (v1) y `SPEC_lastre_sensor_v2_cambios.md` (v2; gana
ante conflictos). Geometría de referencia: `modelos/analisis_vol_int.ork` (OpenRocket 24.12).

| Fase | Contenido | Estado |
|---|---|---|
| A2 | lector del `.ork` sin JVM, perfiles de transición, aleta freeform, tests 1–4, 6 y 8 | ✅ |
| B2 | script 2 con aletas freeform, popa abierta, $SM_\infty$, figuras nuevas, test 7 | ✅ |
| C2 | `or_bridge.py` y script 1 sobre el `.ork` real, exportación por componente, `--regresion` (test 5) | ✅ |
| D2 | barrido `h_tip` × `extension`, comparación OpenRocket frente al modelo interno | ✅ |

## Instalación

```bash
pip install -e ".[dev]"            # numpy, scipy ≥ 1.12, pandas, pyyaml, matplotlib, pytest
pip install -e ".[openrocket]"     # script 1: orlab + JPype
python -m orlab fetch 24.12        # jar verificado por sha256 (o ORLAB_JAR=/ruta/OpenRocket-24.12.jar)
```

El script 1 requiere JDK 17 o 21. Los tests que usan OpenRocket (`tests/test_openrocket.py`) se
saltan solos si no hay `java` o jar.

## Uso

```bash
python scripts/01_orlab_export.py --regresion    # compara con OpenRocket (v2 §1.4, M = 0.3)
python scripts/01_orlab_export.py                # barrido → data/or_*.csv
python scripts/02_volumen_lastre.py              # usa data/or_*.csv
python scripts/02_volumen_lastre.py --sin-orlab  # perfiles analíticos + Barrowman interno
pytest
```

Opciones comunes: `--solo id1,id2`, `--rellenos plomo_macizo,W90_macizo`, `--sin-figuras`,
`--dir-datos`, `--guardar-ork` (script 1).

### Configuración (`config/config.yaml`)

- `barrido.parametros` acepta **cualquier ruta** del YAML: primero se busca en las secciones y
  luego dentro de `geometria_base`, p. ej. `aletas.h_tip.valor` o `electronica.masa_g`. El modo
  puede ser `cartesiano` o `zip`, y el id queda como `L350_D65_h_tip30_extension20`.
- `configuraciones.base_ork` reproduce el `.ork`. El script 1 verifica contra ella que la
  geometría leída por OpenRocket coincide con la configuración antes de modificar nada.
- `longitud_referencia: cuerpo` interpreta `L_mm` como longitud del cuerpo (350 mm). La longitud
  total con aletas se reporta aparte (`L_total_mm` = 390 mm en la base).
- Los datos que el equipo debe confirmar están marcados `# PENDIENTE` (ver más abajo).

## Salidas (`data/`)

| Archivo | Contenido |
|---|---|
| `or_resumen.csv` | script 1: geometría aplicada, r_LE, r_tip, CP, C_Nα (y por diferencias), CD, masa y CG en vacío, envolvente, advertencias, versiones, sha256 |
| `or_componentes.csv` | C_Nα, CP, masa y CG por componente (nariz, cuerpo, cola, aletas, total) |
| `or_aletas.csv` | polígono global de la aleta tal como lo aceptó OpenRocket (`origen` = punto \| superficie) |
| `or_perfiles.csv` | r_e(x) muestreado de OpenRocket |
| `resultados.csv` | una fila por configuración × relleno. `_geo` = límite geométrico; sin sufijo = ℓ_SM,max |
| `optimo.csv` | **ranking por el criterio del equipo: máxima masa de lastre con SM ≥ SM_min**, por relleno, solo entre las configuraciones que caben en la bahía (`optimizacion.exigir`). `limitante_masa` dice si la masa la limita el SM o la geometría |
| `presupuestos.csv` | para cada masa objetivo de lastre: ℓ delantero, ℓ trasero, `no_cabe`, CG, SM, trim y costo |
| `configuraciones.csv` | una fila por configuración: L total, D, **D acostado** y **D parado**, sección en vuelo (h_env, w_env), geometría de aleta, x_CP y ruta del dibujo |
| `masas_capas.csv`, `curvas/`, `ventanas.csv` | detalle por capa, curvas ℓ → CG/SM/trim y todos los intervalos que cumplen el SM |

Figuras (`figs/`):
- `masa_max__*`, `sm_inf__*`, `sm_max__*` y `masa_SM_min__*`: curvas de diseño frente a
  `h_tip`, una serie por `extension` y una figura por valor de los demás parámetros (p. ej.
  `lastre_trasero`); el marcador hueco indica que no cabe en la bahía.
- `perfil__*`: perfil, capas, polígono de aleta, lastre útil (delantero y trasero), electrónica y x_CP.
- `dibujo__*`: **dibujo acotado de cada configuración**. Vista lateral: cuerpo, aletas proyectadas
  según la rotación, lastre y electrónica del óptimo, CG y CP. Vista frontal: D máx. con aletas y
  rectángulo de bahía. El relleno dibujado se elige en `salida.relleno_dibujo`.

**Dimensiones guardado.** El sensor se guarda siempre **acostado con las aletas en diagonal**
(`envolvente.rotacion_guardado_deg` = 45°), sin importar la rotación con que vuela. Por eso cada
configuración reporta:

- `D_acostado_mm`: lado del cuadrado que ocupa guardado, $\max(D, \sqrt2\,r_{tip})$ para 4 aletas.
  **Contra este valor se chequea la bahía** (`cabe_alto`, `cabe_ancho`).
- `D_parado_mm`: círculo que tocan las puntas, $2\max(R, r_{tip})$; es lo que ocupa de pie.
- `h_env_mm`, `w_env_mm`: sección con la rotación de **vuelo**, solo como referencia.

Barrowman da el mismo CP para 4 aletas en + y en X, así que la rotación de vuelo no cambia el SM
ni la masa. En el barrido, las filas a 0° y a 45° salen iguales salvo h_env y w_env.

**Costo del relleno** (`costo_usd_kg`, `costo_relleno_usd` en `resultados.csv` y `optimo.csv`):
m_relleno × precio por kg de `costos_usd_kg`. Para el perdigón con epoxy se usa el promedio
ponderado por fracción de masa. Son **precios aproximados de material** (2026, compra minorista de
pocos kg). No incluyen mecanizado, moldes ni envío, así que ajústalos con cotizaciones reales.
- `cg_sm__*`: x_CG(ℓ) y SM(ℓ).

### Banderas

| Bandera | Significado |
|---|---|
| `inviable_geo` | ℓ_geo ≤ 0 |
| `SM_inalcanzable` | ningún ℓ ∈ [0, ℓ_geo] cumple el SM con este relleno |
| `SM_inalcanzable_por_geometria` | $SM_\infty < SM_{min}$: ni un lastre infinitamente denso basta. La solución son las aletas o el CP, no el lastre. No se reporta ℓ_SM |
| `margen_SM_bajo` | $SM_\infty - SM_{min}$ < `umbral_margen_bajo_cal` |
| `no_unimodal` | x_CG(ℓ) tiene más de un valle o la ventana está partida (ver `ventanas.csv`) |
| `inestable_respecto_remolque` | x_CP ≤ x_T |
| `trim_no_lineal` | \|α_trim\| > `alpha_lineal_max_deg`: el amarre está tan lejos del CG que la fórmula lineal no vale |
| `dif_CP_alta` | \|x_CP,OR − x_CP,interno\| > 5 % de L_cuerpo |
| `dif_masa_alta` | la masa de algún componente difiere de OpenRocket en más de 3 % (aletas: 2 %). Detecta, por ejemplo, una transición marcada *Filled* |

`cabe_largo`, `cabe_alto` y `cabe_ancho` quedan vacíos si falta el dato de la bahía.

## Criterio de optimización: máxima masa con SM ≥ 1

El criterio del equipo es **maximizar la masa de lastre sin que el SM baje de `SM_min_cal` = 1.0**.
Para cada configuración y relleno, el óptimo se construye así:

1. **Tapón delantero** desde la nariz hasta ℓ_SM,max: lo que permitan el SM o la geometría
   (la electrónica va detrás y debe terminar antes de la cola). Este lastre suma masa **y**
   adelanta el CG.
2. **Tapón trasero** (`lastre.lastre_trasero: true`), solo si el delantero llegó a ℓ_geo:
   detrás de la electrónica, desde su final hacia popa, hasta que SM = SM_min o se acabe la
   cavidad. Suma masa pero atrasa el CG.

El orden importa. Mientras el tapón delantero pueda crecer, pasar volumen de atrás hacia adelante
sube el SM con la misma masa; por eso en el óptimo el delantero está lleno o no hay trasero. Como
el tapón trasero queda detrás del CG, SM(ℓ₂) es monótona y el límite se halla con un `brentq`.

## Resultados (OpenRocket 24.12, M = 0.0897, SM_min = 1.0, plomo en la cola permitido)

La bahía **no es un dato**: el ejercicio sirve para fijarla. `optimo.csv` reporta para cada
diseño el alto y el ancho que exige (`h_env_mm`, `w_env_mm`, sección con aletas incluidas). Solo
se filtra por el alto conocido (122 mm). Barrido: `h_tip` × `extension` × rotación (+ a 0°, X a 45°).

**Frontera masa / SM / tamaño de bahía (plomo):**

| Diseño (h_tip / ext. / rotación) | Masa de lastre | SM | Bahía (alto × ancho) | Comentario |
|---|---|---|---|---|
| 30 / 40 / X | 6 816 g | 1.00 | **85.9 × 85.9 mm** | bahía mínima; sin margen de SM |
| **40 / 40 / X** | **6 944 g** (cavidad llena) | **1.46** | **100.0 × 100.0 mm** | **recomendado**: masa máxima con margen |
| 50 / 40 / X | 6 944 g | 1.78 | 114.2 × 114.2 mm | puesto 1 del ranking; la bahía extra solo compra SM |
| 30 / 40 / + | 6 816 g | 1.00 | 121.5 × 121.5 mm | dominado por la versión X |
| 20 / 40 / + (base) | 4 116 g | 1.00 | 101.5 × 101.5 mm | aleta chica: el SM limita |

Con W90, los mismos diseños dan 10 272 g (30/40) y 10 410 g (cavidad llena).

Costo aproximado del relleno del diseño recomendado (40/40/X):

| Relleno | Masa | US$/kg | Costo material |
|---|---|---|---|
| Plomo | 6.94 kg | 4 | ≈ US$ 28 |
| Perdigón + epoxy | 4.57 kg | 7.4 | ≈ US$ 34 |
| Bismuto | 5.99 kg | 25 | ≈ US$ 150 |
| W90 | 10.41 kg | 150 | ≈ US$ 1 560 (más mecanizado) |
| Acero | 4.81 kg | 3 | ≈ US$ 14 |

El diámetro máximo con aletas (circunferencia de las puntas) es $2 r_{tip}$ y no depende de la
rotación: 141.5 mm con h_tip = 40. En X la bahía necesaria es menor, 100 × 100 mm, porque las
puntas quedan en las esquinas del rectángulo.

1. **La masa máxima físicamente posible es la cavidad llena:** 612 cm³, es decir 6.94 kg de plomo
   o 10.4 kg de W90. Se alcanza con h_tip ≥ 40 mm y tapón trasero. Por encima de eso, más aleta
   solo aumenta el SM.
2. **Las aletas en X dominan a las de +:** con la misma aleta, a 45° la sección mide
   $\sqrt2\,r_{tip}$ en vez de $2 r_{tip}$ (−29 %). Con la misma bahía caben aletas más grandes.
3. **Para fijar el ancho:** una bahía cuadrada de **100 mm** (más la holgura de montaje y el
   medio espesor de aleta, 0.9 mm) admite el diseño 40/40/X. Da la masa máxima con SM 1.46, un
   margen de ~0.46 cal frente a la incertidumbre del CP. Bajar a 86 mm cuesta solo 128 g
   (−1.8 %) pero deja SM = 1.00 justo.
4. **Diseño base (h_tip = 20):** el SM limita a 4.1 kg con plomo y no admite tapón trasero.
5. **Longitud:** L_total va de 350 a 390 mm, más que el `largo_max_mm` = 340 de la v1. El largo
   de la bahía también sale de aquí: 390 mm con extensión 40, 370 mm con extensión 20.

### Advertencias sobre el óptimo

- **SM = 1.00 exacto no deja margen para la incertidumbre del CP.** El CP de OpenRocket es optimista
  (advertencia *jagged*) y el Barrowman interno difiere 0.5–3 mm. Con el tapón trasero, cada 0.1 cal
  de margen cuesta ≈ 230 g de plomo. Para reservarlo, basta subir `SM_min_cal`, por ejemplo a 1.15.
- **Trim y punto de remolque.** Con 7–10 kg, el peso domina la aerodinámica: $m g$ ≈ 72–106 N
  frente a $q S_{ref} C_{N\alpha}$ ≈ 11 N/rad (6.7 a 9.9 veces menos). Por eso el amarre va
  **en el CG** por defecto (`remolque.x_T: en_CG`), con trim estático 0°. Lo que importa es
  cuánto error admite su posición:
  $$|x_{CG}-x_T| \le \frac{\alpha_{max}\,q\,S_{ref}\,C_{N\alpha}\,(x_{CP}-x_{CG})}{m\,g}$$
  Ese valor es `tol_amarre_mm`, con `alpha_trim_max_deg` = 5°. A 30 m/s sale **±1.24 mm** con plomo
  y ±0.85 mm con W90 en el diseño 40/40/X, y ±1.96 / ±1.34 mm en 50/40/X. La tolerancia crece con
  $V^2$ y con el margen $(x_{CP}-x_{CG})$, y cae como $1/m$. El amarre debe ser **ajustable**
  (riel o perforaciones a ±5 mm de `x_T_trim_cero_mm` ≈ 152–156 mm) y calibrarse pesando el
  sensor ya armado. Si se fija x_T lejos del CG, `trim_no_lineal` avisa cuando |α| supera
  `alpha_lineal_max_deg` (15°), porque fuera de ese rango la fórmula lineal no da ángulos físicos.
- **Cargas:** un sensor de 7–10 kg multiplica las cargas en el cable, el winch, la compuerta y el
  casco de PLA durante el despliegue. Hay que verificarlas contra el presupuesto de masa (v2 §10.4).

### OpenRocket frente al modelo interno

- **Masas por componente:** difieren menos de 0.3 % en todo el barrido. Las pequeñas diferencias
  vienen de que la pared es bicapa y la del `.ork` es equivalente.
- **Perfil:** los perfiles analíticos coinciden con OpenRocket en menos de 0.05 mm, en las 6
  formas con y sin recorte.
- **CP:** difieren como máximo 11.6 mm (h_tip 20, extensión 0); con h_tip ≥ 30 mm, menos de 3 mm.
  Siempre menos que el 5 % de L. Ver `dx_CP_or_vs_interno_mm`.

## Hipótesis y limitaciones

- El CP de Barrowman/OpenRocket vale para flujo subsónico, ángulos pequeños y flujo libre. No
  incluye la estela del avión, la interferencia del cable ni la catenaria.
- **La advertencia de perfil irregular** (*jagged*) viene del escalón P4–P5. OpenRocket, y el
  método interno que lo replica, cortan la cuerda solo contra los puntos de la aleta, sin la
  franja que sigue la cola. Además la porción detrás de la base trabaja en la estela separada
  de la base. **El CP reportado es optimista.** Para quitar el escalón: P4 → P5 en línea recta
  hasta la esquina de la base, o cerrar la popa.
- El trim es estático y lineal. La dinámica del cuerpo remolcado queda fuera del alcance.
- La electrónica es una masa uniforme de diámetro completo y se coloca "detrás del lastre".
  Con ℓ pequeño eso la pone dentro de la nariz, donde no cabe (ver pendientes).

## Desviaciones respecto a las especificaciones

1. **v1 §3.10:** k_n de la potencia debe ser 2n/(2n+1), y la razón de diámetros de la cola cónica
   estaba invertida. La v2 C2 lo corrige con la forma general.
2. **Transiciones recortadas.** OpenRocket recorta por defecto las transiciones elipsoide,
   potencia y Haack. La tabla de la v2 §3.1 describe solo el perfil sin recortar, así que se
   agregó `cola.recortada` (por defecto `true`, como OpenRocket). La parabólica del `.ork` no se ve afectada.
3. **Aletas en el Barrowman interno.** Las cuerdas se cortan contra la polilínea de puntos, no
   contra el polígono cerrado por la superficie (así lo hace FinSetCalc). Cerrando por la
   superficie, el C_Nα salía hasta 12 % bajo y el CP hasta 22 mm corrido.
4. **Envolvente general.** Se usa $h_{env} = \max(R, r_{tip}\max\cos\varphi) + \max(R, -r_{tip}\min\cos\varphi)$,
   que coincide con la fórmula de la v2 para n par y además vale para n impar.
5. **$SM_\infty$** se evalúa en ℓ* de cada relleno, como indica la spec. Ojo: con ρ → ∞, ℓ* → 0
   y el techo tiende a (x_CP − x_b0)/D, así que $SM_\infty$ depende de la región elegida. El techo
   físico con un material dado es `SM_max_alcanzable_cal`.
6. **Solape de aletas en el eje:** con r_in = 0 son ≈ 0.17 g, no < 0.1 g (sigue siendo < 1 % de las aletas).
7. Las columnas `cabe_diametro` pasan a `cabe_alto` y `cabe_ancho`. `resultados.csv` agrega
   `L_total_mm`, `h_env_mm` y `w_env_mm`.

## Datos pendientes del equipo (v2 §10)

| Dato | Valor por defecto | Qué cambia |
|---|---|---|
| Electrónica (largo, diámetro, masa, posición) | 100 mm, 150 g, detrás del lastre | ℓ_SM,min, ℓ_geo, x_CG. **Conviene exigir que empiece en el tramo cilíndrico** |
| Punto de remolque y herraje | x_T en el CG (ajustable), 15 g en x = 40 mm | `tol_amarre_mm`, `alpha_trim_*`, `x_T_trim_cero_mm` |
| SM_min | 1.0 cal (criterio: máxima masa con SM ≥ 1) | óptimo de masa, ventana, banderas. Conviene agregar margen por incertidumbre del CP |
| Presupuesto de masa | 250–1500 g | `presupuestos.csv` |
| Material y espesor de aletas | MAT_PARED_EQ, 1.8 mm | m_aletas, CG en vacío |
| Popa abierta o cerrada | abierta | mamparo de base, masa |
| Ancho y largo de la bahía | **salida del ejercicio** (`w_env_mm`, `h_env_mm`, `L_total_mm`) | recomendado: 100 × 100 mm de sección y 390 mm de largo (diseño 40/40/X) |

## Notas de fabricación

- **Plomo fundido (327 °C) dentro de PLA o PETG no es viable:** su transición vítrea está en
  ~60–80 °C. El lastre se funde o mecaniza aparte y luego se inserta. Alternativas: tungsteno
  W90 (más masa en menos longitud) o perdigón con epoxy en capas delgadas. Manipular plomo con
  guantes y no lijarlo.
- **El lastre debe quedar retenido estructuralmente** frente al despliegue, la retracción y el
  tirón del cable.
