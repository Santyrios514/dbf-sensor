# sensor-lastre · Barracuda DBF 2026-27 (UPB)

Herramienta de dimensionamiento del **lastre del sensor remolcado**. A partir de la geometría del
sensor (un `.ork` de OpenRocket más un barrido de parámetros definido en YAML) calcula cuánto
lastre de plomo cabe, dónde va y cómo queda la estabilidad. Para cada configuración entrega:

- el volumen interior útil para lastre, descontadas las paredes;
- la **masa máxima de lastre con margen estático SM ≥ SM_min** y cómo se reparte entre un tapón
  delantero y uno trasero;
- el CG, el CP, el SM y el techo de estabilidad por geometría $SM_\infty$;
- la ventana de longitudes de lastre que cumple el SM;
- el trim respecto al punto de remolque y la tolerancia de posición del amarre;
- las dimensiones exteriores con aletas: largo total, D acostado y D parado;
- el costo aproximado del relleno, un ranking de configuraciones y un dibujo acotado de cada una.

Sirve para elegir la geometría exterior (largo, diámetro, aletas) y para fijar las dimensiones
de la bahía donde se guarda el sensor.

## Cómo funciona

```
config/*.yaml ──► 01_orlab_export.py ──► data/or_*.csv ──► 02_volumen_lastre.py ──► data/*.csv, figs/
   (barrido)       (OpenRocket: CP,        (perfil, aletas,     (volumen, lastre, CG,        │
                    C_Nα, masas, perfil)    CP por componente)   SM, trim, ranking)          ▼
                                                                          03_comparar_LD.py (resumen L × D)
```

1. **Geometría.** El YAML define nariz, cuerpo, cola (transición) y aletas freeform de cola.
   `barrido.parametros` genera las configuraciones (producto cartesiano o `zip`) y
   `configuraciones` agrega casos explícitos.
2. **Aerodinámica (script 1).** Cada configuración se aplica sobre el `.ork` base mediante
   OpenRocket 24.12 (orlab + JPype). Se exportan el perfil exterior, el polígono de aleta que
   OpenRocket aceptó, el CP y el $C_{N\alpha}$ total y por componente, y las masas en vacío. Antes
   de modificar nada, se verifica que la configuración `base_ork` coincide con el `.ork`.
3. **Volumen útil.** El perfil exterior se erosiona capa por capa (fibra + PLA, espesores del
   YAML) para obtener el radio interior $r_i(x)$. Las integrales acumuladas de área y momento se
   calculan por Simpson:
   $$\forall(\ell)=\int_{x_{b0}}^{x_{b0}+\ell}\pi r_i^2\,dx,\qquad \Phi_1(\ell)=\int_{x_{b0}}^{x_{b0}+\ell}\pi r_i^2\,x\,dx$$
4. **Masa y CG.** Casco, aletas, mamparos, electrónica y masas puntuales dan la masa en vacío
   $m_0$ y su momento $M_0$. Con un relleno de densidad $\rho_b$ y la electrónica detrás del lastre:
   $$x_{CG}(\ell)=\frac{M_0+\rho_b\,\Phi_1(\ell)+m_e\,\bar x_e(\ell)}{m_0+\rho_b\,\forall(\ell)+m_e},\qquad SM=\frac{x_{CP}-x_{CG}}{D}$$
5. **Optimización.** Se busca la máxima masa de lastre con $SM \ge SM_{min}$ (ver abajo).
6. **Trim y amarre.** Con el remolque en $x_T$:
   $$\alpha_{trim}\approx\frac{m g\,(x_{CG}-x_T)}{q\,S_{ref}\,C_{N\alpha}\,(x_{CP}-x_T)}$$
   Por defecto el amarre va en el CG ($\alpha_{trim}=0$). Se reporta el error de posición
   admisible para un trim máximo $\alpha_{max}$:
   $$\text{tol}_{amarre}=\frac{\alpha_{max}\,q\,S_{ref}\,C_{N\alpha}\,(x_{CP}-x_{CG})}{m g}$$

Sin OpenRocket (`--sin-orlab`), el perfil se genera analíticamente con las mismas funciones de
forma de OpenRocket (incluido el recorte de elipsoide, potencia y Haack). El CP sale de un
Barrowman interno: forma general para nariz y cola, y aletas freeform integradas por franjas
como `FinSetCalc`.

## Criterio de optimización: máxima masa con SM ≥ SM_min

Para cada configuración:

1. **Tapón delantero.** Crece desde la nariz ($x_{b0}$, primer punto con $r_i \ge$
   `r_min_util_mm`) hasta donde lo permitan el SM o la geometría (fin de la cavidad, cola,
   electrónica, zonas prohibidas, `fraccion_max_L`). Suma masa y adelanta el CG.
2. **Tapón trasero** (`lastre.lastre_trasero: true`). Solo se agrega si el delantero llegó a su
   límite geométrico. Ocupa el espacio detrás de la electrónica, hacia popa, hasta que
   $SM = SM_{min}$ o se acaba la cavidad. Suma masa pero atrasa el CG. Como $SM(\ell_2)$ es
   monótona, el límite se resuelve con `brentq`.

3. **Tope de masa del sensor** (`lastre.masa_max_sensor_g`). Si la masa total del óptimo supera
   el tope, el lastre se recorta a $m_{max} - m_{vacío}$ y se ubica igual (delantero primero,
   luego trasero). Así el SM queda por encima de SM_min, y si con ese lastre no llega, se marca
   `SM_inalcanzable_con_masa_max`. El tope sale del presupuesto de peso del avión (MTOW menos
   vacío menos margen), documentado en el YAML.

`optimo.csv` filtra las configuraciones que cumplen las banderas de `optimizacion.exigir` y las
ordena por masa total del sensor. Los empates al gramo (tope de masa o cavidad llena) se resuelven por mayor SM, luego
menor D acostado y luego menor D parado. `limitante_masa` indica si la masa la limitó el SM, la
geometría o el tope de masa (`masa_max_sensor`).

## Instalación

```bash
pip install -e ".[dev]"            # numpy, scipy ≥ 1.12, pandas, pyyaml, matplotlib, pytest
pip install -e ".[openrocket]"     # script 1: orlab + JPype
python -m orlab fetch 24.12        # jar verificado por sha256 (o ORLAB_JAR=/ruta/OpenRocket-24.12.jar)
```

El script 1 requiere JDK 17 o 21. Los resultados son deterministas: con el mismo YAML y
OpenRocket 24.12 se obtienen los mismos números.

## Uso

```bash
python scripts/01_orlab_export.py                # OpenRocket → data/or_*.csv
python scripts/02_volumen_lastre.py              # lastre, SM, trim, ranking y figuras
python scripts/02_volumen_lastre.py --sin-orlab  # sin OpenRocket: perfil analítico + Barrowman interno
python scripts/01_orlab_export.py --regresion    # verifica la base contra los valores de referencia
pytest
```

Opciones comunes: `--config` (otro YAML), `--solo id1,id2`, `--rellenos plomo_macizo`,
`--sin-figuras`, `--dir-datos`, `--dir-figuras` y `--guardar-ork` (script 1, guarda cada variante
como `.ork`).

**Comparación de largo y diámetro.** `config/config_barrido_LD.yaml` barre el largo del cuerpo y
el diámetro con nariz y cola de longitud fija. Así cambia solo el tubo, y los radios coinciden en
las uniones.

```bash
python scripts/01_orlab_export.py  --config config/config_barrido_LD.yaml
python scripts/02_volumen_lastre.py --config config/config_barrido_LD.yaml
python scripts/03_comparar_LD.py    --config config/config_barrido_LD.yaml
```

El script 3 resume cada par (L, D): la mejor configuración de aletas, cuántas variantes son
factibles, la caja mínima acostada ($L_{tot}\,D_{ac}^2$) y los kg de lastre por litro de caja.
Sirve para elegir el par con el que se corre el barrido fino.

## Configuración (`config/config.yaml`)

Es la única fuente de parámetros. Unidades: mm, g, kg/m³ y grados.

| Sección | Qué define |
|---|---|
| `materiales`, `costos_usd_kg` | densidades y precio aproximado por kg |
| `rellenos` | material de lastre (plomo macizo); admite rellenos granulares con matriz |
| `geometria_base` | nariz, cuerpo, cola y aletas. Las longitudes pueden ser absolutas o relativas a D, L o L_t |
| `pared` | capas de afuera hacia adentro, por componente o por defecto |
| `mamparos`, `electronica`, `masas_puntuales` | masas internas y su posición |
| `lastre` | inicio, radio mínimo útil, factor de llenado, zonas prohibidas, tapón trasero, tope de masa del sensor |
| `estabilidad` | `SM_min_cal`, `SM_max_cal`, umbral de margen bajo |
| `condiciones_vuelo` | velocidad, altitud (ISA) y Mach (`auto` = V/a) |
| `remolque` | posición del amarre (`en_CG`, absoluta o relativa), $\alpha_{max}$ y límite lineal |
| `envolvente` | largo y alto de la bahía (`null` = sin límite) y rotación de guardado |
| `barrido`, `configuraciones` | parámetros a barrer (cualquier ruta del YAML) y casos explícitos |
| `openrocket` | `.ork` base, jar y nombres de componentes |
| `verificacion`, `numerico`, `salida` | tolerancias, discretización y carpetas de salida |

- `barrido.parametros` acepta cualquier ruta del YAML. Primero se busca en las secciones y luego
  dentro de `geometria_base`, por ejemplo `aletas.h_tip.valor` o `electronica.masa_g`. El id de
  cada caso se arma con los valores, por ejemplo `L350_D65_h_tip30_extension20_rotacion_deg0`.
- `longitud_referencia: cuerpo` interpreta `L_mm` como la longitud del cuerpo (nariz + tubo +
  cola). El largo con aletas se reporta aparte en `L_total_mm`.
- `configuraciones.base_ork` reproduce el `.ork` y sirve de referencia de regresión.
- Los valores marcados `# PENDIENTE` son supuestos que el equipo debe confirmar (electrónica,
  herraje de remolque, material de aletas, popa, presupuesto de masa).

## Salidas

**Datos (`data/`):**

| Archivo | Contenido |
|---|---|
| `or_resumen.csv` | script 1: geometría aplicada, CP, $C_{N\alpha}$, CD, masa y CG en vacío, advertencias, versiones y sha256 |
| `or_componentes.csv` | $C_{N\alpha}$, CP, masa y CG por componente |
| `or_aletas.csv` | polígono de aleta tal como lo aceptó OpenRocket |
| `or_perfiles.csv` | radio exterior $r_e(x)$ muestreado |
| `resultados.csv` | una fila por configuración × relleno. Sufijo `_geo` = límite geométrico; sin sufijo = óptimo con SM ≥ SM_min |
| `optimo.csv` | ranking por masa total del sensor entre las configuraciones factibles |
| `configuraciones.csv` | dimensiones por configuración: L total, D acostado, D parado, sección en vuelo, aleta, x_CP y ruta del dibujo |
| `presupuestos.csv` | para cada masa objetivo: longitudes de tapón, `no_cabe`, CG, SM, trim y costo |
| `ventanas.csv` | todos los intervalos de ℓ que cumplen el SM |
| `masas_capas.csv`, `curvas/` | masa por capa y curvas ℓ → CG, SM y trim |
| `comparacion_LD*.csv` | script 3: resumen por par (L, D) y tabla L × D de masa máxima |

**Figuras (`figs/`):**

- `dibujo__*`: dibujo acotado de cada configuración. La vista lateral muestra el cuerpo, las
  aletas, el lastre y la electrónica del óptimo, el CG y el CP. La vista frontal muestra la
  posición de guardado con D acostado y D parado. El relleno dibujado se elige en
  `salida.relleno_dibujo`.
- `perfil__*`: perfil, capas de pared, polígono de aleta, lastre útil, electrónica y x_CP.
- `cg_sm__*`: $x_{CG}(\ell)$ y $SM(\ell)$.
- `masa_max__*`, `masa_SM_min__*`, `sm_max__*`, `sm_inf__*`: curvas de diseño frente a los
  parámetros del barrido.
- `comparacion_LD__*` (script 3): mapa L × D de masa máxima y dispersión de masa frente a D acostado.

**Dimensiones de guardado.** El sensor se guarda acostado con las aletas en diagonal
(`envolvente.rotacion_guardado_deg` = 45°), sin importar la rotación con la que vuela:

- `D_acostado_mm`: lado del cuadrado que ocupa guardado, $\max(D,\sqrt2\,r_{tip})$ para 4
  aletas. Es el valor que se compara con la bahía (`cabe_alto`, `cabe_ancho`).
- `D_parado_mm`: círculo que tocan las puntas, $2\max(R, r_{tip})$.
- `h_env_mm`, `w_env_mm`: sección con la rotación de vuelo,
  $h_{env}=\max(R,\,r_{tip}\max\cos\varphi)+\max(R,\,-r_{tip}\min\cos\varphi)$.

Con 4 aletas, Barrowman da el mismo CP en + y en X, así que la rotación de vuelo no cambia ni el
SM ni la masa.

**Costo.** `costo_relleno_usd` = masa de relleno × `costos_usd_kg`. Es un precio aproximado de
material, sin mecanizado, moldes ni envío.

### Banderas

| Bandera | Significado |
|---|---|
| `inviable_geo` | no hay longitud útil para lastre ($\ell_{geo} \le 0$) |
| `SM_inalcanzable` | ninguna longitud de lastre cumple el SM con este relleno |
| `SM_inalcanzable_por_geometria` | $SM_\infty < SM_{min}$: ni un lastre infinitamente denso basta. Hay que cambiar aletas o CP, no lastre |
| `margen_SM_bajo` | $SM_\infty - SM_{min}$ < `umbral_margen_bajo_cal` |
| `no_unimodal` | $x_{CG}(\ell)$ tiene más de un valle o la ventana está partida (ver `ventanas.csv`) |
| `SM_inalcanzable_con_masa_max` | con el lastre que permite el tope de masa del sensor, el SM no llega a SM_min |
| `inestable_respecto_remolque` | $x_{CP} \le x_T$ |
| `trim_no_lineal` | $\lvert\alpha_{trim}\rvert$ > `alpha_lineal_max_deg`: la fórmula lineal deja de valer |
| `dif_CP_alta` | $\lvert x_{CP,OR} - x_{CP,interno}\rvert$ > `dif_CP_max_frac_L` · L |
| `dif_masa_alta` | la masa de algún componente difiere de OpenRocket más de lo tolerado (detecta, por ejemplo, una transición marcada *Filled*) |

`cabe_largo`, `cabe_alto` y `cabe_ancho` quedan vacíos cuando la dimensión de la bahía es `null`.

## Verificación

- `01_orlab_export.py --regresion` compara la base con los valores de referencia de
  `verificacion.regresion_or` (masa, CG, CP y estabilidad a M = 0.3).
- En cada corrida se contrastan las masas por componente y el CP de OpenRocket con el modelo
  interno (`dif_masa_or_pct`, `dx_CP_or_vs_interno_mm`).
- Los tests (`pytest`) cubren casos analíticos (cilindros, conos, integrales), la geometría de
  las transiciones y las aletas, la validación del YAML, el modelo de CG, SM y trim, los esquemas
  de salida y la integración con OpenRocket. Estos últimos se saltan si no hay `java` ni jar.

## Hipótesis y limitaciones

- El CP de Barrowman/OpenRocket supone flujo subsónico, ángulos pequeños y flujo libre. No
  incluye la estela del avión, la interferencia del cable ni la catenaria.
- El escalón entre P4 y P5 de la aleta genera la advertencia *jagged* de OpenRocket, y la parte
  de la aleta detrás de la base trabaja en estela separada. **El CP reportado es optimista.**
  Conviene exigir un margen sobre SM_min (por ejemplo `SM_min_cal: 1.15`).
- El trim es estático y lineal. La dinámica del cuerpo remolcado queda fuera del alcance.
- La electrónica se modela como una masa uniforme de diámetro completo, colocada detrás del
  lastre o en una posición fija.
- El lastre no se funde dentro del casco: el plomo funde a 327 °C y el PLA/PETG se ablanda a
  ~60–80 °C. Se funde o mecaniza aparte, se inserta y debe quedar retenido estructuralmente
  frente al despliegue y el tirón del cable.

## Estructura

```
config/     config.yaml (base) y config_barrido_LD.yaml (barrido de largo × diámetro)
modelos/    analisis_vol_int.ork (geometría de referencia, OpenRocket 24.12)
scripts/    01_orlab_export.py, 02_volumen_lastre.py, 03_comparar_LD.py
src/sensor_lastre/
  config.py        carga y validación del YAML, barrido
  perfiles.py      funciones de forma de nariz y transición (OpenRocket)
  aletas.py        aleta freeform: construcción, área, envolvente, Barrowman
  geometria.py     perfil exterior e interior (erosión por capas), integrales acumuladas
  masas.py         masas por componente y masa en vacío
  barrowman.py     CP interno
  estabilidad.py   CG(ℓ), SM, ventana de lastre, tapón trasero, trim
  analisis.py      análisis completo por configuración y relleno
  optimizacion.py  ranking por masa total del sensor
  figuras.py       figuras y dibujos
  or_bridge.py     puente con OpenRocket (único módulo que usa Java)
  ork_xml.py       lector del .ork sin JVM
  esquemas.py      columnas de cada CSV
  materiales.py, atmosfera.py, verificacion.py
tests/      pytest
docs/       api_openrocket.md (API de OpenRocket usada por or_bridge)
```
