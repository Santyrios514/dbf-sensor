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
| `presupuestos.csv` | para cada masa objetivo de lastre: ℓ, `no_cabe`, CG, SM y trim. **Es la salida de diseño** |
| `masas_capas.csv`, `curvas/`, `ventanas.csv` | detalle por capa, curvas ℓ → CG/SM/trim y todos los intervalos que cumplen el SM |

Figuras (`figs/`):
- `sm_inf__*`, `sm_max__*` y `masa_SM_min__*`: curvas de diseño frente a `h_tip`, una serie por
  `extension`; el marcador hueco indica que no cabe en la bahía.
- `perfil__*`: perfil, capas, polígono de aleta, lastre útil y x_CP.
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
| `dif_CP_alta` | \|x_CP,OR − x_CP,interno\| > 5 % de L_cuerpo |
| `dif_masa_alta` | la masa de algún componente difiere de OpenRocket en más de 3 % (aletas: 2 %). Detecta, por ejemplo, una transición marcada *Filled* |

`cabe_largo`, `cabe_alto` y `cabe_ancho` quedan vacíos si falta el dato de la bahía.

## Resultados del barrido (OpenRocket 24.12, M = 0.0897, SM_min = 1.5)

Cifras de `data/resultados.csv` con la configuración actual:

| h_tip [mm] | extension [mm] | x_CP [mm] | $SM_\infty$ (Pb) | SM máx. Pb | SM máx. W90 | h_env [mm] |
|---|---|---|---|---|---|---|
| 20 | 40 (**base**) | 157.0 | 1.82 | **1.42** ✗ | 1.54 | 101.5 |
| 20 | 20 | 142.2 | 1.61 | 1.23 ✗ | 1.34 ✗ | 101.5 |
| 30 | 20 | 201.3 | 2.52 | 2.12 | 2.24 | 121.5 |
| 30 | 40 | 217.1 | 2.74 | 2.33 | 2.45 | 121.5 |
| 40 | 0 | 182.3 | 2.24 | 1.86 | 1.97 | 141.5 ✗ |
| 40 | 40 | 250.8 | 3.25 | 2.83 | 2.96 | 141.5 ✗ |

1. **El diseño base no cumple SM = 1.5 con plomo.** Llega a 1.42 como máximo; con W90 alcanza
   1.54, en una ventana estrecha de ℓ entre 39 y 69 mm. La causa es la aleta, no el volumen.
2. **La palanca es `h_tip`.** El C_Nα de la aleta escala con el cuadrado del span. Pasar de 20 a
   30 mm sube el SM máximo de 1.42 a 2.1–2.3 cal.
3. **La bahía limita `h_tip` a 30.3 mm en configuración +.** $h_{env} = 2(r_{LE} + h) \le 122$ mm
   da $h \le 30.26$ mm; con 30 mm quedan 0.5 mm de holgura. **Girar las aletas 45° (X)** reduce
   el alto a $\sqrt2\,r_{tip}$ y permite $h$ ≈ 55 mm, si el ancho de la bahía lo admite.
   Para evaluarlo: `aletas.rotacion_deg: [0, 45]` en el barrido, y definir `envolvente.ancho_max_mm`.
4. **La extensión detrás de la base ayuda mucho.** Sin ella, el CP se adelanta entre 64 y 79 mm.
5. **Longitud:** L_total va de 350 a 390 mm, más que el `largo_max_mm` = 340 de la v1. Hay que
   confirmar el largo real de la bahía.
6. **Trim:** con x_T = 0.15 L, α_trim sale de decenas de grados. El amarre tiene que ir
   prácticamente en el CG (`x_T_trim_cero_mm`, ≈ 105 mm con lastre).

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
| Punto de remolque y herraje | x_T = 0.15 L, 15 g en x = 40 mm | trim (`alpha_trim_*`, `x_T_trim_cero_mm`) |
| SM_min | 1.5 cal | ventana, banderas, figuras |
| Presupuesto de masa | 250–1500 g | `presupuestos.csv` |
| Material y espesor de aletas | MAT_PARED_EQ, 1.8 mm | m_aletas, CG en vacío |
| Popa abierta o cerrada | abierta | mamparo de base, masa |
| Ancho de la bahía (y largo real) | null (340 mm de la v1) | `cabe_ancho`, `cabe_largo`, viabilidad de las aletas a 45° |

## Notas de fabricación

- **Plomo fundido (327 °C) dentro de PLA o PETG no es viable:** su transición vítrea está en
  ~60–80 °C. El lastre se funde o mecaniza aparte y luego se inserta. Alternativas: tungsteno
  W90 (más masa en menos longitud) o perdigón con epoxy en capas delgadas. Manipular plomo con
  guantes y no lijarlo.
- **El lastre debe quedar retenido estructuralmente** frente al despliegue, la retracción y el
  tirón del cable.
