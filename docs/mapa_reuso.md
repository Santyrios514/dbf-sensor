# Mapa de reúso: `sensor_opt` sobre `sensor_lastre`

`sensor_opt` (spec v3) **importa** estas funciones y clases de `sensor_lastre`; no copia ni parchea
ninguna. Si algo no cubre un caso nuevo, `sensor_opt` lo envuelve o extiende (última columna).

## Configuración y materiales

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/config.py` | `cargar(ruta: str \| Path \| dict) -> Config` | Valida el YAML, resuelve cada caso a SI (`Caso`) y construye la aleta | §2.2: cada cuerpo de la malla se traduce a un YAML de sensor_lastre (`sensor_opt.config.raw_cuerpo`) y se resuelve con esta función. Así la pared, la electrónica, las masas puntuales, el lastre, el remolque y el vuelo salen del código actual sin reinterpretarlos |
| `sensor_lastre/config.py` | `deep_merge(base, override) -> dict` | Fusión en profundidad de dicts | §5: herencia de `base_config` |
| `sensor_lastre/config.py` | `ConfigError`, `FORMAS_COLA` | Error con lista de mensajes; formas válidas | §5: validación de `optimizacion.yaml` |
| `sensor_lastre/config.py` | `Caso.con_aletas(g_aleta) -> Caso` | Copia del caso con otra aleta | §3.2: se reemplaza la aleta marcador por la de cada candidato |
| `sensor_lastre/materiales.py` | `Capa`, `Relleno(nombre, rho_b, costo_usd_kg)` | Capas de pared y relleno | §2.2: pared fibra + PLA, relleno plomo |
| `sensor_lastre/atmosfera.py` | `isa(h) -> EstadoISA(T, p, rho, a)` | Atmósfera estándar | §3.2.1: P y a del flutter |

## Perfiles

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/perfiles.py` | `radio_exterior(x, g) -> ndarray` | r_e(x) analítico: nariz, cilindro y cola con las funciones de forma de OpenRocket | §3.1: perfil del cuerpo (siempre analítico en la optimización) |
| `sensor_lastre/perfiles.py` | `radio_superficie(x, g) -> ndarray` | r_e prolongado detrás de la base | §3.2: r_LE, r_TE,raíz y cierre P3 → P0 |
| `sensor_lastre/perfiles.py` | `radio_cola(u, forma, R, Ra, Lt, parametro, recortada)` | Radio de la transición | §8 T7: comparación con OpenRocket |
| `sensor_lastre/geometria.py` | `perfil_desde_muestras(x, x_s, r_s)` | Interpola un perfil de OpenRocket a la malla | §8 T7 |

## Erosión, volúmenes y momentos

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/geometria.py` | `construir_cavidad(caso, r_e=None, x=None) -> Cavidad` | Erosión multicapa por estación y r_i(x) | §4.2: una vez por cuerpo (caché) |
| `sensor_lastre/geometria.py` | `Cavidad.vol(a, b)`, `Cavidad.mom(a, b)`, `C0`, `C1` | ∀ y 𝓜 acumulados (Simpson) | §3.5 a través de `ModeloLastre` |

## Masas

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/masas.py` | `masa_vacia(caso, cav) -> MasaVacia` | Capas, mamparos, aletas y masas puntuales | §3.4: una vez por cuerpo; por candidato se reemplazan `m_aletas` y `x_aletas` con `dataclasses.replace` |
| `sensor_lastre/masas.py` | `masas_capas(caso, cav) -> list[MasaCapa]` | Masa y momento por capa y estación | §8 T3 |
| `sensor_lastre/aletas.py` | `GeomAleta.masa`, `GeomAleta.x_cg`, `area_centroide(P)` | n ρ φ t A por shoelace | §3.2: masa y centroide de la aleta (ρ = Onyx) |
| `sensor_lastre/aletas.py` | `poligono_global(puntos, x_LE, r_LE, r_sup, n_sup, tol)`, `es_simple(P)` | Polígono global con el cierre por la superficie; autointersección | §3.2: validación de la aleta contenida |

## Llenado de lastre, CG, SM, ventana y trim

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/estabilidad.py` | `x_inicio_lastre(caso, cav)`, `limite_geometrico(caso, cav, x_b0) -> LimiteGeo` | Inicio del lastre y ℓ_geo (fin de cavidad, cola, electrónica detrás, fracción de L) | §3.7.1: electrónica cabiendo según la lógica actual; una vez por cuerpo |
| `sensor_lastre/estabilidad.py` | `ModeloLastre`, `ventana_SM`, `alpha_trim` | CG(ℓ), ventana de SM, tapón trasero, trim | §3.5 a través de `_analizar_relleno` |
| `sensor_lastre/analisis.py` | `ResultadoCaso(...)`, `_analizar_relleno(res, rel) -> ResultadoRelleno`, `ResultadoRelleno.ell_delantero` | **Llenado actual**: tapón delantero desde la nariz hasta ℓ_SM,max; si llega a ℓ_geo, tapón trasero detrás de la electrónica hasta SM_min o cavidad llena; si la masa total supera `lastre.masa_max_sensor_g`, recorta el lastre (delantero primero). Trim, x_T y tolerancia de amarre | §3.5: se usa tal cual (`sensor_opt.lastre.llenar`), con el CP del sustituto inyectado como `x_CP`, `CNa` de `ResultadoCaso` y `SM_max = None` para que no retrase el lastre |
| `sensor_lastre/analisis.py` | `banderas_envolvente(caso)`, `tolerancia_amarre(...)` | Dimensiones con aletas; tolerancia de amarre | §3.3, §7 (dentro de `_analizar_relleno`) |

`_analizar_relleno` es privada en `sensor_lastre`; `sensor_opt` depende de su firma. Un cambio
de firma allí rompe `tests_opt/` (no los tests viejos), que es la alarma buscada.

## Barrowman interno

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/barrowman.py` | `cuerpo_revolucion(x, r, x0, x1, A_ref, nombre) -> Contribucion \| None` | C_Nα = 2 (A_a − A_f)/A_ref y x_CP de la forma general | §3.6: nariz y cola, una vez por cuerpo |
| `sensor_lastre/aletas.py` | `barrowman_freeform(g, A_ref, mach, n_franjas) -> CPAleta` | Aletas freeform por franjas (como FinSetCalc), K_TB con r_LE | §3.6: aletas de cada candidato; el span empieza en `GeomAleta.R_a`, que `sensor_opt` fija en r_TE,raíz |
| `sensor_lastre/barrowman.py` | `ResultadoCP`, `Contribucion` | Contenedores | §3.6 |

`cp_interno` no se usa directamente: suma y divide sin guarda. `sensor_opt.sustituto.cp_sustituto`
suma las mismas contribuciones y aplica la guarda Σ C_Nα ≥ ε_CN (v3 §3.6).

## OpenRocket

| Ruta | Firma | Qué hace | Uso en v3 |
|---|---|---|---|
| `sensor_lastre/or_bridge.py` | `PuenteOR(inst, nombres).cargar(ruta)` | Abre el `.ork` y ubica nariz, cuerpo, cola y aletas | §4.6 |
| `sensor_lastre/or_bridge.py` | `PuenteOR.aplicar(caso) -> dict` | Nariz, cuerpo, cola, paredes (ρ_eq) y una aleta freeform con extensión | §4.6: fija el cuerpo; la aleta contenida la escribe la extensión `sensor_opt.verificacion.PuenteOpt` |
| `sensor_lastre/or_bridge.py` | `PuenteOR.aero(mach, aoa, delta_aoa) -> ResultadoAero` | CP, C_Nα total y por componente, advertencias | §4.6 |
| `sensor_lastre/or_bridge.py` | `PuenteOR.masas() -> ResultadoMasas`, `perfil(n)`, `orl.save_doc` | Masas y CG por componente; perfil muestreado; guardado | §4.6, §7, T7, T10 |

`PuenteOpt(PuenteOR)` (en `sensor_opt/verificacion.py`) añade lo que el puente no cubre: los
cuatro puntos de la aleta contenida, el material Onyx y el lastre y la electrónica como
*Mass components*. Todo acceso a Java de `sensor_opt` vive en esa clase.
