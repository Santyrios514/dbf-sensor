# API de OpenRocket verificada (OpenRocket 24.12, orlab 0.10.0)

Firmas comprobadas con `javap` sobre `OpenRocket-24.12.jar` y ejecutadas por JPype desde
`src/sensor_lastre/or_bridge.py`, que es el único módulo que toca Java. Paquete raíz:
`info.openrocket.core` en 24.12 y `net.sf.openrocket` en ≤ 23.09. orlab lo resuelve y lo expone
como `inst.openrocket`.

## orlab

| Uso | Llamada |
|---|---|
| Arrancar la JVM | `with orlab.OpenRocketInstance(jar_path=None) as inst:` (`None` → `ORLAB_JAR` o el caché de `python -m orlab fetch`) |
| Helper | `orl = orlab.Helper(inst)`; `orl.load_doc(ruta)` → `OpenRocketDocument`; `orl.save_doc(ruta, doc)` |
| Paquete Java | `core = inst.openrocket` (JPackage `info.openrocket.core`) |
| Versión | `core.util.BuildProperties.getVersion()` → `"24.12"`; orlab: `importlib.metadata.version("orlab")` |

## Árbol del cohete

- `doc.getRocket()` → `Rocket`.
- `rocket.iterator(True)` recorre todos los componentes; `getName()` devuelve el nombre y
  `getClass().getSimpleName()` el tipo.
- `rocket.getSelectedConfiguration()` → `FlightConfiguration`. De ella salen
  `getReferenceLength()` (el D de referencia con `referencetype = maximum`) y `getLength()`
  (la longitud total **con** aletas: 390 mm en la base).
- `RocketComponent.getComponentLocations()[0].x` da la posición axial absoluta del inicio del componente.

## Componentes y setters

| Clase | Métodos usados |
|---|---|
| `NoseCone` (hereda de `Transition`) | `setLength(double)`, `setBaseRadius(double)`, `getBaseRadius()`, `setShapeType(Transition.Shape)`, `setShapeParameter(double)`, `getRadius(double x_local)` |
| `BodyTube` | `setLength(double)`, `setOuterRadiusAutomatic(boolean)`, `getOuterRadius()`, `getRadius(double)` |
| `Transition` | `setLength`, `setForeRadiusAutomatic(boolean)`, `setAftRadius(double)`, `setShapeType`, `setShapeParameter`, `setClipped(boolean)`, `isClipped()`, `getRadius(double x_local)` (fuera de [0, L] devuelve el radio extremo) |
| `SymmetricComponent` | `setThickness(double)`, `getThickness()`, `getComponentMass()`, `getComponentCG()` (x local) |
| `ExternalComponent` | `setMaterial(Material)`, `getMaterial().getDensity()` |
| `FreeformFinSet` | `setPoints(Coordinate[])`, `getFinPoints()`, `getFinPointsWithRoot()` (puntos + raíz sobre la superficie; repite P5 y P0), `getBodyRadius()` (= r_LE), `getSpan()` (desde el punto más bajo de la raíz hasta la punta), `getPlanformArea()` |
| `FinSet` | `setFinCount(int)`, `setBaseRotation(double rad)`, `setThickness(double)`, `getInstanceAngles()` |
| `RocketComponent` | `setAxialMethod(AxialMethod.BOTTOM)`, `setAxialOffset(double)` |
| `Material` | `Material.newMaterial(Material.Type.BULK, nombre, densidad, userDefined=True)` |
| `Coordinate` | `Coordinate(x, y)`; campos `x, y, z, weight` |

`Transition.Shape`: `CONICAL`, `OGIVE`, `ELLIPSOID`, `POWER`, `PARABOLIC`, `HAACK`. Se verificó
`isClippable()` y `usesParameter()` para cada forma:

| Forma | Recortable | Usa parámetro |
|---|---|---|
| CONICAL | no | no |
| OGIVE | no | sí |
| ELLIPSOID | **sí** | no |
| POWER | **sí** | sí |
| PARABOLIC | no | sí |
| HAACK | **sí** | sí |

Una `Transition` nueva viene con `clipped = true`. El XML solo guarda `<shapeclipped>` si la forma
es recortable; por eso la transición parabólica del `.ork` no lo lleva. Con el recorte,
`r(x) = shape(c + s, R_f, c + L_t)`, donde s es la distancia al extremo fino y c cumple
`shape(c, R_f, c + L_t) = R_a` (`calculateClip`). La implementación analítica está en
`perfiles.radio_cola(..., recortada=True)` y coincide con OpenRocket en menos de 0.05 mm en las
seis formas, con y sin recorte (tests/test_openrocket.py).

**Puntos inválidos:** `FreeformFinSet.setPoints` **no lanza excepción**. Ante un punto dentro del
cuerpo, OpenRocket lo pega a la superficie en silencio. `or_bridge.aplicar()` compara los puntos
aceptados con los pedidos y lanza `ErrorAleta` (→ `estado = error_aleta`) si difieren.

## Aerodinámica

- `calc = core.aerodynamics.BarrowmanCalculator()`.
- `cond = core.aerodynamics.FlightConditions(fc)`; `cond.setMach(double)`; `cond.setAOA(double rad)`.
- `ws = core.logging.WarningSet()` es iterable; `str(w)` da el texto de cada advertencia.
- `calc.getCP(fc, cond, ws)` → `Coordinate`: `x` es el CP total y **`weight` es el C_Nα total**
  (verificado: 3.358 en la base a M = 0.3).
- `calc.getAerodynamicForces(fc, cond, ws)` → `AerodynamicForces`: `getCN()` (a la AoA dada) y
  `getCD()`. El C_Nα por diferencias centradas se obtiene con `getCN()` a AoA ± Δ.
- `calc.getForceAnalysis(fc, cond, ws)` → `Map<RocketComponent, AerodynamicForces>`. Por
  componente, `getCP().weight` es su C_Nα y `getCP().x` su CP. El mapa incluye además la etapa y
  el cohete, que se descartan.

En la base, a M = 0.3: nariz 2.000 @ 21.67 mm; transición −1.280 @ 326.17 mm; aletas 2.638 @
342.41 mm; total 3.358 @ 157.56 mm.

## Masas

- `core.masscalc.MassCalculator.calculateStructure(fc)` → `RigidBody`, con `getMass()` y `getCM().x`.
- Por componente: `getComponentMass()` (incluye todas las instancias de aletas) y
  `getComponentCG().x`, que es local y se suma a `getComponentLocations()[0].x`.

En la base: nariz 26.014 g, cuerpo 105.515 g, transición 27.488 g, aletas 25.317 g; total 184.335 g
con CG en 202.02 mm.

## Método de aletas de OpenRocket (para el Barrowman interno)

`FinSetCalc` corta las cuerdas por franjas en y **solo contra la polilínea de puntos**
P0 → … → Pk, sin el cierre por la superficie de la cola. El span va del punto más bajo de la raíz
a la punta (31.24 mm en la base, no los 20 mm de h). Con esa convención, el método interno
(`aletas.barrowman_freeform`) reproduce el CP de las aletas de OpenRocket en menos de 0.6 mm y su
C_Nα en menos de 7 % en todo el barrido.

## Extensiones de `sensor_opt` (`PuenteOpt`, spec v3)

Verificado con OpenRocket 24.12 desde `src/sensor_opt/verificacion.py`:

| Uso | Llamada |
|---|---|
| Masa concentrada | `MassComponent(length, radius, mass)`; `setName(str)`; `setComponentMass(double)` |
| Colgarla del cuerpo | `BodyTube.addChild(mc)` |
| Posición absoluta | `mc.setAxialMethod(AxialMethod.ABSOLUTE)`; `mc.setAxialOffset(x)` (x = borde delantero desde la punta; puede quedar fuera del tramo del padre) |
| CG propio | `mc.getComponentCG().x` = `length / 2` (local) |
| Masa total | `MassCalculator.calculateStructure` **incluye** los `MassComponent` (comprobado: 1 kg en x = 30.5 mm da el CG esperado a 1e-15) |
| Aleta contenida | `FreeformFinSet.setAxialMethod(AxialMethod.BOTTOM)`; `setAxialOffset(−δ_b)`: el offset positivo mueve el borde de salida de la raíz **hacia popa** de la base, así que una aleta que termina δ_b antes de la base lleva offset negativo |
| Radio de la cola en x local | `Transition.getRadius(x)` para r_LE y r_TE,raíz (P3 cae exactamente sobre la superficie de OpenRocket) |

Los `MassComponent` y el material `ONYX` se guardan en el `.ork` con `orl.save_doc` y al reabrirlo
reproducen masa, CG y CP (T10).
