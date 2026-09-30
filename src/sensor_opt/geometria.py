"""Cuerpo del candidato: perfil de la cola con k y L_t, cavidad, masas del casco y CP del cuerpo
(v3 §3.1, §4.2). Todo lo que no depende de las aletas se calcula aquí una vez por cuerpo.

El perfil es el analítico de `sensor_lastre.perfiles` (r_e = R_a + (R − R_a) f(ξ), R_a = k R),
la cavidad y las masas salen de `sensor_lastre.geometria` y `sensor_lastre.masas`, y la región
de lastre (x_b0, ℓ_geo, electrónica) de `sensor_lastre.estabilidad`, sin cambios.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from sensor_lastre.barrowman import Contribucion, cuerpo_revolucion
from sensor_lastre.config import Caso, ConfigError
from sensor_lastre.config import cargar as cargar_lastre
from sensor_lastre.estabilidad import LimiteGeo, limite_geometrico, x_inicio_lastre
from sensor_lastre.geometria import Cavidad, construir_cavidad
from sensor_lastre.masas import MasaVacia, masa_vacia
from sensor_lastre.perfiles import radio_superficie

from .config import ConfigOpt, CuerpoSpec, raw_cuerpo

MM = 1e-3


def theta_eq(R: float, k: float, L_t: float) -> float:
    """Ángulo equivalente de la cola, arctan(R (1 − k) / L_t) [rad] (v3 §3.1, informativo)."""
    return math.atan2(R * (1.0 - k), L_t)


@dataclass
class Cuerpo:
    """Caché por cuerpo (v3 §4.2). `motivos` no vacío = cuerpo descartado."""

    spec: CuerpoSpec
    motivos: list[str] = field(default_factory=list)
    caso: Caso | None = None  # caso de sensor_lastre con una aleta marcador
    cav: Cavidad | None = None
    mv: MasaVacia | None = None  # casco + mamparos + masas puntuales (la aleta se reemplaza)
    x_b0: float = math.nan
    lim: LimiteGeo | None = None
    partes_cuerpo: tuple[Contribucion, ...] = ()  # nariz y cola (Barrowman general)
    theta_eq: float = math.nan
    L_c: float = math.nan

    @property
    def ok(self) -> bool:
        return not self.motivos

    @property
    def geom(self):
        return self.caso.geom

    @property
    def x_t0(self) -> float:
        return self.caso.geom.x_cola

    def r_sup(self) -> Callable[[np.ndarray], np.ndarray]:
        g = self.caso.geom
        return lambda x: radio_superficie(x, g)


def construir_cuerpo(cfg: ConfigOpt, spec: CuerpoSpec) -> Cuerpo:
    """Resuelve el cuerpo con `sensor_lastre` y precalcula todo lo independiente de las aletas."""
    cu = Cuerpo(spec=spec)
    R = spec.D / 2
    cu.theta_eq = theta_eq(R, spec.k, spec.L_t)
    rest = cfg.restricciones
    if spec.k < rest.k_min - 1e-12:
        cu.motivos.append("k_bajo_minimo")
    if rest.angulo_cola_max is not None and cu.theta_eq > rest.angulo_cola_max + 1e-12:
        cu.motivos.append("angulo_cola")
    try:
        caso = cargar_lastre(raw_cuerpo(cfg, spec)).casos[0]
    except ConfigError as e:
        texto = " ".join(e.errores)
        if "L_n + L_t" in texto:
            cu.motivos.append("L_c_no_positivo")
        elif "Σ t_k" in texto or "R_a" in texto:
            cu.motivos.append("pared_mayor_que_popa")
        else:
            cu.motivos.append("cuerpo_invalido")
        return cu
    cu.caso = caso
    g = caso.geom
    cu.L_c = g.x_cola - g.nariz.Ln
    if not cu.L_c > 0:
        cu.motivos.append("L_c_no_positivo")
        return cu
    cu.cav = construir_cavidad(caso)
    cu.mv = masa_vacia(caso, cu.cav)
    cu.x_b0 = x_inicio_lastre(caso, cu.cav)
    cu.lim = limite_geometrico(caso, cu.cav, cu.x_b0)
    if not cu.lim.ell_geo > 0:
        cu.motivos.append("electronica_no_cabe" if cu.lim.limitante == "electronica_detras" else "sin_region_lastre")
    A_ref = math.pi * g.D**2 / 4
    partes = (cuerpo_revolucion(cu.cav.x, cu.cav.r_e, 0.0, g.nariz.Ln, A_ref, "nariz"),
              cuerpo_revolucion(cu.cav.x, cu.cav.r_e, g.x_cola, g.L, A_ref, "cola"))
    cu.partes_cuerpo = tuple(p for p in partes if p is not None)
    return cu
