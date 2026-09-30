"""Aleta freeform contenida en la cola (`freeform_cola_contenida`, v3 §2.1 y §3.2).

Parametrización (x global desde la punta, L = longitud del cuerpo):

    x_LE = x_t0 + λ L_t,   c_r = μ (L − x_LE),   c_t = γ c_r,   x_s = σ (L − x_LE − c_t)
    h = r_tip − r_LE,      r_LE = r_e(x_LE),     r_TE,raíz = r_e(x_LE + c_r)

Puntos locales de OpenRocket (origen en el borde de ataque de la raíz):

    P0 = (0, 0), P1 = (x_s, h), P2 = (x_s + c_t, h), P3 = (c_r, r_TE,raíz − r_LE)

y el cierre P3 → P0 sigue la superficie de la cola. Por construcción la aleta no sale de la cola.
El polígono, el área, el centroide, la masa y el Barrowman por franjas son los de
`sensor_lastre.aletas`; aquí solo se arman los puntos y se valida.

`GeomAleta.R_a` guarda el radio del punto más bajo de la raíz (r_TE,raíz), que es donde
`barrowman_freeform` y OpenRocket (`FinSet.getSpan`) empiezan el span.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sensor_lastre.aletas import GeomAleta, ParamsAleta, es_simple, poligono_global
from sensor_lastre.atmosfera import isa

from .config import AletaSpec, ConfigOpt
from .geometria import Cuerpo

MM = 1e-3


@dataclass(frozen=True)
class DimAleta:
    """Dimensiones derivadas de un candidato (SI), válidas aunque la aleta se descarte."""

    x_LE: float
    c_r: float
    c_t: float
    x_s: float
    r_LE: float
    r_TE_raiz: float
    r_tip: float
    h: float
    delta_b: float


def dimensiones(cu: Cuerpo, a: AletaSpec, L: float) -> DimAleta:
    x_t0, Lt, R = cu.x_t0, cu.spec.L_t, cu.spec.D / 2
    x_LE = x_t0 + a.lambda_LE * Lt
    c_r = a.mu_cr * (L - x_LE)
    c_t = a.gamma_ct * c_r
    x_s = a.sigma_flecha * (L - x_LE - c_t)
    r = cu.r_sup()
    r_LE, r_TE = (float(v) for v in r(np.array([x_LE, x_LE + c_r])))
    r_tip = a.r_tip_rel_R * R
    return DimAleta(x_LE=x_LE, c_r=c_r, c_t=c_t, x_s=x_s, r_LE=r_LE, r_TE_raiz=r_TE, r_tip=r_tip,
                    h=r_tip - r_LE, delta_b=L - (x_LE + c_r))


def params_sensor_lastre(cfg: ConfigOpt, d: DimAleta) -> ParamsAleta:
    """ParamsAleta de sensor_lastre equivalente (sin extensión; la punta se fija con los puntos)."""
    p = cfg.aleta
    return ParamsAleta(n=p.n, rotacion=p.rotacion, c_r=d.c_r, x_s=d.x_s, h=d.h, r_tip_obj=None, e=0.0,
                       r_in=0.0, delta_b=d.delta_b, eps=0.0, t=p.t, material=p.material, rho=p.rho,
                       phi=p.phi)


def puntos_locales(d: DimAleta) -> np.ndarray:
    return np.array([(0.0, 0.0), (d.x_s, d.h), (d.x_s + d.c_t, d.h), (d.c_r, d.r_TE_raiz - d.r_LE)])


def construir(cfg: ConfigOpt, cu: Cuerpo, a: AletaSpec) -> tuple[GeomAleta | None, DimAleta, list[str]]:
    """(geometría, dimensiones, motivos). Si hay motivos, la geometría es None (sin excepción)."""
    rest = cfg.restricciones
    d = dimensiones(cu, a, cfg.L)
    motivos = []
    if not d.h >= rest.h_min - 1e-12:
        motivos.append("h_menor_minimo")
    if not d.c_r >= rest.c_min - 1e-12:
        motivos.append("c_r_menor_minimo")
    if not d.c_t >= rest.c_min - 1e-12:
        motivos.append("c_t_menor_minimo")
    if motivos:
        return None, d, motivos
    P = puntos_locales(d)
    num = cfg.numerico
    n_sup = int(num.get("n_superficie_aleta", 60))
    tol_sup = float(num.get("tol_superficie_mm", 0.01)) * MM
    y3_sup = float(cu.r_sup()(np.array([d.x_LE + d.c_r]))[0]) - d.r_LE
    if abs(P[3, 1] - y3_sup) >= tol_sup:
        motivos.append("P3_fuera_de_superficie")
    poli, orig = poligono_global(P, d.x_LE, d.r_LE, cu.r_sup(), n_sup, tol=1e-9)
    if not es_simple(poli):
        motivos.append("poligono_no_simple")
    if motivos:
        return None, d, motivos
    g = GeomAleta(params=params_sensor_lastre(cfg, d), x_LE=d.x_LE, r_LE=d.r_LE, R_a=d.r_TE_raiz, h=d.h,
                  puntos=P, poligono=poli, origen=tuple(orig))
    return g, d, []


# --------------------------------------------------------------------------- flutter (v3 §3.2.1)


def velocidad_flutter(h: float, A_f: float, c_r: float, c_t: float, t: float, E: float, nu: float,
                      altitud: float) -> float:
    """Martin (NACA TN 4197) sobre la aleta trapezoidal equivalente [m/s]:

        V_f = a sqrt(G / (1.337 AR³ P (λ + 1) / (2 (AR + 2) (t/c_r)³))),
        AR = h²/A_f, λ = c_t/c_r, G = E / (2 (1 + ν))

    G y P en las mismas unidades (el cociente es adimensional); P y a de la ISA.
    """
    if not (h > 0 and A_f > 0 and c_r > 0 and t > 0 and E > 0):
        return math.nan
    atm = isa(altitud)
    AR = h**2 / A_f
    lam = c_t / c_r
    Gmod = E / (2 * (1 + nu))
    den = 1.337 * AR**3 * atm.p * (lam + 1) / (2 * (AR + 2) * (t / c_r) ** 3)
    return atm.a * math.sqrt(Gmod / den)
