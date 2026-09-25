"""Masas estructurales y sus primeros momentos (sección 3.5)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import ESTACIONES, Aletas, Caso
from .geometria import Cavidad, integrar


@dataclass(frozen=True)
class MasaCapa:
    estacion: str
    capa_idx: int  # 1 = exterior
    material: str
    espesor: float
    phi: float
    m: float
    M: float  # primer momento m·x̄

    @property
    def x_cg(self) -> float:
        return self.M / self.m if self.m > 0 else float("nan")


@dataclass
class MasaVacia:
    capas: list[MasaCapa] = field(default_factory=list)
    m_aletas: float = 0.0
    x_aletas: float = 0.0
    m_mamparos: float = 0.0
    M_mamparos: float = 0.0
    m_puntuales: float = 0.0
    M_puntuales: float = 0.0

    @property
    def m_casco(self) -> float:
        return sum(c.m for c in self.capas)

    @property
    def M_casco(self) -> float:
        return sum(c.M for c in self.capas)

    @property
    def m0(self) -> float:
        return self.m_casco + self.m_mamparos + self.m_aletas + self.m_puntuales

    @property
    def M0(self) -> float:
        return self.M_casco + self.M_mamparos + self.m_aletas * self.x_aletas + self.M_puntuales


def masas_capas(caso: Caso, cav: Cavidad) -> list[MasaCapa]:
    """m_k = ρ_k φ_k π ∫ (r_{k−1}² − r_k²) dx por estación, sin aproximación de pared delgada."""
    out = []
    x = cav.x
    for i, s in enumerate(ESTACIONES):
        mask = cav.estacion == i
        if not mask.any():
            continue
        fr = cav.fronteras[s]
        for k, capa in enumerate(caso.pared[s], start=1):
            dA = np.pi * (fr[k - 1] ** 2 - fr[k] ** 2) * mask
            f = capa.rho * capa.phi
            out.append(MasaCapa(
                estacion=s, capa_idx=k, material=capa.material, espesor=capa.t, phi=capa.phi,
                m=f * integrar(dA, x), M=f * integrar(dA * x, x),
            ))
    return out


def masa_aletas(a: Aletas) -> tuple[float, float]:
    """(m_f, x̄_f) de n aletas trapezoidales."""
    m = a.n * a.rho * a.phi * a.tf * a.s * (a.cr + a.ct) / 2.0
    return m, centroide_aleta(a.cr, a.ct, a.xs, a.x_r0)


def centroide_aleta(cr: float, ct: float, xs: float, x_r0: float = 0.0) -> float:
    """x̄ del área de una aleta trapezoidal (flecha x_s del borde de ataque)."""
    return x_r0 + (cr**2 + cr * ct + ct**2 + xs * (cr + 2 * ct)) / (3.0 * (cr + ct))


def masa_vacia(caso: Caso, cav: Cavidad) -> MasaVacia:
    mv = MasaVacia(capas=masas_capas(caso, cav))
    mv.m_aletas, mv.x_aletas = masa_aletas(caso.geom.aletas)
    for mp in caso.mamparos:
        mv.m_mamparos += mp.rho * float(cav.vol_bruto(mp.x, mp.x + mp.e))
        mv.M_mamparos += mp.rho * float(cav.mom_bruto(mp.x, mp.x + mp.e))
    for p in caso.puntuales:
        mv.m_puntuales += p.m
        mv.M_puntuales += p.m * p.x
    return mv
