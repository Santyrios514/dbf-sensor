"""Masas estructurales y sus primeros momentos (sección 3.5)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .aletas import GeomAleta
from .config import ESTACIONES, Caso
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
    m_puntuales: float = 0.0  # todas, también las que van en el CG
    M_puntuales: float = 0.0  # solo las de posición fija
    m_en_CG: float = 0.0  # masas puntuales ubicadas en el CG del sensor (no mueven el CG)
    V_aloj: float = 0.0  # volumen de lastre desplazado por sus alojamientos (también en el CG)

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


def masa_aletas(a: GeomAleta) -> tuple[float, float]:
    """(m_f, x̄_f) de las n aletas freeform: m_f = n ρ_f φ_f t_f A_f (v2 §3.4)."""
    return a.masa, a.x_cg


def masas_por_componente(caso: Caso, mv: "MasaVacia") -> dict[str, tuple[float, float]]:
    """{componente: (m, x_CG)} con los mismos cortes que OpenRocket (nariz, cuerpo, cola, aletas)."""
    out = {}
    for s in ESTACIONES:
        cs = [c for c in mv.capas if c.estacion == s]
        m = sum(c.m for c in cs)
        out[s] = (m, sum(c.M for c in cs) / m if m > 0 else float("nan"))
    out["aletas"] = (mv.m_aletas, mv.x_aletas)
    return out


def masa_vacia(caso: Caso, cav: Cavidad) -> MasaVacia:
    mv = MasaVacia(capas=masas_capas(caso, cav))
    mv.m_aletas, mv.x_aletas = masa_aletas(caso.geom.aletas)
    for mp in caso.mamparos:
        mv.m_mamparos += mp.rho * float(cav.vol_bruto(mp.x, mp.x + mp.e))
        mv.M_mamparos += mp.rho * float(cav.mom_bruto(mp.x, mp.x + mp.e))
    for p in caso.puntuales:
        mv.m_puntuales += p.m
        if p.en_CG:
            mv.m_en_CG += p.m
            mv.V_aloj += p.V_aloj
        else:
            mv.M_puntuales += p.m * p.x
    return mv
