"""Densidades de materiales, densidad efectiva de rellenos y densidad equivalente de paredes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class Capa:
    """Capa de pared: material, densidad [kg/m³], espesor [m] y fracción sólida (infill)."""

    material: str
    rho: float
    t: float
    phi: float = 1.0


@dataclass(frozen=True)
class Relleno:
    """Material de lastre con su densidad efectiva ρ_b [kg/m³]."""

    nombre: str
    rho_b: float
    costo_usd_kg: float = float("nan")  # costo aproximado por kg de relleno (material, sin mecanizado)


def densidad_granular(rho_grano: float, rho_matriz: float, empaquetamiento: float) -> float:
    """ρ_b = φ ρ_grano + (1 − φ) ρ_matriz."""
    return empaquetamiento * rho_grano + (1.0 - empaquetamiento) * rho_matriz


def costo_granular(phi: float, rho_grano: float, rho_matriz: float, c_grano: float, c_matriz: float) -> float:
    """Costo por kg de un relleno granular: promedio de los costos ponderado por fracción de masa,
    w_grano = φ ρ_grano / ρ_b."""
    rho_b = densidad_granular(rho_grano, rho_matriz, phi)
    w = phi * rho_grano / rho_b
    return w * c_grano + (1 - w) * c_matriz


def densidad_equivalente(capas: Sequence[Capa]) -> float:
    """ρ_eq = Σ ρ_k φ_k t_k / Σ t_k (una sola capa equivalente para OpenRocket)."""
    t_total = sum(c.t for c in capas)
    return sum(c.rho * c.phi * c.t for c in capas) / t_total


def espesor_total(capas: Sequence[Capa]) -> float:
    return sum(c.t for c in capas)
