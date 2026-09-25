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


def densidad_granular(rho_grano: float, rho_matriz: float, empaquetamiento: float) -> float:
    """ρ_b = φ ρ_grano + (1 − φ) ρ_matriz."""
    return empaquetamiento * rho_grano + (1.0 - empaquetamiento) * rho_matriz


def densidad_equivalente(capas: Sequence[Capa]) -> float:
    """ρ_eq = Σ ρ_k φ_k t_k / Σ t_k (una sola capa equivalente para OpenRocket)."""
    t_total = sum(c.t for c in capas)
    return sum(c.rho * c.phi * c.t for c in capas) / t_total


def espesor_total(capas: Sequence[Capa]) -> float:
    return sum(c.t for c in capas)
