"""Atmósfera estándar internacional (ISA), troposfera (0–11 km)."""

from __future__ import annotations

import math
from dataclasses import dataclass

# Constantes ISA (ISO 2533)
T0 = 288.15  # K
P0 = 101325.0  # Pa
LAPSE = 0.0065  # K/m
R_AIRE = 287.05287  # J/(kg K)
GAMMA = 1.4
G_ISA = 9.80665  # m/s^2
H_TROPOPAUSA = 11000.0  # m


@dataclass(frozen=True)
class EstadoISA:
    T: float  # K
    p: float  # Pa
    rho: float  # kg/m^3
    a: float  # m/s


def isa(h: float) -> EstadoISA:
    """Estado ISA a la altitud geopotencial h [m] (troposfera)."""
    if not (-500.0 <= h <= H_TROPOPAUSA):
        raise ValueError(f"altitud {h} m fuera del rango de la troposfera ISA")
    T = T0 - LAPSE * h
    p = P0 * (T / T0) ** (G_ISA / (R_AIRE * LAPSE))
    rho = p / (R_AIRE * T)
    a = math.sqrt(GAMMA * R_AIRE * T)
    return EstadoISA(T=T, p=p, rho=rho, a=a)
