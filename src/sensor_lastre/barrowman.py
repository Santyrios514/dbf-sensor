"""CP interno por Barrowman (sección 3.10), respaldo y verificación cruzada de OpenRocket.

Cuerpos de revolución (nariz y cola) con la forma general de Barrowman que usa OpenRocket
(SymmetricComponentCalc), válida para cualquier perfil:

    C_Nα = 2 (A_a − A_f) / A_ref,        x_CP = x_0 + (ℓ A_a − ∀) / (A_a − A_f)

con A_f, A_a las áreas delantera y trasera del componente, ℓ su longitud y ∀ su volumen.
Reproduce k_n = 2/3 (cono), 1/3 (elipsoide), 2n/(2n+1) (potencia) y ≈ 0.466 (ogiva esbelta),
y para la cola cónica da x = L_t/3 · [1 + (1 − d_F/d_R)/(1 − (d_F/d_R)²)].
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import Geometria
from .geometria import integrar


@dataclass(frozen=True)
class Contribucion:
    nombre: str
    CNa: float
    x: float


@dataclass(frozen=True)
class ResultadoCP:
    x_CP: float
    CNa: float
    partes: tuple[Contribucion, ...]


def cuerpo_revolucion(x: np.ndarray, r: np.ndarray, x0: float, x1: float, A_ref: float,
                      nombre: str) -> Contribucion | None:
    m = (x >= x0) & (x <= x1)
    xs, rs = x[m], r[m]
    if xs.size < 3:
        return None
    A_f, A_a = np.pi * rs[0] ** 2, np.pi * rs[-1] ** 2
    if np.isclose(A_a, A_f, rtol=1e-12, atol=0):
        return None
    V = integrar(np.pi * rs**2, xs)
    ell = xs[-1] - xs[0]
    return Contribucion(nombre, 2.0 * (A_a - A_f) / A_ref, xs[0] + (ell * A_a - V) / (A_a - A_f))


def aletas(g: Geometria) -> Contribucion:
    a, R, D = g.aletas, g.R, g.D
    cr, ct, s, xs = a.cr, a.ct, a.s, a.xs
    lm = np.sqrt(s**2 + (xs + ct / 2 - cr / 2) ** 2)
    CNa = (1 + R / (s + R)) * (4 * a.n * (s / D) ** 2) / (1 + np.sqrt(1 + (2 * lm / (cr + ct)) ** 2))
    xf = a.x_r0 + xs * (cr + 2 * ct) / (3 * (cr + ct)) + (cr + ct - cr * ct / (cr + ct)) / 6
    return Contribucion("aletas", float(CNa), float(xf))


def cp_interno(g: Geometria, x: np.ndarray, r_e: np.ndarray) -> ResultadoCP:
    """x_CP y C_Nα (referidos a S_ref = π D²/4) a partir del perfil exterior en la malla x."""
    A_ref = np.pi * g.D**2 / 4
    partes = [
        cuerpo_revolucion(x, r_e, 0.0, g.nariz.Ln, A_ref, "nariz"),
        cuerpo_revolucion(x, r_e, g.x_cola, g.L, A_ref, "cola"),
        aletas(g),
    ]
    partes = tuple(p for p in partes if p is not None)
    CNa = sum(p.CNa for p in partes)
    return ResultadoCP(x_CP=sum(p.CNa * p.x for p in partes) / CNa, CNa=CNa, partes=partes)
