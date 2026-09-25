"""CP interno por Barrowman (sección 3.10), respaldo y verificación cruzada de OpenRocket.

Cuerpos de revolución (nariz y cola) con la forma general de Barrowman que usa OpenRocket
(SymmetricComponentCalc; Niskanen ecs. 3.19 y 3.28), válida para cualquier perfil:

    C_Nα = 2 (A_a − A_f) / A_ref,        x_CP = x_0 + (ℓ A_a − ∀) / (A_a − A_f)

con A_f, A_a las áreas delantera y trasera del componente, ℓ su longitud y ∀ su volumen.
Reproduce k_n = 2/3 (cono), 1/3 (elipsoide), 2n/(2n+1) (potencia) y ≈ 0.466 (ogiva esbelta),
y para la cola cónica da x = L_t/3 · [1 + (1 − d_F/d_R)/(1 − (d_F/d_R)²)].
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import Geometria
from .aletas import barrowman_freeform
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


def cp_interno(g: Geometria, x: np.ndarray, r_e: np.ndarray, mach: float, n_franjas: int) -> ResultadoCP:
    """x_CP y C_Nα (referidos a S_ref = π D²/4) a partir del perfil exterior en la malla x.

    Cuerpo: forma general (v2 §3.2). Aletas freeform: franjas en y (v2 §3.7, aletas.py).
    """
    A_ref = np.pi * g.D**2 / 4
    partes = [
        cuerpo_revolucion(x, r_e, 0.0, g.nariz.Ln, A_ref, "nariz"),
        cuerpo_revolucion(x, r_e, g.x_cola, g.L, A_ref, "cola"),
    ]
    if g.aletas is not None:
        f = barrowman_freeform(g.aletas, A_ref, mach, n_franjas)
        partes.append(Contribucion("aletas", f.CNa, f.x_CP))
    partes = tuple(p for p in partes if p is not None)
    CNa = sum(p.CNa for p in partes)
    return ResultadoCP(x_CP=sum(p.CNa * p.x for p in partes) / CNa, CNa=CNa, partes=partes)
