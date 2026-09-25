"""Perfiles exteriores analíticos r_e(x) (modo sin OpenRocket), secciones 3.2 (v1) y 3.1 (v2).

Las funciones de forma replican `Transition.Shape.getRadius(x, radius, length, param)` de
OpenRocket (sin recorte, `clipped = false`), donde x se mide desde el extremo de menor radio:

- nariz:      r(x) = shape(x, R, L_n, p)
- transición: r(x) = R_a + shape(L_t − u, R_f − R_a, L_t, p), con u medido desde la unión

x global se mide desde la punta de la nariz, positivo hacia popa. Todo en metros.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .config import Geometria

MM = 1e-3

# forma del YAML -> (forma de OpenRocket, parámetro por defecto)
FORMAS = {
    "conica": ("CONICAL", None),
    "ogiva": ("OGIVE", 1.0),
    "ogiva_tangente": ("OGIVE", 1.0),
    "elipsoide": ("ELLIPSOID", None),
    "potencia": ("POWER", None),
    "parabolica": ("PARABOLIC", 1.0),
    "haack": ("HAACK", 0.0),
}


def forma_or(x, forma: str, radius: float, length: float, param: float | None) -> np.ndarray:
    """Radio de la forma `forma` a la distancia x del extremo fino (convención de OpenRocket)."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, length)
    xi = x / length
    if forma == "conica":
        return radius * xi
    if forma in ("ogiva", "ogiva_tangente"):
        p = 1.0 if param is None else float(param)
        if p < 0.001:
            return radius * xi
        Rc = np.sqrt((length**2 + radius**2) * (((2 - p) * length) ** 2 + (p * radius) ** 2)
                     / (4 * (p * radius) ** 2))
        Lc = length / p
        y0 = np.sqrt(max(Rc**2 - Lc**2, 0.0))
        return np.sqrt(np.clip(Rc**2 - (Lc - x) ** 2, 0.0, None)) - y0
    if forma == "elipsoide":
        return radius * np.sqrt(np.clip(2 * xi - xi**2, 0.0, None))
    if forma == "potencia":
        return radius * xi ** float(param)
    if forma == "parabolica":
        K = 1.0 if param is None else float(param)
        return radius * (2 * xi - K * xi**2) / (2 - K)
    if forma == "haack":
        C = 0.0 if param is None else float(param)
        th = np.arccos(1.0 - 2.0 * xi)
        return radius * np.sqrt(np.clip((th - np.sin(2 * th) / 2 + C * np.sin(th) ** 3) / np.pi, 0.0, None))
    raise ValueError(f"forma desconocida: {forma}")


def radio_nariz(x, forma: str, parametro: float | None, R: float, Ln: float) -> np.ndarray:
    return forma_or(x, forma, R, Ln, parametro)


def radio_cola(u, forma: str, R: float, Ra: float, Lt: float, parametro: float | None = None) -> np.ndarray:
    """Radio de la cola en función de u = x − (L − L_t) ∈ [0, L_t] (se estrecha hacia popa)."""
    u = np.clip(np.asarray(u, dtype=float), 0.0, Lt)
    if Lt == 0 or R == Ra:
        return np.full_like(u, R)
    return Ra + forma_or(Lt - u, forma, R - Ra, Lt, parametro)


def radio_exterior(x, g: "Geometria") -> np.ndarray:
    """r_e(x) del cuerpo en [0, L]; 0 fuera de ese intervalo."""
    x = np.asarray(x, dtype=float)
    r = np.full_like(x, g.R)
    m_n = x < g.nariz.Ln
    r[m_n] = radio_nariz(x[m_n], g.nariz.forma, g.nariz.parametro, g.R, g.nariz.Ln)
    m_t = x > g.x_cola
    r[m_t] = radio_cola(x[m_t] - g.x_cola, g.cola.forma, g.R, g.cola.Ra, g.cola.Lt, g.cola.parametro)
    r[(x < 0) | (x > g.L)] = 0.0
    return r


def radio_superficie(x, g: "Geometria") -> np.ndarray:
    """Como `radio_exterior`, pero prolongando el último radio más allá de la base (la raíz de
    las aletas puede terminar unas centésimas de mm detrás de la popa por el offset *bottom*)."""
    return radio_exterior(np.clip(np.asarray(x, dtype=float), 0.0, g.L), g)


def muestrear(g: "Geometria", n_por_componente: int, config_id: str = "") -> pd.DataFrame:
    """Perfil en formato largo (esquema de or_perfiles.csv): config_id, componente, x_mm, r_ext_mm."""
    tramos = [("nariz", 0.0, g.nariz.Ln), ("cuerpo", g.nariz.Ln, g.x_cola), ("cola", g.x_cola, g.L)]
    filas = []
    for comp, a, b in tramos:
        if b <= a:
            continue
        xs = np.linspace(a, b, n_por_componente)
        if comp == "nariz":
            r = radio_nariz(xs, g.nariz.forma, g.nariz.parametro, g.R, g.nariz.Ln)
        elif comp == "cola":
            r = radio_cola(xs - g.x_cola, g.cola.forma, g.R, g.cola.Ra, g.cola.Lt, g.cola.parametro)
        else:
            r = np.full_like(xs, g.R)
        filas.append(pd.DataFrame({"config_id": config_id, "componente": comp,
                                   "x_mm": xs / MM, "r_ext_mm": r / MM}))
    return pd.concat(filas, ignore_index=True)
