"""Perfiles exteriores analíticos r_e(x) (modo sin OpenRocket), sección 3.2.

x se mide desde la punta de la nariz, positivo hacia popa. Todo en metros.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Geometria

MM = 1e-3


def radio_nariz(x: np.ndarray, forma: str, parametro: float | None, R: float, Ln: float) -> np.ndarray:
    xi = np.clip(np.asarray(x, dtype=float) / Ln, 0.0, 1.0)
    if forma == "conica":
        return R * xi
    if forma == "elipsoide":
        return R * np.sqrt(np.clip(1.0 - (1.0 - xi) ** 2, 0.0, None))
    if forma == "ogiva_tangente":
        rho_o = (R**2 + Ln**2) / (2.0 * R)
        xx = xi * Ln
        return np.sqrt(np.clip(rho_o**2 - (Ln - xx) ** 2, 0.0, None)) + R - rho_o
    if forma == "potencia":
        return R * xi ** float(parametro)
    if forma == "haack":
        C = 0.0 if parametro is None else float(parametro)
        th = np.arccos(1.0 - 2.0 * xi)
        return R / np.sqrt(np.pi) * np.sqrt(np.clip(th - np.sin(2 * th) / 2 + C * np.sin(th) ** 3, 0.0, None))
    raise ValueError(f"forma de nariz desconocida: {forma}")


def radio_cola(u: np.ndarray, forma: str, R: float, Ra: float, Lt: float) -> np.ndarray:
    """Radio de la cola en función de u = x − (L − L_t) ∈ [0, L_t]; tangente al cuerpo en u = 0
    para `ogiva` y `parabolica`."""
    u = np.clip(np.asarray(u, dtype=float), 0.0, Lt)
    if Lt == 0:
        return np.full_like(u, R)
    delta = R - Ra
    if forma == "conica":
        return R + (Ra - R) * u / Lt
    if delta == 0:
        return np.full_like(u, R)
    if forma == "parabolica":
        return R - delta * (u / Lt) ** 2
    if forma == "ogiva":
        rho_o = (delta**2 + Lt**2) / (2.0 * delta)
        return R - (rho_o - np.sqrt(np.clip(rho_o**2 - u**2, 0.0, None)))
    raise ValueError(f"forma de cola desconocida: {forma}")


def radio_exterior(x: np.ndarray, g: Geometria) -> np.ndarray:
    """r_e(x) en [0, L]; 0 fuera de ese intervalo."""
    x = np.asarray(x, dtype=float)
    r = np.full_like(x, g.R)
    m_n = x < g.nariz.Ln
    r[m_n] = radio_nariz(x[m_n], g.nariz.forma, g.nariz.parametro, g.R, g.nariz.Ln)
    m_t = x > g.x_cola
    r[m_t] = radio_cola(x[m_t] - g.x_cola, g.cola.forma, g.R, g.cola.Ra, g.cola.Lt)
    r[(x < 0) | (x > g.L)] = 0.0
    return r


def muestrear(g: Geometria, n_por_componente: int, config_id: str = "") -> pd.DataFrame:
    """Perfil en formato largo (esquema de or_perfiles.csv): config_id, componente, x_mm, r_ext_mm."""
    tramos = [("nariz", 0.0, g.nariz.Ln), ("cuerpo", g.nariz.Ln, g.x_cola), ("cola", g.x_cola, g.L)]
    filas = []
    for comp, a, b in tramos:
        if b <= a:
            continue
        xs = np.linspace(a, b, n_por_componente)
        # evaluar dentro del tramo (el extremo derecho pertenece al tramo, no al siguiente)
        if comp == "nariz":
            r = radio_nariz(xs, g.nariz.forma, g.nariz.parametro, g.R, g.nariz.Ln)
        elif comp == "cola":
            r = radio_cola(xs - g.x_cola, g.cola.forma, g.R, g.cola.Ra, g.cola.Lt)
        else:
            r = np.full_like(xs, g.R)
        filas.append(pd.DataFrame({"config_id": config_id, "componente": comp,
                                   "x_mm": xs / MM, "r_ext_mm": r / MM}))
    return pd.concat(filas, ignore_index=True)
