"""CP sustituto (Barrowman interno) con guarda numérica y calibración lineal contra OpenRocket
(v3 §3.6 y §4.7).

    x_CP = Σ C_Nα,i x_i / Σ C_Nα,i,   con C_Nα,nariz + C_Nα,cola = 2 k²

El CP no diverge por cerrar la cola: con k = 0 la cola aporta −2 y el CP sigue siendo finito
mientras haya aletas. La singularidad aparece cuando 2k² + C_Nα,f → 0 (cola cerrada y aletas
muy pequeñas); por eso se descarta el candidato si Σ C_Nα < ε_CN, antes de dividir.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sensor_lastre.aletas import GeomAleta, barrowman_freeform
from sensor_lastre.barrowman import Contribucion, ResultadoCP

from .geometria import Cuerpo

MOTIVO_CN = "CN_total_pequeno"


def combinar(partes: tuple[Contribucion, ...], eps_CN: float) -> tuple[ResultadoCP | None, str | None]:
    """Suma de contribuciones con la guarda Σ C_Nα ≥ ε_CN. Nunca devuelve inf ni NaN."""
    CNa = sum(p.CNa for p in partes)
    if not (math.isfinite(CNa) and CNa >= eps_CN):
        return None, MOTIVO_CN
    x = sum(p.CNa * p.x for p in partes) / CNa
    if not math.isfinite(x):
        return None, MOTIVO_CN
    return ResultadoCP(x_CP=x, CNa=CNa, partes=partes), None


def cp_sustituto(cu: Cuerpo, g: GeomAleta, mach: float, n_franjas: int,
                 eps_CN: float) -> tuple[ResultadoCP | None, str | None]:
    """CP del candidato: nariz y cola del caché del cuerpo + aletas por franjas (sensor_lastre)."""
    A_ref = math.pi * cu.spec.D**2 / 4
    f = barrowman_freeform(g, A_ref, mach, n_franjas)
    return combinar(cu.partes_cuerpo + (Contribucion("aletas", f.CNa, f.x_CP),), eps_CN)


# --------------------------------------------------------------------------- calibración (v3 §4.7)


@dataclass(frozen=True)
class Calibracion:
    """x_CP,cal = α + β x_CP,sust (m)."""

    alpha: float = 0.0
    beta: float = 1.0
    R2: float = math.nan
    res_max: float = math.nan  # m, máx |x_OR − x_cal| sobre la muestra
    n: int = 0
    iteracion: int = 0

    def aplicar(self, x):
        return self.alpha + self.beta * np.asarray(x, dtype=float)

    @property
    def identidad(self) -> bool:
        return self.alpha == 0.0 and self.beta == 1.0

    def a_dict(self) -> dict:
        return {"alpha_mm": self.alpha * 1e3, "beta": self.beta, "R2": self.R2,
                "residuo_max_mm": self.res_max * 1e3, "n": self.n, "iteracion": self.iteracion}

    @classmethod
    def desde_dict(cls, d: dict) -> "Calibracion":
        return cls(alpha=d["alpha_mm"] * 1e-3, beta=d["beta"], R2=d.get("R2", math.nan),
                   res_max=d.get("residuo_max_mm", math.nan) * 1e-3, n=d.get("n", 0),
                   iteracion=d.get("iteracion", 0))


IDENTIDAD = Calibracion()


def ajustar(x_sust, x_or, iteracion: int = 0) -> Calibracion:
    """Mínimos cuadrados de x_OR ≈ α + β x_sust."""
    x = np.asarray(x_sust, dtype=float)
    y = np.asarray(x_or, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 2 or np.ptp(x) == 0:
        raise ValueError("se necesitan al menos dos x_CP sustitutos distintos para calibrar")
    A = np.column_stack([np.ones_like(x), x])
    (alpha, beta), *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - (alpha + beta * x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    R2 = 1.0 - float((r**2).sum()) / ss_tot if ss_tot > 0 else 1.0
    return Calibracion(alpha=float(alpha), beta=float(beta), R2=R2, res_max=float(np.abs(r).max()),
                       n=int(x.size), iteracion=iteracion)
