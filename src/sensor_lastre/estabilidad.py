"""CG(ℓ), margen estático, ventana de longitudes de lastre y trim (secciones 3.6–3.9, 3.11)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import G0
from .config import Caso
from .geometria import Cavidad
from .masas import MasaVacia


# --------------------------------------------------------------------------- región de lastre


def x_inicio_lastre(caso: Caso, cav: Cavidad) -> float:
    """x_b0: absoluto, o el primer x donde r_i ≥ r_min,util (modo auto)."""
    la = caso.lastre
    if la.x_inicio is not None:
        return la.x_inicio
    idx = np.nonzero(cav.r_i >= la.r_min_util)[0]
    return float(cav.x[idx[0]]) if idx.size else math.nan


def cabe_cavidad_electronica(caso: Caso, cav: Cavidad, tol: float = 1e-9) -> bool | None:
    """¿El cono de la cavidad de electrónica queda dentro de la pared (r_i ≥ r_cavidad en todo su
    largo)? La cavidad real es la intersección del cono con el interior, así que un False solo
    indica que la pared la recorta (es informativo). None si no hay cavidad declarada."""
    el = caso.electronica
    if el.D_cav is None or el.x_inicio is None:
        return None
    x0, x1 = el.x_inicio, el.x_inicio + el.Le
    if x0 < 0 or x1 > cav.x_fin_cavidad + tol:
        return False
    m = (cav.x > x0) & (cav.x < x1)
    xs = np.r_[x0, cav.x[m], x1]
    return bool(np.all(np.interp(xs, cav.x, cav.r_i) >= el.r_cavidad(xs) - tol))


@dataclass(frozen=True)
class LimiteGeo:
    ell_geo: float
    limitante: str
    candidatos: dict[str, float]


def limite_geometrico(caso: Caso, cav: Cavidad, x_b0: float) -> LimiteGeo:
    g, la, el = caso.geom, caso.lastre, caso.electronica
    if not math.isfinite(x_b0):
        return LimiteGeo(0.0, "sin_cavidad_util", {})
    cand: dict[str, float] = {"fin_cavidad": cav.x_fin_cavidad - x_b0}

    zonas = [(z0, z1, f"zona_prohibida[{i}]") for i, (z0, z1) in enumerate(la.zonas_prohibidas)]
    if el.modo == "fija":
        zonas.append((el.x_inicio, el.x_inicio + el.Le, "electronica_fija"))
    if not la.permitir_en_cola:
        zonas.append((g.x_cola - la.margen_cola, math.inf, "cola"))
    for z0, z1, nombre in zonas:
        if z1 <= x_b0:
            continue
        cand[nombre] = max(z0 - x_b0, 0.0)

    if el.modo == "detras_del_lastre":
        cand["electronica_detras"] = (g.x_cola - la.margen_cola) - el.Le - el.holgura - x_b0
    if la.fmax is not None:
        cand["fraccion_max_L"] = la.fmax * g.L - x_b0

    limitante = min(cand, key=cand.get)
    return LimiteGeo(ell_geo=cand[limitante], limitante=limitante, candidatos=cand)


# --------------------------------------------------------------------------- modelo CG(ℓ)


@dataclass
class ModeloLastre:
    """x_CG(ℓ) = (M_0 + ρ_b Φ_1(ℓ) + m_e x̄_e(ℓ)) / (m_0 + ρ_b ∀_b(ℓ) + m_e)."""

    caso: Caso
    cav: Cavidad
    mv: MasaVacia
    x_b0: float
    rho_b: float

    @property
    def fll(self) -> float:
        return self.caso.lastre.fll

    @property
    def ell_cavidad(self) -> float:
        """Máxima longitud que cabe físicamente en la cavidad (ignora zonas prohibidas)."""
        return max(self.cav.x_fin_cavidad - self.x_b0, 0.0)

    def V_b(self, ell):
        return self.fll * self.cav.vol(self.x_b0, self.x_b0 + np.asarray(ell, dtype=float))

    def Phi1(self, ell):
        return self.fll * self.cav.mom(self.x_b0, self.x_b0 + np.asarray(ell, dtype=float))

    def m_b(self, ell):
        return self.rho_b * self.V_b(ell)

    def xbar_b(self, ell):
        """Centroide del lastre, x̄_b(ℓ) = Φ_1/∀_b."""
        return self.Phi1(ell) / self.V_b(ell)

    def SM_inf(self, ell, x_CP: float, D_ref: float):
        """Techo por geometría para una región de lastre dada (v2 §3.6): el SM al que tiende
        SM(ℓ) cuando ρ_b → ∞ con la misma región, (x_CP − x̄_b(ℓ))/D."""
        return (x_CP - self.xbar_b(ell)) / D_ref

    def x_e(self, ell):
        el = self.caso.electronica
        ell = np.asarray(ell, dtype=float)
        if el.modo == "detras_del_lastre":
            return self.x_b0 + ell + el.holgura
        return np.full_like(ell, el.x_inicio)

    def dxe_dell(self) -> float:
        return 1.0 if self.caso.electronica.modo == "detras_del_lastre" else 0.0

    def xbar_e(self, ell):
        return self.caso.electronica.x_cg(self.x_e(ell))

    def m(self, ell):
        return self.mv.m0 + self.m_b(ell) + self.caso.electronica.me

    def x_CG(self, ell):
        me = self.caso.electronica.me
        num = self.mv.M0 + self.rho_b * self.Phi1(ell) + me * self.xbar_e(ell)
        return num / self.m(ell)

    def dxCG_dell(self, ell):
        """Derivada analítica (sección 3.8)."""
        ell = np.asarray(ell, dtype=float)
        xb = self.x_b0 + ell
        A = self.cav.area(xb)
        num = self.rho_b * self.fll * A * (xb - self.x_CG(ell)) + self.caso.electronica.me * self.dxe_dell()
        return num / self.m(ell)

    # --- tapón trasero: detrás de la electrónica, desde x_r0(ℓ) hacia popa
    def x_r0(self, ell) -> float:
        el = self.caso.electronica
        return float(self.x_e(ell)) + el.Le + el.holgura

    def ell2_max(self, ell) -> float:
        return max(self.cav.x_fin_cavidad - self.caso.lastre.margen_popa - self.x_r0(ell), 0.0)

    def V_tras(self, ell, ell2):
        a = self.x_r0(ell)
        return self.fll * self.cav.vol(a, a + np.asarray(ell2, dtype=float))

    def m_con_trasero(self, ell, ell2):
        return self.m(ell) + self.rho_b * self.V_tras(ell, ell2)

    def x_CG_con_trasero(self, ell, ell2):
        a = self.x_r0(ell)
        M_tras = self.rho_b * self.fll * self.cav.mom(a, a + np.asarray(ell2, dtype=float))
        return (self.x_CG(ell) * self.m(ell) + M_tras) / self.m_con_trasero(ell, ell2)

    def ell2_por_SM(self, ell, x_CP: float, D_ref: float, SM_min: float, tol: float) -> float:
        """Máximo ℓ_2 ∈ [0, ℓ_2,max] con SM ≥ SM_min. El tapón trasero está detrás del CG, así
        que SM(ℓ_2) decrece monótonamente y basta un brentq. 0 si ni ℓ_2 = 0 cumple."""
        l2max = self.ell2_max(ell)

        def g(l2):
            return (x_CP - float(self.x_CG_con_trasero(ell, l2))) / D_ref - SM_min

        if l2max <= 0 or g(0.0) < 0:
            return 0.0
        if g(l2max) >= 0:
            return l2max
        return brentq(g, 0.0, l2max, xtol=tol)

    def ell_de_masa(self, m_obj: float, tol: float) -> float:
        """Resuelve ρ_b ∀_b(ℓ) = m_obj en [0, ℓ_cavidad]; NaN si no cabe en la cavidad."""
        ell_c = self.ell_cavidad
        if m_obj <= 0:
            return 0.0
        if self.m_b(ell_c) < m_obj:
            return math.nan
        return brentq(lambda l: float(self.m_b(l)) - m_obj, 0.0, ell_c, xtol=tol)


# --------------------------------------------------------------------------- ventana de SM


@dataclass
class Ventana:
    ell_star: float
    x_CG_min: float
    intervalos: list[tuple[float, float]]
    unimodal: bool
    ell_SM_min: float = math.nan
    ell_SM_max: float = math.nan
    banderas: list[str] = field(default_factory=list)


def _n_valles(y: np.ndarray, tol: float) -> int:
    """Número de mínimos locales interiores de y, más el extremo derecho si la curva llega
    bajando. El extremo izquierdo no cuenta: si allí se cumple el SM y luego deja de
    cumplirse, eso ya aparece como ventana partida en `intervalos`.

    Con la electrónica detrás del lastre, x_CG(ℓ) suele subir unas décimas de mm en los
    primeros milímetros (en la punta A_i es tan pequeña que ρ_b A_i (x_b0 − x_CG) + m_e > 0)
    y luego bajar hasta el mínimo ℓ*. Esa joroba inicial no crea un segundo valle ni parte la
    ventana, así que no se considera pérdida de unimodalidad; sí lo es un segundo valle.
    """
    d = np.diff(y)
    s = np.sign(np.where(np.abs(d) <= tol, 0.0, d))
    s = s[s != 0]
    if s.size == 0:
        return 1
    return int(np.sum((s[:-1] < 0) & (s[1:] > 0))) + int(s[-1] < 0)


def ventana_SM(mod: ModeloLastre, ell_geo: float, x_CP: float, D_ref: float, SM_min: float,
               SM_max: float | None, n: int, tol: float) -> Ventana:
    """Conjunto de ℓ ∈ [0, ℓ_geo] con SM_min ≤ SM(ℓ) ≤ SM_max.

    ℓ_SM,min y ℓ_SM,max son el ínfimo y el supremo de ese conjunto; si no es un intervalo
    único, o x_CG(ℓ) tiene más de un valle, se levanta `no_unimodal` y todos los intervalos
    quedan en `intervalos`.
    """
    ells = np.linspace(0.0, ell_geo, n)
    xcg = mod.x_CG(ells)
    i = int(np.argmin(xcg))
    lo, hi = ells[max(i - 1, 0)], ells[min(i + 1, n - 1)]
    if hi > lo:
        r = minimize_scalar(lambda l: float(mod.x_CG(l)), bounds=(lo, hi), method="bounded",
                            options={"xatol": tol})
        ell_star, xmin = (float(r.x), float(r.fun)) if r.fun <= xcg[i] else (ells[i], xcg[i])
    else:
        ell_star, xmin = float(ells[i]), float(xcg[i])

    x_max = x_CP - SM_min * D_ref
    x_min = -math.inf if SM_max is None else x_CP - SM_max * D_ref

    def h(l):
        x = mod.x_CG(l)
        return np.maximum(x - x_max, x_min - x)

    hv = h(ells)
    ok = hv <= 0
    intervalos: list[tuple[float, float]] = []
    inicio = 0.0 if ok[0] else None
    for k in range(1, n):
        if ok[k] != ok[k - 1]:
            raiz = brentq(lambda l: float(h(l)), ells[k - 1], ells[k], xtol=tol)
            if ok[k]:
                inicio = raiz
            else:
                intervalos.append((inicio, raiz))
                inicio = None
    if inicio is not None:
        intervalos.append((inicio, float(ells[-1])))

    unimodal = _n_valles(xcg, tol * 1e-3) <= 1
    v = Ventana(ell_star=ell_star, x_CG_min=xmin, intervalos=intervalos, unimodal=unimodal)
    if not intervalos:
        v.banderas.append("SM_inalcanzable")
    else:
        v.ell_SM_min, v.ell_SM_max = intervalos[0][0], intervalos[-1][1]
    if not unimodal or len(intervalos) > 1:
        v.banderas.append("no_unimodal")
    return v


# --------------------------------------------------------------------------- trim


def alpha_trim(m, x_CG, x_T, x_CP, q, S_ref, CNa):
    """α_trim ≈ m g (x_CG − x_T) / (q S_ref C_Nα (x_CP − x_T)) [rad]; NaN si x_CP ≤ x_T."""
    brazo = np.asarray(x_CP - x_T, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = m * G0 * (x_CG - x_T) / (q * S_ref * CNa * brazo)
    return np.where(brazo > 0, a, np.nan)
