"""Aleta freeform anclada a la cola (v2 §3.3–3.4): puntos, polígono global, área, centroide,
masa y envolvente.

Coordenadas locales de OpenRocket: origen en el borde de ataque de la raíz, x hacia popa,
y radial medido desde el radio del cuerpo en ese borde de ataque (r_LE). El cierre P5 → P0 lo
hace OpenRocket siguiendo la superficie de la transición; aquí se reproduce muestreando r_e(x).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

MM = 1e-3


@dataclass(frozen=True)
class ParamsAleta:
    """Parámetros de la aleta ya resueltos a SI (salvo `h` si se dio como radio de punta)."""

    n: int
    rotacion: float  # rad, φ_0 medido desde la vertical
    c_r: float
    x_s: float
    h: float | None  # altura de la punta sobre r_LE
    r_tip_obj: float | None  # alternativa: radio absoluto de la punta
    e: float  # extensión detrás de la base
    r_in: float
    delta_b: float  # offset axial *bottom*
    eps: float
    t: float
    material: str
    rho: float
    phi: float
    descontar_solape_eje: bool = False


@dataclass(frozen=True)
class GeomAleta:
    params: ParamsAleta
    x_LE: float
    r_LE: float
    R_a: float
    h: float
    puntos: np.ndarray  # (6, 2) locales P0..P5
    poligono: np.ndarray  # (k, 2) global (x, r), cerrado implícitamente
    origen: tuple[str, ...]  # 'punto' | 'superficie' por vértice del polígono

    @property
    def r_tip(self) -> float:
        return self.r_LE + self.h

    @property
    def area(self) -> float:
        return area_centroide(self.poligono)[0]

    @property
    def x_cg(self) -> float:
        return area_centroide(self.poligono)[1]

    @property
    def x_TE(self) -> float:
        return self.x_LE + self.params.c_r + self.params.e

    @property
    def masa(self) -> float:
        p = self.params
        V = p.n * p.t * self.area - self.volumen_solape()
        return p.rho * p.phi * V

    def volumen_solape(self) -> float:
        """Volumen contado de más donde las n aletas se cruzan cerca del eje (opcional).

        Si r_in < t/2, cada aleta mete una franja de alto (t/2 − r_in) en el núcleo común; el
        núcleo (≈ t² por unidad de longitud) se cuenta una sola vez. Longitud: borde interior P4–P3.
        """
        p = self.params
        if not p.descontar_solape_eje or p.r_in >= p.t / 2:
            return 0.0
        Lx = self.puntos[3, 0] - self.puntos[4, 0]
        return max(p.n * p.t * (p.t / 2 - p.r_in) * Lx - p.t**2 * Lx, 0.0)


# --------------------------------------------------------------------------- polígono


def area_centroide(P: np.ndarray) -> tuple[float, float]:
    """Área (valor absoluto) y x̄ de un polígono simple por la fórmula del polígono (shoelace)."""
    x, y = P[:, 0], P[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cruz = x * y1 - x1 * y
    A = 0.5 * cruz.sum()
    xc = ((x + x1) * cruz).sum() / (6.0 * A)
    return abs(A), float(xc)


def _segs_cortan(p1, p2, p3, p4, tol=1e-12) -> bool:
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    d1, d2 = orient(p3, p4, p1), orient(p3, p4, p2)
    d3, d4 = orient(p1, p2, p3), orient(p1, p2, p4)
    return (d1 * d2 < -tol) and (d3 * d4 < -tol)


def es_simple(P: np.ndarray) -> bool:
    """True si ningún par de aristas no contiguas se cruza."""
    n = len(P)
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            if _segs_cortan(a, b, P[j], P[(j + 1) % n]):
                return False
    return True


def _sin_duplicados(P: np.ndarray, orig: list[str], tol: float) -> tuple[np.ndarray, list[str]]:
    keep = [0]
    for i in range(1, len(P)):
        if np.hypot(*(P[i] - P[keep[-1]])) > tol:
            keep.append(i)
    if len(keep) > 1 and np.hypot(*(P[keep[-1]] - P[keep[0]])) <= tol:
        keep.pop()
    return P[keep], [orig[k] for k in keep]


def poligono_global(puntos_locales: np.ndarray, x_LE: float, r_LE: float,
                    r_sup: Callable[[np.ndarray], np.ndarray], n_sup: int, tol: float):
    """Transforma P0…Pk a (x, r) globales y agrega el tramo de superficie Pk → P0."""
    P = np.column_stack([x_LE + puntos_locales[:, 0], r_LE + puntos_locales[:, 1]])
    orig = ["punto"] * len(P)
    x_fin = P[-1, 0]
    xs = np.linspace(x_fin, x_LE, n_sup)[1:-1]
    S = np.column_stack([xs, r_sup(xs)])
    return _sin_duplicados(np.vstack([P, S]), orig + ["superficie"] * len(S), tol)


# --------------------------------------------------------------------------- construcción


def construir(p: ParamsAleta, x_t0: float, Lt: float, r_sup: Callable[[np.ndarray], np.ndarray],
              n_sup: int, tol_superficie: float, errores: list[str] | None = None) -> GeomAleta:
    """Puntos P0…P5 (v2 §3.3) y polígono global a partir de los parámetros y del perfil r_sup(x).

    r_sup debe prolongar el radio de popa más allá de la base (ver perfiles.radio_superficie).
    """
    err = errores if errores is not None else []
    x_LE = x_t0 + Lt + p.delta_b - p.c_r
    r_LE = float(r_sup(np.array([x_LE]))[0])
    R_a = float(r_sup(np.array([x_t0 + Lt]))[0])
    h = p.h if p.h is not None else p.r_tip_obj - r_LE

    if not p.c_r <= Lt + 1e-12:
        err.append(f"aletas: c_r = {p.c_r / MM:.2f} mm debe ser ≤ L_t = {Lt / MM:.2f} mm")
    if not h > 0:
        err.append(f"aletas: h_tip = {h / MM:.2f} mm debe ser > 0")
    if not 0 <= p.r_in < R_a:
        err.append(f"aletas: r_interior = {p.r_in / MM:.3f} mm debe cumplir 0 ≤ r_in < R_a = {R_a / MM:.2f} mm")
    if not 0 <= p.x_s <= p.c_r + p.e:
        err.append(f"aletas: x_tip_le = {p.x_s / MM:.2f} mm debe estar en [0, c_r + e]")
    if not p.e >= 0:
        err.append("aletas: la extensión debe ser ≥ 0")

    y5 = R_a - r_LE
    y_in = p.r_in - r_LE
    cre = p.c_r + p.e
    if p.e > 0:
        pts = [(0, 0), (p.x_s, h), (cre, h), (cre, y_in), (p.c_r + p.eps, y_in), (p.c_r, y5)]
    else:  # sin extensión: P3, P4 y P5 colapsan en el borde de salida de la raíz
        pts = [(0, 0), (p.x_s, h), (p.c_r, h), (p.c_r, y5)]
    puntos = np.array(pts, dtype=float)
    y5_sup = float(r_sup(np.array([x_LE + p.c_r]))[0]) - r_LE
    if abs(puntos[-1, 1] - y5_sup) >= tol_superficie:
        err.append(f"aletas: P5 fuera de la superficie ({(puntos[-1, 1] - y5_sup) / MM:+.4f} mm)")
    poli, orig = poligono_global(puntos, x_LE, r_LE, r_sup, n_sup, tol=1e-9)
    if not es_simple(poli):
        err.append("aletas: el polígono de la aleta se autointerseca")
    if errores is None and err:
        raise ValueError("; ".join(err))
    return GeomAleta(params=p, x_LE=x_LE, r_LE=r_LE, R_a=R_a, h=h, puntos=puntos, poligono=poli,
                     origen=tuple(orig))


def desde_poligono(p: ParamsAleta, poligono: np.ndarray, origen: list[str], x_LE: float, r_LE: float,
                   R_a: float) -> GeomAleta:
    """GeomAleta a partir del polígono global exportado por OpenRocket (or_aletas.csv)."""
    es_pto = np.array([o == "punto" for o in origen])
    pts = poligono[es_pto] - np.array([x_LE, r_LE])
    h = float(pts[:, 1].max())
    return GeomAleta(params=p, x_LE=x_LE, r_LE=r_LE, R_a=R_a, h=h, puntos=pts, poligono=poligono,
                     origen=tuple(origen))


# --------------------------------------------------------------------------- envolvente


def angulos(n: int, rotacion: float) -> np.ndarray:
    return rotacion + 2 * np.pi * np.arange(n) / n


def envolvente(r_tip: float, R: float, n: int, rotacion: float) -> tuple[float, float]:
    """(alto, ancho) de la sección transversal: cuerpo de radio R y n aletas hasta r_tip.

    Con φ medido desde la vertical, cada aleta llega a r_tip·cos φ en vertical y r_tip·sin φ en
    horizontal. Para n par simétrico se reduce a 2 r_tip max|cos φ| (v2 §3.3); la forma general
    también vale para n impar (p. ej. 3 aletas a 0°: 1.5 r_tip, no 2 r_tip).
    """
    phi = angulos(n, rotacion)
    c, s = r_tip * np.cos(phi), r_tip * np.sin(phi)
    alto = max(R, c.max()) + max(R, -c.min())
    ancho = max(R, s.max()) + max(R, -s.min())
    return float(alto), float(ancho)


# --------------------------------------------------------------------------- Barrowman (Niskanen §3.2.2)


@dataclass(frozen=True)
class CPAleta:
    CNa: float  # del conjunto, referido a A_ref
    x_CP: float  # global
    span: float
    area: float
    mac: float
    x_mac_le: float  # global
    cos_gamma: float


def barrowman_freeform(g: GeomAleta, A_ref: float, mach: float, n_franjas: int) -> CPAleta:
    """CP y C_Nα del conjunto de aletas por franjas en y, como FinSetCalc de OpenRocket.

    El span va desde el punto más bajo de la raíz (P5, sobre la cola) hasta la punta, igual que
    `FinSet.getSpan()` de OpenRocket. La cuerda de cada franja va del borde de ataque al de
    salida (incluye huecos), el área que entra en C_Nα es el área real del polígono y el ángulo
    de flecha es el de la línea de cuerdas medias, promediado por franjas.
    """
    P = g.poligono - np.array([g.x_LE, g.r_LE])  # local
    y_lo = g.R_a - g.r_LE
    y_hi = g.h
    span = y_hi - y_lo
    ys = np.linspace(y_lo, y_hi, n_franjas)
    lead = np.full(n_franjas, np.inf)
    trail = np.full(n_franjas, -np.inf)
    for (x1, y1), (x2, y2) in zip(P, np.roll(P, -1, axis=0)):
        if abs(y2 - y1) < 1e-12:
            continue
        lo, hi = min(y1, y2), max(y1, y2)
        m = (ys >= lo - 1e-12) & (ys <= hi + 1e-12)
        x = x1 + (ys[m] - y1) * (x2 - x1) / (y2 - y1)
        lead[m] = np.minimum(lead[m], x)
        trail[m] = np.maximum(trail[m], x)
    ok = np.isfinite(lead) & np.isfinite(trail)
    c = np.where(ok, trail - lead, 0.0)
    lead = np.where(ok, lead, 0.0)
    dy = span / (n_franjas - 1)
    w = np.full(n_franjas, dy)
    w[[0, -1]] = dy / 2  # trapecio
    A_franjas = float((c * w).sum())
    mac = float((c**2 * w).sum() / A_franjas)
    x_mac_le = float((lead * c * w).sum() / A_franjas)
    xm = lead + c / 2
    tan_g = np.gradient(xm[ok], ys[ok])
    cos_g = float(np.mean(np.cos(np.arctan(tan_g))))
    A_f = g.area
    beta = math.sqrt(max(1 - mach**2, 1e-6))
    CNa1 = 2 * math.pi * span**2 / A_ref / (1 + math.sqrt(1 + (beta * span**2 / (A_f * cos_g)) ** 2))
    n = g.params.n
    factor_n = n / 2 if n >= 3 else 1.0
    K_TB = 1 + g.r_LE / (span + g.r_LE)
    return CPAleta(CNa=CNa1 * factor_n * K_TB, x_CP=g.x_LE + x_mac_le + 0.25 * mac, span=span,
                   area=A_f, mac=mac, x_mac_le=g.x_LE + x_mac_le, cos_gamma=cos_g)
