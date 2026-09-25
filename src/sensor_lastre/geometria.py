"""Erosión multicapa del meridiano, volúmenes y momentos (secciones 3.3–3.4).

La erosión morfológica con un disco de radio T sobre el meridiano {|y| ≤ r_e(x)} da la
superficie interior a espesor normal constante T:

    r_T(x) = max(0, min_{|ζ|≤T} [ r_e(x+ζ) − sqrt(T² − ζ²) ])

Es exacta para superficies de revolución y compone: erosionar por t1 y luego por t2 equivale
a erosionar por t1 + t2 (la suma de Minkowski de dos discos es un disco).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import cumulative_simpson, simpson

from .config import ESTACIONES, Caso, Geometria, Mamparo
from .materiales import espesor_total
from .perfiles import radio_exterior


# --------------------------------------------------------------------------- malla y perfil


def malla(L: float, dx: float) -> np.ndarray:
    n = int(round(L / dx)) + 1
    return np.linspace(0.0, L, n)


def perfil_desde_muestras(x: np.ndarray, x_s: np.ndarray, r_s: np.ndarray) -> np.ndarray:
    """Interpola linealmente un perfil muestreado a la malla uniforme.

    La interpolación lineal por tramos no redondea las esquinas: los quiebres de pendiente
    (unión nariz–cuerpo, cuerpo–cola) quedan donde caen las muestras.
    """
    orden = np.argsort(x_s, kind="stable")
    return np.interp(x, np.asarray(x_s)[orden], np.asarray(r_s)[orden], left=0.0, right=0.0)


def erosionar(r: np.ndarray, T: float, dx: float, popa_cerrada: bool = True) -> np.ndarray:
    """Erosión morfológica del perfil r (malla uniforme de paso dx) con un disco de radio T.

    Fuera de la malla, r = 0 a proa (la nariz termina en punta) y a popa si la popa es cerrada;
    con popa abierta se prolonga el último radio para que no aparezca pared de fondo.
    """
    if T <= 0:
        return r.copy()
    m = int(np.floor(T / dx + 1e-9))
    derecha = np.zeros(m) if popa_cerrada else np.full(m, r[-1])
    rp = np.concatenate([np.zeros(m), r, derecha])
    n = r.size
    out = np.full(n, np.inf)
    for j in range(-m, m + 1):
        z = j * dx
        c = np.sqrt(max(T * T - z * z, 0.0))
        np.minimum(out, rp[m + j : m + j + n] - c, out=out)
    return np.maximum(out, 0.0)


def erosion_aproximada(x: np.ndarray, r: np.ndarray, T: float) -> np.ndarray:
    """Aproximación rápida r − T sqrt(1 + r'²) (solo para comparar)."""
    dr = np.gradient(r, x)
    return np.maximum(r - T * np.sqrt(1 + dr**2), 0.0)


# --------------------------------------------------------------------------- cavidad


@dataclass
class Cavidad:
    """Perfil exterior, fronteras de capas por estación y cavidad interior con sus integrales."""

    x: np.ndarray
    dx: float
    r_e: np.ndarray
    r_i: np.ndarray
    estacion: np.ndarray  # índice 0/1/2 (nariz/cuerpo/cola) de cada nodo
    fronteras: dict[str, list[np.ndarray]]  # estación -> [r_0 = r_e, r_1, ..., r_N = r_i]
    mamparos: tuple[Mamparo, ...] = ()
    C0: np.ndarray = field(init=False)  # ∫_0^x π r_i² dx
    C1: np.ndarray = field(init=False)  # ∫_0^x π r_i² x dx

    def __post_init__(self):
        A = np.pi * self.r_i**2
        self.C0 = cumulative_simpson(A, x=self.x, initial=0.0)
        self.C1 = cumulative_simpson(A * self.x, x=self.x, initial=0.0)

    # --- integrales en la cavidad (sin descontar mamparos)
    def vol_bruto(self, a, b):
        a = np.clip(a, self.x[0], self.x[-1])
        b = np.clip(b, self.x[0], self.x[-1])
        return np.interp(b, self.x, self.C0) - np.interp(a, self.x, self.C0)

    def mom_bruto(self, a, b):
        a = np.clip(a, self.x[0], self.x[-1])
        b = np.clip(b, self.x[0], self.x[-1])
        return np.interp(b, self.x, self.C1) - np.interp(a, self.x, self.C1)

    # --- descontando mamparos
    def _solape(self, a, b, f):
        tot = 0.0
        for mp in self.mamparos:
            lo = np.maximum(a, mp.x)
            hi = np.minimum(b, mp.x + mp.e)
            tot = tot + np.where(hi > lo, f(lo, hi), 0.0)
        return tot

    def vol(self, a, b):
        """Volumen libre de la cavidad en [a, b]."""
        return self.vol_bruto(a, b) - self._solape(a, b, self.vol_bruto)

    def mom(self, a, b):
        """Primer momento ∫ π r_i² x dx del volumen libre en [a, b]."""
        return self.mom_bruto(a, b) - self._solape(a, b, self.mom_bruto)

    def area(self, x):
        """Área libre A_c(x) = π r_i(x)² (0 dentro de un mamparo)."""
        x = np.asarray(x, dtype=float)
        A = np.pi * np.interp(x, self.x, self.r_i) ** 2
        for mp in self.mamparos:
            A = np.where((x >= mp.x) & (x < mp.x + mp.e), 0.0, A)
        return A

    # --- volúmenes globales
    @property
    def V_ext(self) -> float:
        return float(simpson(np.pi * self.r_e**2, x=self.x))

    @property
    def V_int_bruto(self) -> float:
        return float(simpson(np.pi * self.r_i**2, x=self.x))

    @property
    def V_mamparos(self) -> float:
        return float(sum(self.vol_bruto(m.x, m.x + m.e) for m in self.mamparos))

    @property
    def V_int(self) -> float:
        return self.V_int_bruto - self.V_mamparos

    @property
    def x_fin_cavidad(self) -> float:
        idx = np.nonzero(self.r_i > 0)[0]
        return float(self.x[idx[-1]]) if idx.size else 0.0


def indice_estacion(x: np.ndarray, g: Geometria) -> np.ndarray:
    return np.where(x < g.nariz.Ln, 0, np.where(x < g.x_cola, 1, 2))


def construir_cavidad(caso: Caso, r_e: np.ndarray | None = None, x: np.ndarray | None = None) -> Cavidad:
    """Erosión por estaciones con sus propias paredes.

    Cada estación se erosiona sobre el perfil completo con su espesor acumulado T_k; cada nodo
    toma la erosión de su estación. En una ventana ±T_max alrededor de cada unión, la cavidad
    es el mínimo de las erosiones de ambas estaciones.
    """
    g = caso.geom
    dx_nom = caso.numerico.dx
    if x is None:
        x = malla(g.L, dx_nom)
    dx = float(x[1] - x[0])
    if r_e is None:
        r_e = radio_exterior(x, g)
    cerrada = g.cola.popa_cerrada
    est = indice_estacion(x, g)

    # erosiones acumuladas por estación: E[s][k] con k = 0..N_s
    E: dict[str, list[np.ndarray]] = {}
    for s in ESTACIONES:
        capas = caso.pared[s]
        T = np.cumsum([c.t for c in capas])
        E[s] = [r_e] + [erosionar(r_e, Tk, dx, cerrada) for Tk in T]

    r_i = np.empty_like(r_e)
    for i, s in enumerate(ESTACIONES):
        r_i[est == i] = E[s][-1][est == i]
    T_max = max(espesor_total(caso.pared[s]) for s in ESTACIONES)
    for xu, (s1, s2) in ((g.nariz.Ln, ("nariz", "cuerpo")), (g.x_cola, ("cuerpo", "cola"))):
        w = np.abs(x - xu) <= T_max
        r_i[w] = np.minimum(E[s1][-1][w], E[s2][-1][w])

    fronteras = {}
    for s in ESTACIONES:
        fr = [r_e] + [np.maximum(Ek, r_i) for Ek in E[s][1:-1]] + [r_i]
        fronteras[s] = fr
    return Cavidad(x=x, dx=dx, r_e=r_e, r_i=r_i, estacion=est, fronteras=fronteras,
                   mamparos=caso.mamparos)


def integrar(y: np.ndarray, x: np.ndarray) -> float:
    return float(simpson(y, x=x))
