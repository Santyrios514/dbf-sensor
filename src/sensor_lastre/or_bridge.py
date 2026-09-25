"""Todo el acceso a OpenRocket (orlab + JPype) vive aquí. El resto del código es Python puro.

Firmas verificadas contra OpenRocket 24.12 (paquete raíz `info.openrocket.core`) y
documentadas en docs/api_openrocket.md. Con OpenRocket ≤ 23.09 el paquete raíz es
`net.sf.openrocket`; orlab lo resuelve y lo expone como `inst.openrocket`.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import ESTACIONES, Caso
from .materiales import densidad_equivalente, espesor_total
from .perfiles import FORMAS
from .verificacion import FORMA_OR_A_YAML, ValoresOrk

MM = 1e-3
COMPONENTES = ("nariz", "cuerpo", "cola", "aletas")


class ErrorAleta(RuntimeError):
    """OpenRocket no aceptó los puntos de la aleta tal como se pidieron."""


@dataclass
class ResultadoAero:
    x_CP: float
    CNa: float
    CNa_dif: float
    CD: float
    por_componente: dict[str, tuple[float, float]]  # nombre -> (CNa, x_CP global)
    warnings: list[str]


@dataclass
class ResultadoMasas:
    m_total: float
    x_CG_total: float
    por_componente: dict[str, tuple[float, float]] = field(default_factory=dict)  # (m, x_CG global)


def sha256(ruta: str | Path) -> str:
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class PuenteOR:
    """Envoltorio de una instancia de OpenRocket abierta con orlab.

    Uso:
        with orlab.OpenRocketInstance(jar_path=...) as inst:
            pr = PuenteOR(inst, nombres)
            pr.cargar("modelos/analisis_vol_int.ork")
    """

    def __init__(self, inst, nombres: dict[str, str]):
        import jpype
        import orlab

        self._jpype = jpype
        self.inst = inst
        self.orl = orlab.Helper(inst)
        self.core = inst.openrocket
        self.nombres = dict(nombres)
        self.doc = self.rocket = None
        self.comp: dict[str, object] = {}

    # ------------------------------------------------------------------ metadatos
    @property
    def version_or(self) -> str:
        return str(self.core.util.BuildProperties.getVersion())

    @staticmethod
    def version_orlab() -> str:
        try:
            return importlib.metadata.version("orlab")
        except importlib.metadata.PackageNotFoundError:
            return ""

    # ------------------------------------------------------------------ carga
    def cargar(self, ruta: str | Path):
        self.doc = self.orl.load_doc(str(ruta))
        self.rocket = self.doc.getRocket()
        encontrados = {}
        it = self.rocket.iterator(True)
        while it.hasNext():
            c = it.next()
            encontrados.setdefault(str(c.getName()), c)
        faltan = [n for n in self.nombres.values() if n not in encontrados]
        if faltan:
            raise KeyError(f"componentes {faltan} no están en {ruta}; hay {sorted(encontrados)}")
        self.comp = {k: encontrados[v] for k, v in self.nombres.items()}
        tipos = {"nariz": "NoseCone", "cuerpo": "BodyTube", "cola": "Transition", "aletas": "FreeformFinSet"}
        for k, t in tipos.items():
            real = str(self.comp[k].getClass().getSimpleName())
            if real != t:
                raise TypeError(f"'{self.nombres[k]}' es {real}, se esperaba {t}")
        return self

    @property
    def configuracion(self):
        return self.rocket.getSelectedConfiguration()

    def x_abs(self, clave: str) -> float:
        """Posición axial absoluta del inicio del componente (primera instancia)."""
        return float(self.comp[clave].getComponentLocations()[0].x)

    # ------------------------------------------------------------------ lectura de la base
    def valores(self) -> ValoresOrk:
        n, c, t, f = (self.comp[k] for k in COMPONENTES)

        def param(comp):
            sh = comp.getShapeType()
            return float(comp.getShapeParameter()) if sh.usesParameter() else None

        def forma(comp):
            return FORMA_OR_A_YAML[str(comp.getShapeType().name()).lower()]

        pts = [(float(p.x), float(p.y)) for p in f.getFinPoints()]
        cerrada = bool(t.isAftShoulderCapped()) and float(t.getAftShoulderLength()) > 0
        return ValoresOrk(
            Ln=float(n.getLength()), R_nariz=float(n.getBaseRadius()), t_nariz=float(n.getThickness()),
            forma_nariz=forma(n), param_nariz=param(n),
            Lc=float(c.getLength()), R_cuerpo=float(c.getOuterRadius()), t_cuerpo=float(c.getThickness()),
            Lt=float(t.getLength()), Rf=float(t.getForeRadius()), Ra=float(t.getAftRadius()),
            t_cola=float(t.getThickness()), forma_cola=forma(t), param_cola=param(t),
            popa_abierta=not cerrada,
            recortada=bool(t.isClipped()) if t.getShapeType().isClippable() else None, fin_n=int(f.getFinCount()), fin_t=float(f.getThickness()),
            fin_rot=float(f.getBaseRotation()), fin_offset_bottom=float(f.getAxialOffset()),
            fin_puntos=pts,
            rho={k: float(self.comp[k].getMaterial().getDensity()) for k in COMPONENTES},
        )

    # ------------------------------------------------------------------ escritura
    def _material(self, nombre: str, rho: float):
        M = self.core.material.Material
        return M.newMaterial(M.Type.BULK, nombre, float(rho), True)

    def _forma(self, forma: str):
        return self.core.rocketcomponent.Transition.Shape.valueOf(FORMAS[forma][0])

    def aplicar(self, caso: Caso, eps_puntos: float = 1e-7) -> dict:
        """Fija la geometría del caso en el cohete cargado. Devuelve r_LE, x_LE y los puntos.

        r_LE se lee de OpenRocket (radio de la transición en x_LE local) para que P5 caiga
        exactamente sobre la superficie que usa OpenRocket. Lanza ErrorAleta si OpenRocket
        modifica los puntos (lo hace en silencio con puntos inválidos, en lugar de rechazarlos).
        """
        g, pa = caso.geom, caso.params_aleta
        n, c, t, f = (self.comp[k] for k in COMPONENTES)

        n.setShapeType(self._forma(g.nariz.forma))
        if g.nariz.parametro is not None and n.getShapeType().usesParameter():
            n.setShapeParameter(float(g.nariz.parametro))
        n.setLength(float(g.nariz.Ln))
        n.setBaseRadius(float(g.R))
        c.setLength(float(g.x_cola - g.nariz.Ln))
        c.setOuterRadiusAutomatic(True)
        t.setShapeType(self._forma(g.cola.forma))
        if g.cola.parametro is not None and t.getShapeType().usesParameter():
            t.setShapeParameter(float(g.cola.parametro))
        t.setClipped(bool(g.cola.recortada))
        t.setLength(float(g.cola.Lt))
        t.setForeRadiusAutomatic(True)
        t.setAftRadius(float(g.cola.Ra))

        for k, est in (("nariz", "nariz"), ("cuerpo", "cuerpo"), ("cola", "cola")):
            capas = caso.pared[est]
            self.comp[k].setThickness(float(espesor_total(capas)))
            self.comp[k].setMaterial(self._material(f"PARED_EQ_{est}", densidad_equivalente(capas)))

        # aleta: r_LE desde OpenRocket
        Lt = float(t.getLength())
        x_LE_loc = Lt + pa.delta_b - pa.c_r
        r_LE = float(t.getRadius(x_LE_loc))
        R_a = float(t.getRadius(Lt))
        h = pa.h if pa.h is not None else pa.r_tip_obj - r_LE
        y5, y_in, cre = R_a - r_LE, pa.r_in - r_LE, pa.c_r + pa.e
        if pa.e > 0:
            pts = [(0, 0), (pa.x_s, h), (cre, h), (cre, y_in), (pa.c_r + pa.eps, y_in), (pa.c_r, y5)]
        else:
            pts = [(0, 0), (pa.x_s, h), (pa.c_r, h), (pa.c_r, y5)]
        C = self.core.util.Coordinate
        f.setPoints(self._jpype.JArray(C)([C(float(x), float(y)) for x, y in pts]))
        f.setFinCount(int(pa.n))
        f.setBaseRotation(float(pa.rotacion))
        f.setThickness(float(pa.t))
        f.setMaterial(self._material("ALETA", pa.rho * pa.phi))
        f.setAxialMethod(self.core.rocketcomponent.position.AxialMethod.BOTTOM)
        f.setAxialOffset(float(pa.delta_b))

        aceptados = [(float(p.x), float(p.y)) for p in f.getFinPoints()]
        if len(aceptados) != len(pts) or np.max(np.abs(np.array(aceptados) - np.array(pts, float))) > eps_puntos:
            raise ErrorAleta(f"OpenRocket modificó los puntos de la aleta: pedidos {np.round(np.array(pts) / MM, 3).tolist()}, "
                             f"aceptados {np.round(np.array(aceptados) / MM, 3).tolist()} (mm)")
        return {"r_LE": r_LE, "R_a": R_a, "x_LE": self.x_abs("aletas"), "h": h, "puntos": aceptados}

    # ------------------------------------------------------------------ cálculos
    def _condiciones(self, mach: float, aoa: float):
        cond = self.core.aerodynamics.FlightConditions(self.configuracion)
        cond.setMach(float(mach))
        cond.setAOA(float(aoa))
        return cond

    def aero(self, mach: float, aoa: float = 0.0, delta_aoa: float | None = None) -> ResultadoAero:
        fc = self.configuracion
        calc = self.core.aerodynamics.BarrowmanCalculator()
        ws = self.core.logging.WarningSet()
        cond = self._condiciones(mach, aoa)
        cp = calc.getCP(fc, cond, ws)
        fuerzas = calc.getAerodynamicForces(fc, cond, self.core.logging.WarningSet())
        CD = float(fuerzas.getCD())
        CNa_dif = math.nan
        if delta_aoa:
            CN = []
            for s in (1, -1):
                f = calc.getAerodynamicForces(fc, self._condiciones(mach, aoa + s * delta_aoa),
                                              self.core.logging.WarningSet())
                CN.append(float(f.getCN()))
            CNa_dif = (CN[0] - CN[1]) / (2 * delta_aoa)
        inv = {v: k for k, v in self.nombres.items()}
        por = {}
        for comp, f in calc.getForceAnalysis(fc, cond, self.core.logging.WarningSet()).items():
            nombre = str(comp.getName())
            if nombre in inv:
                c = f.getCP()
                x = float(c.x) if float(c.weight) != 0 else math.nan
                por[inv[nombre]] = (float(c.weight), x)
        return ResultadoAero(x_CP=float(cp.x), CNa=float(cp.weight), CNa_dif=CNa_dif, CD=CD,
                             por_componente=por, warnings=[str(w) for w in ws])

    def masas(self) -> ResultadoMasas:
        rb = self.core.masscalc.MassCalculator.calculateStructure(self.configuracion)
        out = ResultadoMasas(m_total=float(rb.getMass()), x_CG_total=float(rb.getCM().x))
        for k in COMPONENTES:
            c = self.comp[k]
            out.por_componente[k] = (float(c.getComponentMass()), self.x_abs(k) + float(c.getComponentCG().x))
        return out

    def perfil(self, n: int) -> list[tuple[str, float, float]]:
        """(componente, x global, r) muestreando getRadius(x local) de cada componente simétrico."""
        filas = []
        for k in ("nariz", "cuerpo", "cola"):
            c = self.comp[k]
            x0, L = self.x_abs(k), float(c.getLength())
            for xl in np.linspace(0.0, L, n):
                filas.append((k, x0 + xl, float(c.getRadius(float(xl)))))
        return filas

    def poligono_aleta(self) -> tuple[np.ndarray, list[str], float, float]:
        """Polígono global (x, r) de la aleta con el tramo de superficie, tal como lo aceptó
        OpenRocket. Devuelve (polígono, origen, x_LE, r_LE)."""
        f = self.comp["aletas"]
        x_LE = self.x_abs("aletas")
        r_LE = float(f.getBodyRadius())
        pts = [(float(p.x), float(p.y)) for p in f.getFinPoints()]
        con_raiz = [(float(p.x), float(p.y)) for p in f.getFinPointsWithRoot()]
        P, orig = [], []
        for i, (x, y) in enumerate(con_raiz):
            if P and abs(x - P[-1][0]) < 1e-12 and abs(y - P[-1][1]) < 1e-12:
                continue  # OpenRocket repite el último punto antes de la raíz
            es_pto = i < len(pts)
            P.append((x, y))
            orig.append("punto" if es_pto else "superficie")
        if len(P) > 1 and abs(P[-1][0] - P[0][0]) < 1e-12 and abs(P[-1][1] - P[0][1]) < 1e-12:
            P.pop()
            orig.pop()
        P = np.array(P) + np.array([x_LE, r_LE])
        return P, orig, x_LE, r_LE


def estaciones_T_rho(caso: Caso) -> dict[str, tuple[float, float]]:
    return {e: (espesor_total(caso.pared[e]), densidad_equivalente(caso.pared[e])) for e in ESTACIONES}
