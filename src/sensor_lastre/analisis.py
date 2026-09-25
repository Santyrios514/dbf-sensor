"""Análisis completo de una configuración (flujo de la sección 6.2), independiente de la E/S."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .barrowman import ResultadoCP, cp_interno
from .config import Caso
from .estabilidad import (LimiteGeo, ModeloLastre, Ventana, alpha_trim, limite_geometrico,
                          ventana_SM, x_inicio_lastre)
from .geometria import Cavidad, construir_cavidad, malla, perfil_desde_muestras
from .masas import MasaVacia, masa_vacia
from .materiales import Relleno

MM, G, CM3 = 1e-3, 1e-3, 1e-6
NAN = math.nan


@dataclass(frozen=True)
class AeroOR:
    """x_CP y C_Nα exportados por OpenRocket (script 1)."""

    x_CP: float  # m
    CNa: float


@dataclass
class ResultadoRelleno:
    relleno: Relleno
    modelo: ModeloLastre
    ventana: Ventana | None
    fila: dict
    curva: pd.DataFrame
    presupuestos: list[dict]


@dataclass
class ResultadoCaso:
    caso: Caso
    cav: Cavidad
    mv: MasaVacia
    cp_int: ResultadoCP
    x_CP: float
    CNa: float
    x_CP_fuente: str
    x_b0: float
    lim: LimiteGeo
    banderas: list[str]
    rellenos: list[ResultadoRelleno] = field(default_factory=list)

    def filas_masas_capas(self) -> list[dict]:
        return [{
            "config_id": self.caso.id, "estacion": c.estacion, "capa_idx": c.capa_idx,
            "material": c.material, "espesor_mm": c.espesor / MM, "fraccion_solida": c.phi,
            "m_g": c.m / G, "x_cg_mm": c.x_cg / MM,
        } for c in self.mv.capas]


def banderas_envolvente(caso: Caso) -> tuple[bool, bool]:
    g, a = caso.geom, caso.geom.aletas
    borde_salida = max(a.x_r0 + a.cr, a.x_r0 + a.xs + a.ct)
    delta_popa = max(0.0, borde_salida - g.L)
    cabe_largo = g.L + delta_popa <= caso.envolvente.largo_max
    cabe_diametro = g.D + 2 * a.s <= caso.envolvente.diametro_max
    return bool(cabe_largo), bool(cabe_diametro)


def _txt(banderas: list[str]) -> str:
    return ";".join(dict.fromkeys(banderas))


def analizar_caso(caso: Caso, rellenos: list[Relleno] | None = None,
                  perfil: pd.DataFrame | None = None, aero_or: AeroOR | None = None,
                  tol_dif_CP_frac_L: float | None = None) -> ResultadoCaso:
    """perfil: filas de or_perfiles.csv de este caso (None → perfil analítico)."""
    g, nu = caso.geom, caso.numerico
    x = malla(g.L, nu.dx)
    r_e = None
    if perfil is not None and len(perfil):
        r_e = perfil_desde_muestras(x, perfil["x_mm"].to_numpy() * MM, perfil["r_ext_mm"].to_numpy() * MM)
    cav = construir_cavidad(caso, r_e=r_e, x=x)
    mv = masa_vacia(caso, cav)
    cp_int = cp_interno(g, cav.x, cav.r_e)

    banderas: list[str] = []
    if aero_or is not None and math.isfinite(aero_or.x_CP) and math.isfinite(aero_or.CNa):
        x_CP, CNa, fuente = aero_or.x_CP, aero_or.CNa, "openrocket"
        if tol_dif_CP_frac_L is not None and abs(x_CP - cp_int.x_CP) > tol_dif_CP_frac_L * g.L:
            banderas.append("dif_CP_alta")
    else:
        x_CP, CNa, fuente = cp_int.x_CP, cp_int.CNa, "interno"

    x_b0 = x_inicio_lastre(caso, cav)
    lim = limite_geometrico(caso, cav, x_b0)
    if not lim.ell_geo > 0:
        banderas.append("inviable_geo")

    res = ResultadoCaso(caso=caso, cav=cav, mv=mv, cp_int=cp_int, x_CP=x_CP, CNa=CNa,
                        x_CP_fuente=fuente, x_b0=x_b0, lim=lim, banderas=banderas)
    for rel in (rellenos if rellenos is not None else list(caso.rellenos)):
        res.rellenos.append(_analizar_relleno(res, rel))
    return res


def _analizar_relleno(res: ResultadoCaso, rel: Relleno) -> ResultadoRelleno:
    caso, cav, mv = res.caso, res.cav, res.mv
    g, nu, vu = caso.geom, caso.numerico, caso.vuelo
    D_ref, S_ref = g.D, math.pi * g.D**2 / 4
    x_CP, CNa = res.x_CP, res.CNa
    x_b0 = res.x_b0 if math.isfinite(res.x_b0) else 0.0
    mod = ModeloLastre(caso=caso, cav=cav, mv=mv, x_b0=x_b0, rho_b=rel.rho_b)
    ell_geo = res.lim.ell_geo
    banderas = list(res.banderas)
    cabe_largo, cabe_diam = banderas_envolvente(caso)

    def SM(xcg):
        return (x_CP - xcg) / D_ref

    def x_T_de(xcg):
        return xcg if caso.remolque.modo == "en_CG" else caso.remolque.x_T

    def trim_deg(ell):
        xcg = mod.x_CG(ell)
        return np.degrees(alpha_trim(mod.m(ell), xcg, x_T_de(xcg), x_CP, vu.q, S_ref, CNa))

    V_ext, V_int = cav.V_ext, cav.V_int
    fila = {
        "config_id": caso.id, "relleno": rel.nombre, "rho_b_kg_m3": rel.rho_b,
        "L_mm": g.L / MM, "D_mm": g.D / MM, "V_ext_cm3": V_ext / CM3, "V_int_cm3": V_int / CM3,
        "x_b0_mm": res.x_b0 / MM, "ell_geo_mm": ell_geo / MM,
        "x_CP_mm": x_CP / MM, "x_CP_fuente": res.x_CP_fuente,
        "dx_CP_or_vs_interno_mm": (x_CP - res.cp_int.x_CP) / MM if res.x_CP_fuente == "openrocket" else NAN,
        "CN_alpha_rad": CNa,
        "m_casco_g": mv.m_casco / G, "m_aletas_g": mv.m_aletas / G,
        "m_mamparos_g": mv.m_mamparos / G, "m_electronica_g": caso.electronica.me / G,
        "m_puntuales_g": mv.m_puntuales / G, "cabe_largo": cabe_largo, "cabe_diametro": cabe_diam,
    }

    ventana = None
    curva = pd.DataFrame(columns=["ell_mm", "V_b_cm3", "m_b_g", "m_total_g", "x_CG_mm", "SM_cal",
                                  "alpha_trim_deg"])
    if ell_geo > 0:
        Vg = float(mod.V_b(ell_geo))
        xcg_g = float(mod.x_CG(ell_geo))
        fila.update({
            "V_util_geo_cm3": Vg / CM3, "eta_int_geo": Vg / V_int, "eta_ext_geo": Vg / V_ext,
            "m_relleno_geo_g": rel.rho_b * Vg / G, "m_total_geo_g": float(mod.m(ell_geo)) / G,
            "x_CG_geo_mm": xcg_g / MM, "SM_geo_cal": SM(xcg_g),
            "alpha_trim_geo_deg": float(trim_deg(ell_geo)),
        })
        ventana = ventana_SM(mod, ell_geo, x_CP, D_ref, caso.estabilidad.SM_min,
                             caso.estabilidad.SM_max, nu.n_ell, nu.tol)
        banderas += ventana.banderas
        fila.update({
            "ell_star_mm": ventana.ell_star / MM, "x_CG_min_mm": ventana.x_CG_min / MM,
            "SM_max_alcanzable_cal": SM(ventana.x_CG_min),
            "ell_SM_min_mm": ventana.ell_SM_min / MM, "ell_SM_max_mm": ventana.ell_SM_max / MM,
        })
        ell_u = ventana.ell_SM_max
        if math.isfinite(ell_u):
            Vu = float(mod.V_b(ell_u))
            xcg_u = float(mod.x_CG(ell_u))
            fila.update({
                "V_util_cm3": Vu / CM3, "eta_int": Vu / V_int, "eta_ext": Vu / V_ext,
                "m_relleno_g": rel.rho_b * Vu / G, "m_total_g": float(mod.m(ell_u)) / G,
                "x_CG_mm": xcg_u / MM, "SM_cal": SM(xcg_u), "x_T_mm": x_T_de(xcg_u) / MM,
                "alpha_trim_deg": float(trim_deg(ell_u)), "x_T_trim_cero_mm": xcg_u / MM,
            })
        ells = np.linspace(0.0, ell_geo, nu.n_ell)
        xcg = mod.x_CG(ells)
        curva = pd.DataFrame({
            "ell_mm": ells / MM, "V_b_cm3": mod.V_b(ells) / CM3, "m_b_g": mod.m_b(ells) / G,
            "m_total_g": mod.m(ells) / G, "x_CG_mm": xcg / MM, "SM_cal": SM(xcg),
            "alpha_trim_deg": trim_deg(ells),
        })
    if caso.remolque.modo == "absoluto":
        fila.setdefault("x_T_mm", caso.remolque.x_T / MM)
        if x_CP <= caso.remolque.x_T:
            banderas.append("inestable_respecto_remolque")

    presupuestos = []
    for m_obj in caso.presupuesto:
        ell = mod.ell_de_masa(m_obj, nu.tol)
        no_cabe = (not math.isfinite(ell)) or ell > ell_geo
        p = {"config_id": caso.id, "relleno": rel.nombre, "m_lastre_obj_g": m_obj / G,
             "ell_mm": ell / MM, "no_cabe": bool(no_cabe)}
        if math.isfinite(ell):
            xcg = float(mod.x_CG(ell))
            p.update({"m_total_g": float(mod.m(ell)) / G, "x_CG_mm": xcg / MM, "SM_cal": SM(xcg),
                      "alpha_trim_deg": float(trim_deg(ell))})
        presupuestos.append(p)

    fila["banderas"] = _txt(banderas)
    return ResultadoRelleno(relleno=rel, modelo=mod, ventana=ventana, fila=fila, curva=curva,
                            presupuestos=presupuestos)
