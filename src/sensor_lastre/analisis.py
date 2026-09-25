"""Análisis completo de una configuración (v1 §6.2 + v2 §6), independiente de la E/S."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from . import G0
from .aletas import envolvente
from .barrowman import ResultadoCP, cp_interno
from .config import Caso
from .esquemas import CURVA
from .estabilidad import (LimiteGeo, ModeloLastre, Ventana, alpha_trim, limite_geometrico,
                          ventana_SM, x_inicio_lastre)
from .geometria import Cavidad, construir_cavidad, malla, perfil_desde_muestras
from .masas import MasaVacia, masa_vacia, masas_por_componente
from .materiales import Relleno

MM, G, CM3, MM2 = 1e-3, 1e-3, 1e-6, 1e-6
NAN = math.nan


@dataclass(frozen=True)
class AeroOR:
    """x_CP y C_Nα exportados por OpenRocket (script 1)."""

    x_CP: float  # m
    CNa: float


@dataclass(frozen=True)
class Tolerancias:
    dif_CP_frac_L: float | None = None
    dif_masa_frac: float | None = None
    dif_masa_aletas_frac: float | None = None

    @classmethod
    def desde(cls, verif: dict) -> "Tolerancias":
        return cls(verif.get("dif_CP_max_frac_L"), verif.get("dif_masa_or_max_frac"),
                   verif.get("dif_masa_aletas_max_frac"))


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
    dif_masa: dict[str, float] = field(default_factory=dict)  # componente -> fracción
    dif_masa_total: float = NAN
    rellenos: list[ResultadoRelleno] = field(default_factory=list)

    def filas_masas_capas(self) -> list[dict]:
        return [{
            "config_id": self.caso.id, "estacion": c.estacion, "capa_idx": c.capa_idx,
            "material": c.material, "espesor_mm": c.espesor / MM, "fraccion_solida": c.phi,
            "m_g": c.m / G, "x_cg_mm": c.x_cg / MM,
        } for c in self.mv.capas]


def banderas_envolvente(caso: Caso) -> dict:
    """Dimensiones exteriores con aletas y si caben en la bahía (None = sin dato).

    - h_env/w_env: sección con la rotación de VUELO (la del SM).
    - D_acostado: sensor guardado acostado con las aletas en diagonal (`rotacion_guardado_deg`),
      el mayor de alto y ancho de esa sección. Es lo que ocupa en la bahía.
    - D_parado: círculo que tocan las puntas, 2·max(R, r_tip); no depende de la rotación.
    Los chequeos cabe_alto/cabe_ancho usan la posición de guardado.
    """
    g, a = caso.geom, caso.geom.aletas
    env = caso.envolvente
    h_env, w_env = envolvente(a.r_tip, g.R, a.params.n, a.params.rotacion)
    h_g, w_g = envolvente(a.r_tip, g.R, a.params.n, env.rot_guardado)

    def cabe(v, lim):
        return None if lim is None else bool(v <= lim + 1e-12)

    return {"L_total_mm": g.L_total / MM, "D_acostado_mm": max(h_g, w_g) / MM,
            "D_parado_mm": 2 * max(g.R, a.r_tip) / MM, "h_env_mm": h_env / MM, "w_env_mm": w_env / MM,
            "cabe_largo": cabe(g.L_total, env.largo_max), "cabe_alto": cabe(h_g, env.alto_max),
            "cabe_ancho": cabe(w_g, env.ancho_max)}


def tolerancia_amarre(m: float, x_CG: float, x_CP: float, q: float, S_ref: float, CNa: float,
                      alpha_max: float) -> float:
    """Error admisible |x_CG − x_T| para que |α_trim| ≤ α_max con el amarre cerca del CG.

    De α ≈ m g (x_CG − x_T) / (q S C_Nα (x_CP − x_T)) con x_T ≈ x_CG en el brazo:
        |x_CG − x_T| ≤ α_max q S C_Nα (x_CP − x_CG) / (m g)
    """
    brazo = x_CP - x_CG
    return alpha_max * q * S_ref * CNa * brazo / (m * G0) if brazo > 0 else math.nan


def _txt(banderas: list[str]) -> str:
    return ";".join(dict.fromkeys(banderas))


def analizar_caso(caso: Caso, rellenos: list[Relleno] | None = None,
                  perfil: pd.DataFrame | None = None, aero_or: AeroOR | None = None,
                  masas_or: dict[str, float] | None = None,
                  tol: Tolerancias = Tolerancias()) -> ResultadoCaso:
    """perfil: filas de or_perfiles.csv de este caso (None → perfil analítico).
    masas_or: {componente: m [kg]} de OpenRocket (nariz, cuerpo, cola, aletas)."""
    g, nu = caso.geom, caso.numerico
    x = malla(g.L, nu.dx)
    r_e = None
    if perfil is not None and len(perfil):
        r_e = perfil_desde_muestras(x, perfil["x_mm"].to_numpy() * MM, perfil["r_ext_mm"].to_numpy() * MM)
    cav = construir_cavidad(caso, r_e=r_e, x=x)
    mv = masa_vacia(caso, cav)
    cp_int = cp_interno(g, cav.x, cav.r_e, caso.vuelo.mach, nu.n_franjas_aleta)

    banderas: list[str] = []
    if aero_or is not None and math.isfinite(aero_or.x_CP) and math.isfinite(aero_or.CNa):
        x_CP, CNa, fuente = aero_or.x_CP, aero_or.CNa, "openrocket"
        if tol.dif_CP_frac_L is not None and abs(x_CP - cp_int.x_CP) > tol.dif_CP_frac_L * g.L:
            banderas.append("dif_CP_alta")
    else:
        x_CP, CNa, fuente = cp_int.x_CP, cp_int.CNa, "interno"

    dif_masa, dif_total = {}, NAN
    if masas_or:
        propias = masas_por_componente(caso, mv)
        for comp, m_or in masas_or.items():
            if comp in propias and m_or > 0:
                dif_masa[comp] = propias[comp][0] / m_or - 1
        comunes = [c for c in masas_or if c in propias]
        m_or_tot = sum(masas_or[c] for c in comunes)
        if m_or_tot > 0:
            dif_total = sum(propias[c][0] for c in comunes) / m_or_tot - 1
        for comp, d in dif_masa.items():
            lim = tol.dif_masa_aletas_frac if comp == "aletas" else tol.dif_masa_frac
            if lim is not None and abs(d) > lim:
                banderas.append("dif_masa_alta")

    x_b0 = x_inicio_lastre(caso, cav)
    lim = limite_geometrico(caso, cav, x_b0)
    if not lim.ell_geo > 0:
        banderas.append("inviable_geo")

    res = ResultadoCaso(caso=caso, cav=cav, mv=mv, cp_int=cp_int, x_CP=x_CP, CNa=CNa,
                        x_CP_fuente=fuente, x_b0=x_b0, lim=lim, banderas=banderas,
                        dif_masa=dif_masa, dif_masa_total=dif_total)
    for rel in (rellenos if rellenos is not None else list(caso.rellenos)):
        res.rellenos.append(_analizar_relleno(res, rel))
    return res


def _analizar_relleno(res: ResultadoCaso, rel: Relleno) -> ResultadoRelleno:
    caso, cav, mv = res.caso, res.cav, res.mv
    g, nu, vu, es = caso.geom, caso.numerico, caso.vuelo, caso.estabilidad
    D_ref, S_ref = g.D, math.pi * g.D**2 / 4
    x_CP, CNa = res.x_CP, res.CNa
    x_b0 = res.x_b0 if math.isfinite(res.x_b0) else 0.0
    mod = ModeloLastre(caso=caso, cav=cav, mv=mv, x_b0=x_b0, rho_b=rel.rho_b)
    ell_geo = res.lim.ell_geo
    banderas = list(res.banderas)

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
        "costo_usd_kg": rel.costo_usd_kg, "L_mm": g.L / MM, "D_mm": g.D / MM,
        "V_ext_cm3": V_ext / CM3, "V_int_cm3": V_int / CM3,
        "x_b0_mm": res.x_b0 / MM, "ell_geo_mm": ell_geo / MM,
        "x_CP_mm": x_CP / MM, "x_CP_fuente": res.x_CP_fuente,
        "dx_CP_or_vs_interno_mm": (x_CP - res.cp_int.x_CP) / MM if res.x_CP_fuente == "openrocket" else NAN,
        "CN_alpha_rad": CNa,
        "m_casco_g": mv.m_casco / G, "A_aleta_mm2": g.aletas.area / MM2, "m_aletas_g": mv.m_aletas / G,
        "x_CG_aletas_mm": mv.x_aletas / MM, "m_mamparos_g": mv.m_mamparos / G,
        "m_electronica_g": caso.electronica.me / G, "m_puntuales_g": mv.m_puntuales / G,
        "dif_masa_or_pct": 100 * res.dif_masa_total,
        **banderas_envolvente(caso),
    }

    ventana = None
    curva = pd.DataFrame(columns=CURVA)
    if ell_geo > 0:
        Vg = float(mod.V_b(ell_geo))
        xcg_g = float(mod.x_CG(ell_geo))
        fila.update({
            "V_util_geo_cm3": Vg / CM3, "eta_int_geo": Vg / V_int, "eta_ext_geo": Vg / V_ext,
            "m_relleno_geo_g": rel.rho_b * Vg / G, "m_total_geo_g": float(mod.m(ell_geo)) / G,
            "x_CG_geo_mm": xcg_g / MM, "SM_geo_cal": SM(xcg_g),
            "alpha_trim_geo_deg": float(trim_deg(ell_geo)),
        })
        ventana = ventana_SM(mod, ell_geo, x_CP, D_ref, es.SM_min, es.SM_max, nu.n_ell, nu.tol)
        SM_inf = float(mod.SM_inf(ventana.ell_star, x_CP, D_ref)) if ventana.ell_star > 0 else NAN
        fila.update({
            "ell_star_mm": ventana.ell_star / MM, "x_CG_min_mm": ventana.x_CG_min / MM,
            "SM_max_alcanzable_cal": SM(ventana.x_CG_min),
            "SM_inf_cal": SM_inf, "margen_SM_inf_cal": SM_inf - es.SM_min,
        })
        por_geometria = math.isfinite(SM_inf) and SM_inf < es.SM_min
        if por_geometria:
            # no se proponen longitudes de lastre: el techo lo ponen las aletas y el CP
            banderas.append("SM_inalcanzable_por_geometria")
            ventana.intervalos, ventana.ell_SM_min, ventana.ell_SM_max = [], NAN, NAN
            if "SM_inalcanzable" not in ventana.banderas:
                ventana.banderas.append("SM_inalcanzable")
        elif math.isfinite(SM_inf) and SM_inf - es.SM_min < es.umbral_margen_bajo:
            banderas.append("margen_SM_bajo")
        banderas += ventana.banderas
        fila.update({"ell_SM_min_mm": ventana.ell_SM_min / MM, "ell_SM_max_mm": ventana.ell_SM_max / MM})
        if math.isfinite(ventana.ell_SM_min):
            # no contractual (no va a resultados.csv): alimenta la figura de masa necesaria
            fila["m_lastre_SM_min_g"] = float(mod.m_b(ventana.ell_SM_min)) / G
        ell_u = ventana.ell_SM_max
        if math.isfinite(ell_u):
            # Óptimo de masa con SM ≥ SM_min. Primero el tapón delantero (suma masa y adelanta el
            # CG); solo si llega a ℓ_geo, el tapón trasero detrás de la electrónica hasta SM = SM_min.
            ell2 = 0.0
            lleno = ell_u >= ell_geo - nu.tol
            if caso.lastre.trasero and lleno:
                ell2 = mod.ell2_por_SM(ell_u, x_CP, D_ref, es.SM_min, nu.tol)
            V1 = float(mod.V_b(ell_u))
            V2 = float(mod.V_tras(ell_u, ell2))
            Vu = V1 + V2
            m_u = float(mod.m_con_trasero(ell_u, ell2))
            xcg_u = float(mod.x_CG_con_trasero(ell_u, ell2))
            if not lleno:
                limitante = "SM"
            elif caso.lastre.trasero:
                limitante = "SM" if ell2 < mod.ell2_max(ell_u) - nu.tol else "geometria"
            else:
                limitante = "geometria"
            alpha_u = float(np.degrees(alpha_trim(m_u, xcg_u, x_T_de(xcg_u), x_CP, vu.q, S_ref, CNa)))
            if not (abs(alpha_u) <= np.degrees(caso.remolque.alpha_lineal)):  # también NaN
                banderas.append("trim_no_lineal")
            fila.update({
                "V_util_cm3": Vu / CM3, "eta_int": Vu / V_int, "eta_ext": Vu / V_ext,
                "m_relleno_g": rel.rho_b * Vu / G, "m_total_g": m_u / G,
                "x_CG_mm": xcg_u / MM, "SM_cal": SM(xcg_u), "x_T_mm": x_T_de(xcg_u) / MM,
                "alpha_trim_deg": alpha_u,
                "x_T_trim_cero_mm": xcg_u / MM,
                "tol_amarre_mm": tolerancia_amarre(m_u, xcg_u, x_CP, vu.q, S_ref, CNa,
                                                   caso.remolque.alpha_max) / MM,
                "ell_trasero_mm": ell2 / MM,
                "costo_relleno_usd": rel.rho_b * Vu * rel.costo_usd_kg,
                "m_relleno_trasero_g": rel.rho_b * V2 / G, "limitante_masa": limitante,
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
    m_del_max = float(mod.m_b(ell_geo)) if ell_geo > 0 else 0.0
    for m_obj in caso.presupuesto:
        ell, ell2 = mod.ell_de_masa(m_obj, nu.tol), 0.0
        if caso.lastre.trasero and ell_geo > 0 and m_obj > m_del_max:
            # el delantero se llena hasta ℓ_geo y el resto va detrás de la electrónica
            resto = m_obj - m_del_max
            l2max = mod.ell2_max(ell_geo)
            if rel.rho_b * float(mod.V_tras(ell_geo, l2max)) >= resto:
                ell = ell_geo
                ell2 = brentq(lambda l2: rel.rho_b * float(mod.V_tras(ell_geo, l2)) - resto, 0.0, l2max,
                              xtol=nu.tol)
            else:
                ell = math.nan
        no_cabe = (not math.isfinite(ell)) or ell > ell_geo + nu.tol
        p = {"config_id": caso.id, "relleno": rel.nombre, "m_lastre_obj_g": m_obj / G,
             "ell_mm": ell / MM, "ell_trasero_mm": ell2 / MM, "no_cabe": bool(no_cabe)}
        if math.isfinite(ell):
            m_p = float(mod.m_con_trasero(ell, ell2))
            xcg = float(mod.x_CG_con_trasero(ell, ell2))
            p.update({"m_total_g": m_p / G, "x_CG_mm": xcg / MM, "SM_cal": SM(xcg),
                      "costo_usd": m_obj * rel.costo_usd_kg,
                      "alpha_trim_deg": float(np.degrees(alpha_trim(m_p, xcg, x_T_de(xcg), x_CP, vu.q,
                                                                    S_ref, CNa)))})
        presupuestos.append(p)

    fila["banderas"] = _txt(banderas)
    return ResultadoRelleno(relleno=rel, modelo=mod, ventana=ventana, fila=fila, curva=curva,
                            presupuestos=presupuestos)
