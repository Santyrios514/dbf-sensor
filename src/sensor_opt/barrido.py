"""Barrido en malla con caché por cuerpo, en paralelo; ranking, Pareto y refinamiento (v3 §4).

Unidad de trabajo: un cuerpo con su lista de aletas. El cuerpo (perfil, cavidad, masas del
casco, x_b0, ℓ_geo, CP de nariz y cola) se construye una vez; cada aleta solo agrega su
polígono, su masa, su C_Nα y el llenado. Los resultados se devuelven en el mismo orden en que se
pidieron, así que el ranking es reproducible (v3 T9) sin importar el número de procesos.
"""

from __future__ import annotations

import itertools
import math
import multiprocessing as mp
from collections import OrderedDict

import numpy as np
import pandas as pd

from .aletas import construir, velocidad_flutter
from .config import VARS_ALETA, AletaSpec, ConfigOpt, CuerpoSpec, cand_id
from .esquemas import RANKING
from .geometria import construir_cuerpo
from .lastre import llenar
from .objetivo import J, f_valores, pareto
from .sustituto import IDENTIDAD, Calibracion, cp_sustituto

MM, G = 1e-3, 1e-3
Grupos = "OrderedDict[CuerpoSpec, list[AletaSpec]]"


# --------------------------------------------------------------------------- evaluación


def _fila_base(cfg: ConfigOpt, c: CuerpoSpec, a: AletaSpec) -> dict:
    return {
        "cand_id": cand_id(c, a), "cuerpo_id": c.id, "factible": False, "motivos": "",
        "D_mm": c.D_mm, "cola_forma": c.forma, "cola_parametro": c.parametro,
        "L_t_rel_D": c.L_t_rel_D, "L_t_mm": c.L_t / MM, "k": c.k,
        "lambda_LE": a.lambda_LE, "mu_cr": a.mu_cr, "r_tip_rel_R": a.r_tip_rel_R, "gamma_ct": a.gamma_ct,
        "sigma_flecha": a.sigma_flecha, "verificado_or": False, "en_pareto": False,
    }


def evaluar_cuerpo(cfg: ConfigOpt, c: CuerpoSpec, aletas: list[AletaSpec],
                   cal: Calibracion = IDENTIDAD) -> list[dict]:
    """Filas del ranking para todas las aletas de un cuerpo (factibles o no, con motivos)."""
    cu = construir_cuerpo(cfg, c)
    filas = []
    if cu.caso is None or not cu.ok:
        for a in aletas:
            f = _fila_base(cfg, c, a)
            f.update({"motivos": ";".join(cu.motivos), "theta_eq_deg": math.degrees(cu.theta_eq)})
            filas.append(f)
        return filas
    rest, pa = cfg.restricciones, cfg.aleta
    vuelo = cu.caso.vuelo
    n_franjas = cu.caso.numerico.n_franjas_aleta
    R = c.D / 2
    comunes = {"theta_eq_deg": math.degrees(cu.theta_eq), "L_c_mm": cu.L_c / MM, "x_b0_mm": cu.x_b0 / MM,
               "ell_geo_mm": cu.lim.ell_geo / MM, "m_casco_g": cu.mv.m_casco / G}
    for a in aletas:
        f = _fila_base(cfg, c, a)
        f.update(comunes)
        g, d, motivos = construir(cfg, cu, a)
        f.update({"x_LE_mm": d.x_LE / MM, "c_r_mm": d.c_r / MM, "c_t_mm": d.c_t / MM, "x_s_mm": d.x_s / MM,
                  "r_LE_mm": d.r_LE / MM, "r_TE_raiz_mm": d.r_TE_raiz / MM, "r_tip_mm": d.r_tip / MM,
                  "h_mm": d.h / MM, "D_ap_mm": 2 * max(R, d.r_tip) / MM,
                  "h_45_mm": max(c.D, math.sqrt(2) * d.r_tip) / MM})
        if g is None:
            f["motivos"] = ";".join(motivos)
            filas.append(f)
            continue
        V_f = velocidad_flutter(d.h, g.area, d.c_r, d.c_t, pa.t, pa.E, pa.nu, vuelo.altitud)
        f.update({"A_aleta_mm2": g.area / MM**2, "m_aletas_g": g.masa / G, "V_flutter_m_s": V_f,
                  "flutter_margen_bajo": bool(V_f < pa.factor_flutter * vuelo.V)})
        cp, motivo = cp_sustituto(cu, g, vuelo.mach, n_franjas, rest.eps_CN)
        if cp is None:
            f["motivos"] = motivo
            filas.append(f)
            continue
        x_sust = cp.x_CP
        x_cal = float(cal.aplicar(x_sust))
        cp_uso = cp if cal.identidad else type(cp)(x_CP=x_cal, CNa=cp.CNa, partes=cp.partes)
        Ll = llenar(cu, g, cp_uso, rest.SM_max, "sustituto" if cal.identidad else "calibrado")
        r = Ll.fila
        f.update({"x_CP_sust_mm": x_sust / MM, "x_CP_cal_mm": x_cal / MM, "CN_alpha_total": cp.CNa,
                  "restriccion_activa": Ll.restriccion_activa, "m_electronica_g": r.get("m_electronica_g")})
        if math.isfinite(Ll.ell):
            f.update({"m_total_g": r["m_total_g"], "m_lastre_g": r["m_relleno_g"],
                      "m_lastre_delantero_g": Ll.m_delantero / G,
                      "m_lastre_trasero_g": r.get("m_relleno_trasero_g", 0.0),
                      "ell_mm": Ll.ell / MM, "ell_trasero_mm": r.get("ell_trasero_mm", 0.0),
                      "x_CG_mm": r["x_CG_mm"], "SM_cal": r["SM_cal"], "alpha_trim_deg": r.get("alpha_trim_deg"),
                      "x_T_mm": r.get("x_T_mm"), "x_T_trim_cero_mm": r.get("x_T_trim_cero_mm"),
                      "tol_amarre_mm": r.get("tol_amarre_mm"), "banderas_lastre": r.get("banderas", "")})
        f["motivos"] = ";".join(Ll.motivos)
        f["factible"] = Ll.ok
        filas.append(f)
    return filas


_CFG: ConfigOpt | None = None
_CAL: Calibracion = IDENTIDAD


def _init(cfg: ConfigOpt, cal: Calibracion):
    global _CFG, _CAL
    _CFG, _CAL = cfg, cal


def _tarea(args):
    c, aletas = args
    return evaluar_cuerpo(_CFG, c, aletas, _CAL)


def evaluar(cfg: ConfigOpt, grupos: Grupos, cal: Calibracion = IDENTIDAD, procesos: int | None = None,
            progreso=None) -> pd.DataFrame:
    """Evalúa {cuerpo: [aletas]} y devuelve las filas en el mismo orden, con f, J y el ranking."""
    procesos = procesos or cfg.procesos
    tareas = list(grupos.items())
    filas: list[dict] = []
    if procesos <= 1 or len(tareas) <= 1:
        _init(cfg, cal)
        for i, t in enumerate(tareas):
            filas += _tarea(t)
            if progreso:
                progreso(i + 1, len(tareas))
    else:
        with mp.get_context("fork").Pool(procesos, initializer=_init, initargs=(cfg, cal)) as pool:
            for i, r in enumerate(pool.imap(_tarea, tareas, chunksize=1)):
                filas += r
                if progreso:
                    progreso(i + 1, len(tareas))
    return completar(cfg, pd.DataFrame(filas))


def completar(cfg: ConfigOpt, df: pd.DataFrame) -> pd.DataFrame:
    """Objetivos, J, puesto y frente de Pareto; columnas en el orden del esquema."""
    df = df.reindex(columns=list(dict.fromkeys(RANKING + list(df.columns))))
    ok = df["factible"].astype(bool)
    m_total = df["m_total_g"].astype(float) * G
    SM = df["SM_cal"].astype(float)
    F = f_valores(cfg, m_total, df["D_ap_mm"].astype(float) * MM, df["k"].astype(float), SM)
    valido = np.isfinite(m_total) & np.isfinite(SM)
    for j, col in enumerate(("f1", "f2", "f3", "f4")):
        df[col] = np.where(valido, F[:, j], np.nan)
    df["J"] = np.where(valido, J(cfg, np.nan_to_num(F)), np.nan)
    df["en_pareto"] = False
    df = ordenar(df)
    ok = df["factible"].astype(bool)
    if ok.any():
        fac = df[ok]
        tol = float(cfg.numerico.get("tol_masa_max_g", 0.5)) * G
        alcanza = fac["m_total_g"] * G >= cfg.m_max - tol
        base = fac[alcanza] if alcanza.any() else fac[fac["f1"] <= fac["f1"].min() + tol / cfg.m_max]
        nd = pareto(base[["f2", "f3", "f4"]].to_numpy())
        df.loc[base.index[nd], "en_pareto"] = True
    return df


def ordenar(df: pd.DataFrame) -> pd.DataFrame:
    """Factibles por J (desempate por cand_id), luego los infactibles; puesto = 1..n factibles."""
    df = df.assign(_fac=~df["factible"].astype(bool)).sort_values(
        ["_fac", "J", "cand_id"], na_position="last", kind="stable").drop(columns="_fac")
    df = df.reset_index(drop=True)
    n_ok = int(df["factible"].astype(bool).sum())
    df["puesto"] = np.where(np.arange(len(df)) < n_ok, np.arange(1, len(df) + 1), np.nan)
    return df


# --------------------------------------------------------------------------- malla y refinamiento


def grupos_malla(cfg: ConfigOpt, malla: dict | None = None) -> Grupos:
    aletas = cfg.aletas(malla)
    return OrderedDict((c, list(aletas)) for c in cfg.cuerpos(malla))


def grupos_desde_ranking(df: pd.DataFrame) -> Grupos:
    """Reconstruye {cuerpo: [aletas]} a partir de las columnas del ranking (para recalcular)."""
    out: OrderedDict = OrderedDict()
    for r in df.itertuples(index=False):
        par = None if pd.isna(r.cola_parametro) else float(r.cola_parametro)
        c = CuerpoSpec(float(r.D_mm), r.cola_forma, par, float(r.L_t_rel_D), float(r.k))
        out.setdefault(c, []).append(AletaSpec(float(r.lambda_LE), float(r.mu_cr), float(r.r_tip_rel_R),
                                                float(r.gamma_ct), float(r.sigma_flecha)))
    return out


def _vecinos(v: float, lista: list[float], lo: float, hi: float) -> list[float]:
    """v y los puntos a medio paso hacia sus vecinos de la malla original, dentro de [lo, hi]."""
    xs = sorted(set(float(x) for x in lista))
    out = {v}
    menores = [x for x in xs if x < v - 1e-12]
    mayores = [x for x in xs if x > v + 1e-12]
    if menores:
        out.add(v - (v - menores[-1]) / 2)
    if mayores:
        out.add(v + (mayores[0] - v) / 2)
    return sorted(round(x, 10) for x in out if lo - 1e-12 <= x <= hi + 1e-12)


def grupos_refinamiento(cfg: ConfigOpt, ranking: pd.DataFrame, top_K: int) -> Grupos:
    """Malla a medio paso alrededor de los top_K factibles (v3 §4.5), sin repetir candidatos."""
    m, rest = cfg.malla, cfg.restricciones
    fac = ranking[ranking["factible"].astype(bool)].head(top_K)
    vistos = set(ranking["cand_id"])
    limites = {
        "D_mm": (min(m["D_mm"]), max(m["D_mm"])), "L_t_rel_D": (min(m["L_t_rel_D"]), max(m["L_t_rel_D"])),
        "k": (max(min(m["k"]), rest.k_min), max(m["k"])),
        "lambda_LE": (min(m["lambda_LE"]), min(max(m["lambda_LE"]), 1 - 1e-9)),
        "mu_cr": (min(m["mu_cr"]), min(max(m["mu_cr"]), 1.0)), "r_tip_rel_R": (min(m["r_tip_rel_R"]), max(m["r_tip_rel_R"])),
        "gamma_ct": (min(m["gamma_ct"]), min(max(m["gamma_ct"]), 1.0)),
        "sigma_flecha": (max(min(m["sigma_flecha"]), 0.0), min(max(m["sigma_flecha"]), 1.0)),
    }
    out: OrderedDict = OrderedDict()
    for r in fac.itertuples(index=False):
        par = None if pd.isna(r.cola_parametro) else float(r.cola_parametro)
        vals = {k: _vecinos(float(getattr(r, k)), m[k], *limites[k]) for k in limites}
        for D in vals["D_mm"]:
            for lt in vals["L_t_rel_D"]:
                for k in vals["k"]:
                    c = CuerpoSpec(D, r.cola_forma, par, lt, k)
                    lista = out.setdefault(c, [])
                    ya = {a.id for a in lista}
                    for combo in itertools.product(*(vals[v] for v in VARS_ALETA)):
                        a = AletaSpec(*combo)
                        if a.id in ya or cand_id(c, a) in vistos:
                            continue
                        ya.add(a.id)
                        lista.append(a)
    return OrderedDict((c, a) for c, a in out.items() if a)


def optimizar(cfg: ConfigOpt, refinar: bool = True, procesos: int | None = None,
              cal: Calibracion = IDENTIDAD, progreso=None) -> pd.DataFrame:
    """Barrido completo (+ refinamiento opcional) → ranking con todos los candidatos."""
    df = evaluar(cfg, grupos_malla(cfg), cal, procesos, progreso)
    ref = cfg.ejecucion.get("refinamiento") or {}
    if refinar and ref.get("activar", True):
        gr = grupos_refinamiento(cfg, df, int(ref.get("top_K", 10)))
        if gr:
            nuevo = evaluar(cfg, gr, cal, procesos, progreso)
            nuevo["refinamiento"] = True
            df = completar(cfg, pd.concat([df.assign(refinamiento=False), nuevo], ignore_index=True))
    else:
        df["refinamiento"] = False
    return df
