"""Verificación y calibración con OpenRocket (v3 §4.6–4.8). Único módulo de sensor_opt con Java.

`PuenteOpt` extiende `sensor_lastre.or_bridge.PuenteOR` sin modificarlo: el cuerpo y las paredes
se fijan con `PuenteOR.aplicar` (con la aleta marcador del cuerpo) y luego se escribe la aleta
contenida (4 puntos, posición *bottom* con δ_b hacia proa), el material Onyx y el lastre, la
electrónica y las masas puntuales como *Mass components* (API en docs/api_openrocket.md).

Flujo (`verificar_y_calibrar`):
1. se verifican los N_verif mejores y una muestra de N_cal estratificada por (D, k);
2. se ajusta x_OR ≈ α + β x_sust; si el residuo máximo supera tol_cp_mm se recalcula toda la
   malla con el CP calibrado y se verifican los nuevos mejores, hasta max_iter_calibracion o
   hasta que los N_verif mejores no cambien;
3. cada verificado se vuelve a llenar con el CP de OpenRocket (SM_or, J_or); gana el mejor J_or
   con SM_or en [SM_min, SM_max].
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from sensor_lastre.or_bridge import ErrorAleta, PuenteOR

from .barrido import completar, evaluar, grupos_desde_ranking
from .config import ConfigOpt
from .exportar import Detalle, detalle, specs_de_fila
from .objetivo import J, f_valores
from .sustituto import IDENTIDAD, Calibracion, ajustar

MM, G = 1e-3, 1e-3
L_MASA = 1e-4  # m, largo de los Mass components (masas casi puntuales en su CG)


@dataclass
class MasaOR:
    nombre: str
    m: float  # kg
    x: float  # CG global, m


def masas_del_llenado(det: Detalle) -> list[MasaOR]:
    """Lastre delantero y trasero, electrónica y masas puntuales del llenado, con su CG."""
    Ll = det.llenado
    rr, caso = Ll.rr, Ll.res.caso
    mod = rr.modelo
    out = []
    x_cg = rr.fila["x_CG_mm"] * MM
    if Ll.m_delantero > 0 and Ll.ell > 0:
        # el alojamiento del herraje se descuenta del tapón delantero en el CG (allí está el agujero)
        M_del = mod.rho_b * float(mod.Phi1(Ll.ell)) - mod.m_aloj * x_cg
        out.append(MasaOR("lastre_delantero", Ll.m_delantero, M_del / Ll.m_delantero))
    ell2 = (rr.fila.get("ell_trasero_mm") or 0.0) * MM
    if ell2 > 0:
        a = mod.x_r0(Ll.ell)
        V2 = float(mod.fll * mod.cav.vol(a, a + ell2))
        M2 = float(mod.fll * mod.cav.mom(a, a + ell2))
        out.append(MasaOR("lastre_trasero", mod.rho_b * V2, M2 / V2))
    el = caso.electronica
    out.append(MasaOR("electronica", el.me, float(mod.xbar_e(Ll.ell))))
    for p in caso.puntuales:
        out.append(MasaOR(p.nombre, p.m, x_cg if p.en_CG else p.x))
    return out


class PuenteOpt(PuenteOR):
    """PuenteOR + aleta contenida, material Onyx y Mass components."""

    def __init__(self, inst, nombres: dict[str, str], ork_base: str | Path):
        super().__init__(inst, nombres)
        self.ork_base = str(ork_base)

    def aplicar_candidato(self, cfg: ConfigOpt, det: Detalle) -> dict:
        """Recarga la base, fija el cuerpo y escribe la aleta contenida. Devuelve r_LE, x_LE, puntos."""
        self.cargar(self.ork_base)
        cu, g = det.cu, det.llenado.res.caso.geom.aletas
        self.aplicar(cu.caso)  # cuerpo, paredes y una aleta marcador válida
        f, t = self.comp["aletas"], self.comp["cola"]
        pa = g.params
        Lt = float(t.getLength())
        c_t = float(g.puntos[2, 0] - g.puntos[1, 0])
        x_LE_loc = Lt - pa.delta_b - pa.c_r
        r_LE = float(t.getRadius(x_LE_loc))
        r_TE = float(t.getRadius(Lt - pa.delta_b))
        h = float(g.r_tip - r_LE)
        pts = [(0.0, 0.0), (pa.x_s, h), (pa.x_s + c_t, h), (pa.c_r, r_TE - r_LE)]
        AM = self.core.rocketcomponent.position.AxialMethod
        f.setAxialMethod(AM.BOTTOM)
        f.setAxialOffset(-float(pa.delta_b))  # OpenRocket: offset positivo = hacia popa
        C = self.core.util.Coordinate
        f.setPoints(self._jpype.JArray(C)([C(float(x), float(y)) for x, y in pts]))
        f.setFinCount(int(pa.n))
        f.setBaseRotation(float(pa.rotacion))
        f.setThickness(float(pa.t))
        f.setMaterial(self._material(cfg.aleta.material.upper(), pa.rho * pa.phi))
        aceptados = [(float(p.x), float(p.y)) for p in f.getFinPoints()]
        if len(aceptados) != len(pts) or np.max(np.abs(np.array(aceptados) - np.array(pts))) > 1e-7:
            raise ErrorAleta(f"OpenRocket modificó los puntos de la aleta contenida: {pts} → {aceptados}")
        return {"r_LE": r_LE, "x_LE": self.x_abs("aletas"), "puntos": aceptados}

    def agregar_masas(self, masas: list[MasaOR]):
        MC = self.core.rocketcomponent.MassComponent
        AM = self.core.rocketcomponent.position.AxialMethod
        padre = self.comp["cuerpo"]
        for m in masas:
            mc = MC(L_MASA, 0.001, float(m.m))
            mc.setName(m.nombre)
            padre.addChild(mc)
            mc.setAxialMethod(AM.ABSOLUTE)
            mc.setAxialOffset(float(m.x - L_MASA / 2))

    def guardar(self, ruta: Path):
        ruta.parent.mkdir(parents=True, exist_ok=True)
        self.orl.save_doc(str(ruta), self.doc)


# --------------------------------------------------------------------------- verificación de un candidato


def verificar_candidato(cfg: ConfigOpt, pr: PuenteOpt, fila: pd.Series, cal: Calibracion,
                        guardar: Path | None = None) -> dict:
    """CP, masas y CG de OpenRocket para un candidato; SM y J con el CP de OpenRocket."""
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a, cal)
    out = {"cand_id": fila["cand_id"], "x_CP_sust_mm": fila.get("x_CP_sust_mm"),
           "x_CP_cal_mm": float(cal.aplicar(fila["x_CP_sust_mm"] * MM)) / MM}
    if det.llenado is None:
        out.update({"factible_or": False, "advertencias": ";".join(det.motivos)})
        return out
    pr.aplicar_candidato(cfg, det)
    aero = pr.aero(cfg_mach(det))
    x_or = aero.x_CP
    # llenado final con el CP de OpenRocket (v3 §4.8)
    det_or = detalle(cfg, c, a, x_CP=x_or)
    Ll = det_or.llenado
    out.update({"x_CP_or_mm": x_or / MM, "dx_or_sust_mm": x_or / MM - out["x_CP_sust_mm"],
                "dx_or_cal_mm": x_or / MM - out["x_CP_cal_mm"], "CNa_or": aero.CNa,
                "CNa_sust": fila.get("CN_alpha_total"), "advertencias": ";".join(aero.warnings)})
    for comp, (cna, x) in aero.por_componente.items():
        out[f"CNa_or_{comp}"], out[f"x_CP_or_{comp}_mm"] = cna, x / MM
    D = det.cu.spec.D
    if Ll is None or not math.isfinite(Ll.ell):
        out.update({"factible_or": False, "SM_or_cal": math.nan})
        out["advertencias"] = ";".join(filter(None, [out["advertencias"]] + (det_or.motivos or [])))
        return out
    f = Ll.fila
    pr.agregar_masas(masas_del_llenado(det_or))
    mas = pr.masas()
    out.update({"m_or_g": mas.m_total / G, "x_CG_or_mm": mas.x_CG_total / MM, "m_modelo_g": f["m_total_g"],
                "x_CG_modelo_mm": f["x_CG_mm"], "dif_masa_or_pct": 100 * (mas.m_total / G / f["m_total_g"] - 1),
                "SM_or_cal": (x_or - f["x_CG_mm"] * MM) / D, "m_total_or_g": f["m_total_g"],
                "restriccion_activa_or": Ll.restriccion_activa, "factible_or": Ll.ok,
                "motivos_or": ";".join(Ll.motivos)})
    for comp, (m, _) in mas.por_componente.items():
        out[f"m_or_{comp}_g"] = m / G
    F = f_valores(cfg, [f["m_total_g"] * G], [fila["D_ap_mm"] * MM], [fila["k_efectivo"]], [out["SM_or_cal"]])
    out["J_or"] = float(J(cfg, F)[0])
    if guardar is not None:
        pr.guardar(guardar)
    return out


def cfg_mach(det: Detalle) -> float:
    return det.cu.caso.vuelo.mach


# --------------------------------------------------------------------------- muestra y calibración


def muestra_estratificada(df: pd.DataFrame, n: int, excluir: set[str], semilla: int) -> list[str]:
    """n candidatos factibles repartidos por estratos (D, k), en ronda, con semilla fija."""
    rng = np.random.default_rng(semilla)
    fac = df[df["factible"].astype(bool) & ~df["cand_id"].isin(excluir)]
    estratos = [list(g["cand_id"]) for _, g in fac.groupby(["D_mm", "k"], sort=True)]
    for e in estratos:
        rng.shuffle(e)
    out: list[str] = []
    while len(out) < n and any(estratos):
        for e in estratos:
            if e and len(out) < n:
                out.append(e.pop())
    return out


@dataclass
class ResultadoVerificacion:
    ranking: pd.DataFrame
    verificacion: pd.DataFrame
    calibracion: Calibracion
    historial: list[dict] = field(default_factory=list)
    ganador: str | None = None


def verificar_y_calibrar(cfg: ConfigOpt, pr: PuenteOpt, ranking: pd.DataFrame, procesos: int | None = None,
                         log=print) -> ResultadoVerificacion:
    vo = cfg.ejecucion.get("verificacion_or") or {}
    N_verif, N_cal = int(vo.get("N_verif", 20)), int(vo.get("N_cal", 30))
    tol = float(vo.get("tol_cp_mm", 2.0)) * MM
    max_iter = int(vo.get("max_iter_calibracion", 3))
    semilla = int(cfg.ejecucion.get("semilla", 12345))
    cal = IDENTIDAD
    df = ranking
    hechos: dict[str, dict] = {}
    historial = []

    def verificar(ids: list[str], rol: str, it: int):
        filas = df.set_index("cand_id")
        for cid in ids:
            if cid in hechos:
                hechos[cid]["rol"] += f"+{rol}" if rol not in hechos[cid]["rol"] else ""
                continue
            fila = filas.loc[cid].copy()
            fila["cand_id"] = cid
            try:
                r = verificar_candidato(cfg, pr, fila, IDENTIDAD)
            except ErrorAleta as e:
                r = {"cand_id": cid, "factible_or": False, "advertencias": str(e)}
            r.update({"rol": rol, "iteracion": it})
            hechos[cid] = r
            log(f"    {rol:13s} {cid}  x_sust {r.get('x_CP_sust_mm', math.nan):7.2f}  "
                f"x_OR {r.get('x_CP_or_mm', math.nan):7.2f} mm  SM_OR {r.get('SM_or_cal', math.nan):.2f}")

    top = list(df[df["factible"].astype(bool)].head(N_verif)["cand_id"])
    log(f"  Iteración 0: {len(top)} mejores + muestra de calibración de {N_cal}")
    verificar(top, "verificacion", 0)
    verificar(muestra_estratificada(df, N_cal, set(top), semilla), "calibracion", 0)
    for it in range(1, max_iter + 1):
        v = pd.DataFrame(hechos.values())
        v = v[np.isfinite(v.get("x_CP_or_mm", pd.Series(dtype=float)).astype(float))]
        cal = ajustar(v["x_CP_sust_mm"] * MM, v["x_CP_or_mm"] * MM, iteracion=it)
        historial.append(cal.a_dict())
        log(f"  Calibración {it}: α = {cal.alpha / MM:+.3f} mm, β = {cal.beta:.5f}, R² = {cal.R2:.5f}, "
            f"residuo máx. {cal.res_max / MM:.3f} mm")
        if cal.res_max <= tol:
            log(f"  Residuo ≤ {tol / MM:g} mm: el sustituto basta, no se recalcula la malla.")
            break
        df = _recalcular(cfg, df, cal, procesos, log)
        nuevo_top = list(df[df["factible"].astype(bool)].head(N_verif)["cand_id"])
        if set(nuevo_top) <= set(hechos):
            log("  Los mejores no cambiaron: fin de la calibración.")
            break
        verificar([c for c in nuevo_top if c not in hechos], "verificacion", it)

    ver = pd.DataFrame(hechos.values())
    ver["x_CP_cal_mm"] = cal.aplicar(ver["x_CP_sust_mm"] * MM) / MM
    ver["dx_or_cal_mm"] = ver["x_CP_or_mm"] - ver["x_CP_cal_mm"]
    df = _anotar(df, ver)
    ok = ver[ver["factible_or"].fillna(False).astype(bool)] if "factible_or" in ver else ver.iloc[0:0]
    ganador = None
    if not ok.empty:
        ganador = str(ok.sort_values(["J_or", "cand_id"]).iloc[0]["cand_id"])
    return ResultadoVerificacion(df, ver, cal, historial, ganador)


def _recalcular(cfg, df, cal, procesos, log):
    log(f"  Recalculando {len(df):,} candidatos con el CP calibrado…")
    nuevo = evaluar(cfg, grupos_desde_ranking(df), cal, procesos)
    if "refinamiento" in df:
        nuevo = nuevo.merge(df[["cand_id", "refinamiento"]], on="cand_id", how="left")
    return nuevo


def _anotar(df: pd.DataFrame, ver: pd.DataFrame) -> pd.DataFrame:
    cols = {"x_CP_or_mm": "x_CP_or_mm", "SM_or_cal": "SM_or_cal", "dif_masa_or_pct": "dif_masa_or_pct",
            "J_or": "J_or", "m_total_or_g": "m_total_or_g", "x_CG_or_mm": "x_CG_or_mm"}
    df = df.drop(columns=[c for c in cols.values() if c in df]).merge(
        ver[["cand_id"] + [c for c in cols if c in ver]], on="cand_id", how="left")
    df["verificado_or"] = df["cand_id"].isin(ver["cand_id"])
    return df


# --------------------------------------------------------------------------- perfiles (T7)


def comparar_perfiles(cfg: ConfigOpt, pr: PuenteOpt, n: int = 400) -> pd.DataFrame:
    """Máx |r_OR − r_analítico| de la cola para cada forma de la malla, con k y L_t extremos."""
    from .geometria import construir_cuerpo
    from sensor_lastre.perfiles import radio_exterior
    filas = []
    m = cfg.malla
    for forma in m["cola_forma"]:
        for k in (min(m["k"]), max(m["k"])):
            for lt in (min(m["L_t_rel_D"]), max(m["L_t_rel_D"])):
                from .config import CuerpoSpec
                spec = CuerpoSpec(float(max(m["D_mm"])), forma["forma"],
                                  None if forma.get("parametro") is None else float(forma["parametro"]), float(lt), float(k))
                cu = construir_cuerpo(cfg, spec)
                if cu.caso is None:
                    continue
                pr.cargar(pr.ork_base)
                pr.aplicar(cu.caso)
                per = [(x, r) for comp, x, r in pr.perfil(n) if comp == "cola"]
                x = np.array([p[0] for p in per])
                r_or = np.array([p[1] for p in per])
                r_an = radio_exterior(x, cu.caso.geom)
                filas.append({"cuerpo_id": spec.id, "dif_max_mm": float(np.abs(r_or - r_an).max() / MM)})
    return pd.DataFrame(filas)
