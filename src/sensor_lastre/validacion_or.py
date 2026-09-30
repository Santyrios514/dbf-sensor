"""Validación con OpenRocket de los ganadores del modelo propio.

El barrido (script 01) corre con el perfil analítico y el Barrowman interno, que es rápido y no
necesita JVM. Después (script 02) solo los mejores de cada grupo se aplican al `.ork` base y se
vuelven a analizar con el CP, el C_Nα, el perfil, el polígono de aleta y las masas de OpenRocket.

Selección iterativa por grupo (relleno × `validacion_or.agrupar_por`):
1. se validan los N_verif mejores del ranking interno;
2. el ganador provisional es el mejor del ranking con los valores de OpenRocket;
3. si algún candidato sin validar podría superarlo (más masa total, o la misma masa con más SM
   según el modelo propio), se validan los N_verif siguientes de esos; se repite hasta que
   ninguno pueda superarlo o se agote `max_iter`.

La exportación de un caso (`exportar_caso`) también la usa el modo `--todos` del script 02, que
reproduce la exportación completa de antes (or_*.csv para todas las configuraciones).
"""

from __future__ import annotations

import math
import traceback
from dataclasses import dataclass, field

import pandas as pd

from . import aletas as mod_aletas
from .aletas import envolvente
from .analisis import AeroOR, ResultadoCaso, Tolerancias, analizar_caso
from .config import Caso, Config
from .optimizacion import ranking

MM, G = 1e-3, 1e-3


# --------------------------------------------------------------------------- exportación de un caso


def fila_resumen(caso: Caso, meta: dict) -> dict:
    """Datos de entrada del caso para or_resumen.csv (sin resultados de OpenRocket)."""
    from .or_bridge import estaciones_T_rho
    g, pa, vu = caso.geom, caso.params_aleta, caso.vuelo
    T = estaciones_T_rho(caso)
    return {
        "config_id": caso.id, "L_mm": g.L / MM, "D_mm": g.D / MM, "L_cuerpo_mm": g.L / MM,
        "nariz_forma": g.nariz.forma, "nariz_parametro": g.nariz.parametro,
        "L_nariz_mm": g.nariz.Ln / MM, "L_cola_mm": g.cola.Lt / MM, "D_popa_mm": 2 * g.cola.Ra / MM,
        "cola_forma": g.cola.forma, "cola_parametro": g.cola.parametro, "popa": g.cola.popa,
        "aletas_n": pa.n, "aletas_t_mm": pa.t / MM, "aleta_c_r_mm": pa.c_r / MM,
        "aleta_x_tip_le_mm": pa.x_s / MM, "aleta_extension_mm": pa.e / MM,
        "aleta_r_interior_mm": pa.r_in / MM, "aleta_rotacion_deg": math.degrees(pa.rotacion),
        **{f"T_pared_{e}_mm": T[e][0] / MM for e in T}, **{f"rho_eq_{e}_kg_m3": T[e][1] for e in T},
        "mach": vu.mach, "mach_usado": vu.mach, "aoa_deg": math.degrees(vu.aoa),
        "S_ref_mm2": math.pi * g.D**2 / 4 / MM**2, "D_ref_mm": g.D / MM, **meta,
    }


@dataclass
class ExportOR:
    """Salida de OpenRocket de un caso, con el formato de los or_*.csv."""

    resumen: dict
    perfiles: list[dict] = field(default_factory=list)
    componentes: list[dict] = field(default_factory=list)
    aletas: list[dict] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.resumen.get("estado") == "ok"


def exportar_caso(pr, caso: Caso, ork, orc: dict, meta: dict, delta_aoa: float | None,
                  guardar=None, verbose: bool = False) -> ExportOR:
    """Aplica el caso al .ork base y lee CP, C_Nα, masas, perfil y polígono de aleta.
    Un fallo no lanza excepción: queda en resumen['estado'] y resumen['mensaje']."""
    from .or_bridge import ErrorAleta
    g, pa, vu = caso.geom, caso.params_aleta, caso.vuelo
    ex = ExportOR(resumen=fila_resumen(caso, meta))
    try:
        pr.cargar(ork)  # cada variante parte del .ork base
        geo = pr.aplicar(caso)
        aero = pr.aero(vu.mach, vu.aoa, delta_aoa)
        m = pr.masas()
        P, orig, x_LE, r_LE = pr.poligono_aleta()
        h_env, w_env = envolvente(r_LE + geo["h"], g.R, pa.n, pa.rotacion)
        ex.resumen.update({
            "estado": "ok", "mensaje": "", "L_total_mm": float(pr.configuracion.getLength()) / MM,
            "aletas_x_r0_mm": x_LE / MM, "aleta_h_tip_mm": geo["h"] / MM, "r_LE_mm": r_LE / MM,
            "r_tip_mm": (r_LE + geo["h"]) / MM, "x_CP_mm": aero.x_CP / MM, "CN_alpha_rad": aero.CNa,
            "CN_alpha_dif_rad": aero.CNa_dif, "CD": aero.CD, "m_vacio_or_g": m.m_total / G,
            "x_CG_vacio_or_mm": m.x_CG_total / MM, "h_env_mm": h_env / MM, "w_env_mm": w_env / MM,
            "or_warnings": " | ".join(aero.warnings),
        })
        ex.perfiles = [{"config_id": caso.id, "componente": k, "x_mm": x / MM, "r_ext_mm": r / MM}
                       for k, x, r in pr.perfil(int(orc.get("n_muestras_perfil", 800)))]
        for k in ("nariz", "cuerpo", "cola", "aletas"):
            cna, xcp = aero.por_componente.get(k, (math.nan, math.nan))
            mk, xk = m.por_componente[k]
            ex.componentes.append({"config_id": caso.id, "componente": k, "CN_alpha_rad": cna,
                                   "x_CP_mm": xcp / MM, "m_g": mk / G, "x_CG_mm": xk / MM})
        ex.componentes.append({"config_id": caso.id, "componente": "total", "CN_alpha_rad": aero.CNa,
                               "x_CP_mm": aero.x_CP / MM, "m_g": m.m_total / G, "x_CG_mm": m.x_CG_total / MM})
        ex.aletas = [{"config_id": caso.id, "orden": i, "x_mm": x / MM, "r_mm": r / MM, "origen": o}
                     for i, ((x, r), o) in enumerate(zip(P, orig))]
        if guardar is not None:
            guardar.parent.mkdir(parents=True, exist_ok=True)
            pr.orl.save_doc(str(guardar), pr.doc)
    except ErrorAleta as e:
        ex.resumen.update({"estado": "error_aleta", "mensaje": str(e)})
    except Exception as e:  # noqa: BLE001 - un fallo no detiene el resto
        ex.resumen.update({"estado": "error", "mensaje": f"{type(e).__name__}: {e}"})
        if verbose:
            traceback.print_exc()
    return ex


# --------------------------------------------------------------------------- análisis con datos de OpenRocket


def caso_con_aleta_or(caso: Caso, fila_res, pol: pd.DataFrame) -> Caso:
    """Sustituye la aleta analítica por el polígono que OpenRocket aceptó (or_aletas.csv)."""
    if pol is None or pol.empty:
        return caso
    pol = pol.sort_values("orden")
    P = pol[["x_mm", "r_mm"]].to_numpy() * MM
    g = mod_aletas.desde_poligono(caso.params_aleta, P, list(pol["origen"]),
                                  x_LE=float(fila_res["aletas_x_r0_mm"]) * MM,
                                  r_LE=float(fila_res["r_LE_mm"]) * MM,
                                  R_a=float(fila_res["D_popa_mm"]) * MM / 2)
    return caso.con_aletas(g)


def entradas_or(caso: Caso, resumen, perfiles: pd.DataFrame, aletas: pd.DataFrame,
                componentes: pd.DataFrame):
    """(caso con la aleta de OpenRocket, perfil, AeroOR, masas por componente) para analizar_caso."""
    aero = AeroOR(x_CP=float(resumen["x_CP_mm"]) * MM, CNa=float(resumen["CN_alpha_rad"]))
    caso = caso_con_aleta_or(caso, resumen, aletas)
    masas = {r.componente: r.m_g * G for r in componentes.itertuples()
             if r.componente in ("nariz", "cuerpo", "cola", "aletas")}
    return caso, perfiles, aero, masas


def analizar_con_or(caso: Caso, ex: ExportOR, rellenos, tol: Tolerancias) -> ResultadoCaso:
    c, perfil, aero, masas = entradas_or(caso, ex.resumen, pd.DataFrame(ex.perfiles), pd.DataFrame(ex.aletas),
                                         pd.DataFrame(ex.componentes))
    return analizar_caso(c, rellenos=rellenos, perfil=perfil, aero_or=aero, masas_or=masas, tol=tol)


# --------------------------------------------------------------------------- selección y validación


def valores_barrido(cfg: Config) -> pd.DataFrame:
    """config_id → L_mm, D_mm y los valores de cada ruta del barrido (para agrupar)."""
    return pd.DataFrame([{"config_id": c.id, "L_mm": c.geom.L / MM, "D_mm": c.geom.D / MM, **c.barrido}
                         for c in cfg.casos])


@dataclass
class Validacion:
    exportes: dict[str, ExportOR]
    resultados_or: pd.DataFrame  # esquema RESULTADOS, con el CP de OpenRocket
    tabla: pd.DataFrame  # interno frente a OpenRocket por (config, relleno) validado
    ganadores: pd.DataFrame  # mejor por grupo con los valores de OpenRocket
    resultados: dict[str, ResultadoCaso] = field(default_factory=dict)


def _grupos(df: pd.DataFrame, claves: list[str]):
    if not claves:
        yield (), df
        return
    for k, d in df.groupby(claves, sort=True, dropna=False):
        yield (k if isinstance(k, tuple) else (k,)), d


def _puede_superar(fila_int, mejor_or) -> bool:
    m_i, m_o = round(float(fila_int["m_total_g"])), round(float(mejor_or["m_total_g"]))
    return m_i > m_o or (m_i == m_o and float(fila_int["SM_cal"]) > float(mejor_or["SM_cal"]))


def validar(cfg: Config, resultados: pd.DataFrame, exportar, n_verif: int, agrupar: list[str],
            max_iter: int, rellenos=None, log=print) -> Validacion:
    """`exportar(caso) -> ExportOR` corre OpenRocket (inyectado para poder probar sin JVM)."""
    exigir = (cfg.raw.get("optimizacion") or {}).get("exigir") or []
    tol = Tolerancias.desde(cfg.verificacion)
    rellenos = list(rellenos or cfg.rellenos)
    por_nombre = {r.nombre: r for r in rellenos}
    vals = valores_barrido(cfg)
    faltan = [k for k in agrupar if k not in vals.columns]
    if faltan:
        raise KeyError(f"validacion_or.agrupar_por: {faltan} no son L_mm, D_mm ni rutas del barrido")
    res_int = resultados[resultados["relleno"].isin(por_nombre)]
    rank_int = ranking(res_int, exigir).merge(vals, on="config_id", how="left")
    exportes: dict[str, ExportOR] = {}
    analisis: dict[str, ResultadoCaso] = {}
    filas_or: dict[tuple[str, str], dict] = {}

    def asegurar(cids):
        for cid in cids:
            if cid in exportes:
                continue
            caso = cfg.caso(cid)
            ex = exportar(caso)
            exportes[cid] = ex
            if not ex.ok:
                log(f"    {cid}: OpenRocket falló ({ex.resumen.get('estado')}: {ex.resumen.get('mensaje')})")
                continue
            res = analizar_con_or(caso, ex, rellenos, tol)
            analisis[cid] = res
            for rr in res.rellenos:
                filas_or[(cid, rr.relleno.nombre)] = rr.fila
            f0 = res.rellenos[0].fila if res.rellenos else {}
            log(f"    {cid}: x_CP OR {res.x_CP / MM:7.2f} mm (interno {res.cp_int.x_CP / MM:7.2f}) · "
                f"m_total {f0.get('m_total_g', math.nan):7.0f} g · SM {f0.get('SM_cal', math.nan):.2f}")

    tabla, ganadores = [], []
    claves = ["relleno"] + list(agrupar)
    for clave, d_int in _grupos(rank_int, claves):
        etiqueta = ", ".join(f"{k}={v}" for k, v in zip(claves, clave))
        rel = clave[0]
        validados: list[str] = []
        pendientes = list(d_int["config_id"])
        mejor = None
        for it in range(max_iter + 1):
            if it == 0:
                lote = pendientes[:n_verif]
            else:
                candidatos = [c for c in pendientes if c not in validados
                              and _puede_superar(d_int.set_index("config_id").loc[c], mejor)]
                if not candidatos:
                    break
                lote = candidatos[:n_verif]
            log(f"  [{etiqueta}] iteración {it}: valida {len(lote)}")
            asegurar(lote)
            validados += [c for c in lote if c not in validados]
            d_or = pd.DataFrame([filas_or[(c, rel)] for c in validados if (c, rel) in filas_or])
            r_or = ranking(d_or, exigir) if len(d_or) else d_or
            if len(r_or) == 0:
                mejor = None
                if it >= max_iter:
                    break
                mejor = {"m_total_g": -math.inf, "SM_cal": -math.inf}
                continue
            mejor = r_or.iloc[0]
        gan = mejor["config_id"] if mejor is not None and "config_id" in mejor else None
        puestos_int = dict(zip(d_int["config_id"], d_int["puesto"]))
        for cid in validados:
            fi = res_int[(res_int["config_id"] == cid) & (res_int["relleno"] == rel)].iloc[0]
            fo = filas_or.get((cid, rel), {})
            tabla.append({
                "grupo": etiqueta, "config_id": cid, "relleno": rel, "puesto_interno": puestos_int.get(cid),
                "estado_or": exportes[cid].resumen.get("estado"),
                "m_total_int_g": fi.get("m_total_g"), "m_total_or_g": fo.get("m_total_g"),
                "m_relleno_int_g": fi.get("m_relleno_g"), "m_relleno_or_g": fo.get("m_relleno_g"),
                "SM_int_cal": fi.get("SM_cal"), "SM_or_cal": fo.get("SM_cal"),
                "x_CP_int_mm": fi.get("x_CP_mm"), "x_CP_or_mm": fo.get("x_CP_mm"),
                "dx_CP_mm": (fo.get("x_CP_mm", math.nan) - fi.get("x_CP_mm", math.nan)),
                "limitante_int": fi.get("limitante_masa"), "limitante_or": fo.get("limitante_masa"),
                "dif_masa_or_pct": fo.get("dif_masa_or_pct"), "banderas_or": fo.get("banderas"),
                "or_warnings": exportes[cid].resumen.get("or_warnings"), "ganador": cid == gan,
            })
        if gan is not None:
            ganadores.append({"grupo": etiqueta, **dict(zip(claves, clave)), **filas_or[(gan, rel)],
                              "n_validados": len(validados)})
            log(f"  [{etiqueta}] ganador validado: {gan}")
        else:
            log(f"  [{etiqueta}] ningún validado cumple con OpenRocket")
    res_or = pd.DataFrame(list(filas_or.values()))
    return Validacion(exportes, res_or, pd.DataFrame(tabla), pd.DataFrame(ganadores), analisis)
