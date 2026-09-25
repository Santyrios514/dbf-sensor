#!/usr/bin/env python
"""Script 2: volumen utilizable, masa de lastre, CG, SM, ventana de ℓ y trim por configuración.

    python scripts/02_volumen_lastre.py --config config/config.yaml \
        [--or-resumen data/or_resumen.csv --or-perfiles data/or_perfiles.csv \
         --or-aletas data/or_aletas.csv --or-componentes data/or_componentes.csv | --sin-orlab] \
        [--rellenos plomo_macizo,W90_macizo] [--sin-figuras]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import pandas as pd  # noqa: E402

from sensor_lastre import aletas as mod_aletas  # noqa: E402
from sensor_lastre import esquemas  # noqa: E402
from sensor_lastre.analisis import AeroOR, Tolerancias, analizar_caso  # noqa: E402
from sensor_lastre.config import ConfigError, cargar  # noqa: E402

MM = 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config.yaml"))
    p.add_argument("--or-resumen")
    p.add_argument("--or-perfiles")
    p.add_argument("--or-aletas")
    p.add_argument("--or-componentes")
    p.add_argument("--sin-orlab", action="store_true", help="perfiles analíticos + Barrowman interno")
    p.add_argument("--rellenos", help="lista separada por comas (por defecto, todos)")
    p.add_argument("--solo", help="ids de configuración separados por comas")
    p.add_argument("--sin-figuras", action="store_true")
    p.add_argument("--dir-datos", help="sobrescribe salida.dir_datos")
    p.add_argument("--dir-figuras", help="sobrescribe salida.dir_figuras")
    return p.parse_args(argv)


def _ruta(base: Path, valor: str) -> Path:
    p = Path(valor)
    return p if p.is_absolute() else base / p


def _escribir(df: pd.DataFrame, columnas: list[str], ruta: Path):
    df = df.reindex(columns=columnas)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(ruta, index=False, float_format="%.6g")


def _aletas_or(caso, fila_res, pol: pd.DataFrame):
    """Sustituye la aleta analítica por el polígono que OpenRocket aceptó (or_aletas.csv)."""
    if pol.empty:
        return caso
    pol = pol.sort_values("orden")
    P = pol[["x_mm", "r_mm"]].to_numpy() * MM
    g = mod_aletas.desde_poligono(caso.params_aleta, P, list(pol["origen"]),
                                  x_LE=float(fila_res["aletas_x_r0_mm"]) * MM,
                                  r_LE=float(fila_res["r_LE_mm"]) * MM,
                                  R_a=float(fila_res["D_popa_mm"]) * MM / 2)
    return caso.con_aletas(g)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    base = Path(a.config).resolve().parents[1]
    dir_datos = Path(a.dir_datos) if a.dir_datos else _ruta(base, cfg.salida.get("dir_datos", "data"))
    dir_figs = Path(a.dir_figuras) if a.dir_figuras else _ruta(base, cfg.salida.get("dir_figuras", "figs"))
    figuras = cfg.salida.get("figuras", True) and not a.sin_figuras

    rellenos = list(cfg.rellenos)
    if a.rellenos:
        pedidos = [r.strip() for r in a.rellenos.split(",")]
        desconocidos = set(pedidos) - {r.nombre for r in rellenos}
        if desconocidos:
            print(f"rellenos desconocidos: {sorted(desconocidos)}", file=sys.stderr)
            return 2
        rellenos = [r for r in rellenos if r.nombre in pedidos]
    casos = cfg.casos
    if a.solo:
        ids = [s.strip() for s in a.solo.split(",")]
        casos = [c for c in casos if c.id in ids]

    or_res = or_perf = or_ale = or_comp = None
    if not a.sin_orlab:
        rutas = {k: Path(getattr(a, "or_" + k)) if getattr(a, "or_" + k) else dir_datos / f"or_{k}.csv"
                 for k in ("resumen", "perfiles", "aletas", "componentes")}
        faltan = [str(p) for p in rutas.values() if not p.exists()]
        if faltan:
            print(f"No se encuentran {faltan}. Corre el script 1 o usa --sin-orlab.", file=sys.stderr)
            return 2
        or_res = pd.read_csv(rutas["resumen"]).set_index("config_id")
        or_perf = pd.read_csv(rutas["perfiles"])
        or_ale = pd.read_csv(rutas["aletas"])
        or_comp = pd.read_csv(rutas["componentes"])

    tol = Tolerancias.desde(cfg.verificacion)
    filas, presup, capas, ventanas, valores = [], [], [], [], []
    if figuras:
        dir_figs.mkdir(parents=True, exist_ok=True)
    for caso in casos:
        perfil = aero = masas_or = None
        if or_res is not None:
            if caso.id in or_res.index and str(or_res.loc[caso.id, "estado"]) == "ok":
                f = or_res.loc[caso.id]
                aero = AeroOR(x_CP=float(f["x_CP_mm"]) * MM, CNa=float(f["CN_alpha_rad"]))
                perfil = or_perf[or_perf["config_id"] == caso.id]
                caso = _aletas_or(caso, f, or_ale[or_ale["config_id"] == caso.id])
                comp = or_comp[or_comp["config_id"] == caso.id]
                masas_or = {r.componente: r.m_g * 1e-3 for r in comp.itertuples()
                            if r.componente in ("nariz", "cuerpo", "cola", "aletas")}
            else:
                print(f"[{caso.id}] sin resultado válido de OpenRocket: se usa perfil analítico "
                      "y Barrowman interno", file=sys.stderr)
        res = analizar_caso(caso, rellenos=rellenos, perfil=perfil, aero_or=aero, masas_or=masas_or,
                            tol=tol)
        valores.append({"config_id": caso.id, **caso.barrido})
        capas += res.filas_masas_capas()
        for rr in res.rellenos:
            filas.append(rr.fila)
            presup += rr.presupuestos
            _escribir(rr.curva, esquemas.CURVA, dir_datos / "curvas" / f"{caso.id}__{rr.relleno.nombre}.csv")
            if rr.ventana is not None:
                ventanas += [{"config_id": caso.id, "relleno": rr.relleno.nombre,
                              "desde_mm": i0 / MM, "hasta_mm": i1 / MM} for i0, i1 in rr.ventana.intervalos]
            if figuras and rr.ventana is not None:
                from sensor_lastre.figuras import fig_cg_sm
                fig_cg_sm(res, rr, dir_figs / f"cg_sm__{caso.id}__{rr.relleno.nombre}.png")
        if figuras:
            from sensor_lastre.figuras import fig_perfil
            fig_perfil(res, res.rellenos[0] if res.rellenos else None, dir_figs / f"perfil__{caso.id}.png")
        f0 = res.rellenos[0].fila if res.rellenos else {}
        print(f"{caso.id:>14}  V_ext={f0.get('V_ext_cm3', math.nan):7.1f} cm³  "
              f"ℓ_geo={res.lim.ell_geo / MM:6.1f} mm ({res.lim.limitante})  "
              f"x_CP={res.x_CP / MM:6.1f} mm [{res.x_CP_fuente}]  "
              f"banderas={sorted({b for rr in res.rellenos for b in rr.fila['banderas'].split(';') if b})}")

    resultados = pd.DataFrame(filas)
    _escribir(resultados, esquemas.RESULTADOS, dir_datos / "resultados.csv")
    _escribir(pd.DataFrame(capas), esquemas.MASAS_CAPAS, dir_datos / "masas_capas.csv")
    _escribir(pd.DataFrame(presup), esquemas.PRESUPUESTOS, dir_datos / "presupuestos.csv")
    _escribir(pd.DataFrame(ventanas), esquemas.VENTANAS, dir_datos / "ventanas.csv")
    if figuras and len(resultados):
        from sensor_lastre.figuras import fig_barrido
        fig_barrido(resultados, pd.DataFrame(valores), list(cfg.parametros_barrido),
                    casos[0].estabilidad.SM_min, dir_figs)
    print(f"\n{len(casos)} configuraciones × {len(rellenos)} rellenos → {dir_datos}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
