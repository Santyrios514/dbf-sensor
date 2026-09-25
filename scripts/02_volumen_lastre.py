#!/usr/bin/env python
"""Script 2: volumen utilizable, masa de lastre, CG, SM, ventana de ℓ y trim por configuración.

    python scripts/02_volumen_lastre.py --config config/config.yaml \
        [--or-resumen data/or_resumen.csv --or-perfiles data/or_perfiles.csv | --sin-orlab] \
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

from sensor_lastre import esquemas  # noqa: E402
from sensor_lastre.analisis import AeroOR, analizar_caso  # noqa: E402
from sensor_lastre.config import ConfigError, cargar  # noqa: E402

MM = 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config.yaml"))
    p.add_argument("--or-resumen")
    p.add_argument("--or-perfiles")
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

    or_res = or_perf = None
    if not a.sin_orlab:
        r_path = Path(a.or_resumen) if a.or_resumen else dir_datos / "or_resumen.csv"
        p_path = Path(a.or_perfiles) if a.or_perfiles else dir_datos / "or_perfiles.csv"
        if not (r_path.exists() and p_path.exists()):
            print(f"No se encuentran {r_path} y {p_path}. Corre el script 1 o usa --sin-orlab.",
                  file=sys.stderr)
            return 2
        or_res = pd.read_csv(r_path).set_index("config_id")
        or_perf = pd.read_csv(p_path)

    tol_cp = cfg.verificacion.get("dif_CP_max_frac_L")
    filas, presup, capas, ventanas = [], [], [], []
    if figuras:
        dir_figs.mkdir(parents=True, exist_ok=True)
    for caso in casos:
        perfil = aero = None
        if or_res is not None:
            if caso.id in or_res.index and str(or_res.loc[caso.id, "estado"]) == "ok":
                f = or_res.loc[caso.id]
                aero = AeroOR(x_CP=float(f["x_CP_mm"]) * MM, CNa=float(f["CN_alpha_rad"]))
                perfil = or_perf[or_perf["config_id"] == caso.id]
            else:
                print(f"[{caso.id}] sin resultado válido de OpenRocket: se usa perfil analítico "
                      "y Barrowman interno", file=sys.stderr)
        res = analizar_caso(caso, rellenos=rellenos, perfil=perfil, aero_or=aero,
                            tol_dif_CP_frac_L=tol_cp)
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
        from sensor_lastre.figuras import fig_mapas
        fig_mapas(resultados, dir_figs)
    print(f"\n{len(casos)} configuraciones × {len(rellenos)} rellenos → {dir_datos}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
