#!/usr/bin/env python
"""Script 1: barrido con el modelo propio (sin JVM). Volumen utilizable, masa de lastre, CG, SM,
ventana de ℓ, trim y ranking por configuración.

    python scripts/01_volumen_lastre.py --config config/config.yaml [--solo id1,id2]
        [--rellenos plomo_macizo] [--sin-figuras | --sin-dibujos] [--figuras-detalle]
        [--con-orlab | --or-resumen ... --or-perfiles ... --or-aletas ... --or-componentes ...]

Por defecto usa el perfil analítico y el Barrowman interno: es rápido y no necesita OpenRocket.
Los ganadores se validan después con OpenRocket (scripts/02_validar_openrocket.py).
`--con-orlab` lee en cambio los or_*.csv de una exportación completa (script 02 `--todos`).

Figuras: un dibujo por configuración y las curvas del barrido; el perfil y las curvas CG/SM de
cada configuración solo con `--figuras-detalle`.
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
from sensor_lastre.analisis import Tolerancias, analizar_caso  # noqa: E402
from sensor_lastre.config import ConfigError, cargar  # noqa: E402
from sensor_lastre.validacion_or import entradas_or  # noqa: E402

MM = 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config.yaml"))
    p.add_argument("--or-resumen")
    p.add_argument("--or-perfiles")
    p.add_argument("--or-aletas")
    p.add_argument("--or-componentes")
    p.add_argument("--con-orlab", action="store_true",
                   help="lee los or_*.csv de dir_datos (exportación completa del script 02 --todos)")
    p.add_argument("--sin-orlab", action="store_true", help=argparse.SUPPRESS)  # ya es el modo por defecto
    p.add_argument("--rellenos", help="lista separada por comas (por defecto, todos)")
    p.add_argument("--solo", help="ids de configuración separados por comas")
    p.add_argument("--sin-figuras", action="store_true")
    p.add_argument("--sin-dibujos", action="store_true", help="omite el dibujo de cada configuración")
    p.add_argument("--figuras-detalle", action="store_true",
                   help="también el perfil y las curvas CG/SM de cada configuración (lento)")
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
    dibujos = figuras and not a.sin_dibujos
    detalle = figuras and a.figuras_detalle

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
    con_or = a.con_orlab or any(getattr(a, "or_" + k) for k in ("resumen", "perfiles", "aletas", "componentes"))
    if con_or:
        rutas = {k: Path(getattr(a, "or_" + k)) if getattr(a, "or_" + k) else dir_datos / f"or_{k}.csv"
                 for k in ("resumen", "perfiles", "aletas", "componentes")}
        faltan = [str(p) for p in rutas.values() if not p.exists()]
        if faltan:
            print(f"No se encuentran {faltan}. Exporte con scripts/02_validar_openrocket.py --todos "
                  "o corra sin --con-orlab (modelo propio).", file=sys.stderr)
            return 2
        or_res = pd.read_csv(rutas["resumen"]).set_index("config_id")
        or_perf = pd.read_csv(rutas["perfiles"])
        or_ale = pd.read_csv(rutas["aletas"])
        or_comp = pd.read_csv(rutas["componentes"])

    tol = Tolerancias.desde(cfg.verificacion)
    filas, presup, capas, ventanas, valores, confs = [], [], [], [], [], []
    relleno_dibujo = cfg.salida.get("relleno_dibujo", rellenos[0].nombre if rellenos else "")
    if figuras:
        dir_figs.mkdir(parents=True, exist_ok=True)
    for caso in casos:
        perfil = aero = masas_or = None
        if or_res is not None:
            if caso.id in or_res.index and str(or_res.loc[caso.id, "estado"]) == "ok":
                caso, perfil, aero, masas_or = entradas_or(
                    caso, or_res.loc[caso.id], or_perf[or_perf["config_id"] == caso.id],
                    or_ale[or_ale["config_id"] == caso.id], or_comp[or_comp["config_id"] == caso.id])
            else:
                print(f"[{caso.id}] sin resultado válido de OpenRocket: se usa perfil analítico "
                      "y Barrowman interno", file=sys.stderr)
        res = analizar_caso(caso, rellenos=rellenos, perfil=perfil, aero_or=aero, masas_or=masas_or,
                            tol=tol)
        valores.append({"config_id": caso.id, **caso.barrido})
        g_a = caso.geom.aletas
        rr_dib = next((rr for rr in res.rellenos if rr.relleno.nombre == relleno_dibujo),
                      res.rellenos[0] if res.rellenos else None)
        ruta_dib = dir_figs / f"dibujo__{caso.id}.png"
        f0 = res.rellenos[0].fila if res.rellenos else {}
        confs.append({
            "config_id": caso.id, "L_mm": caso.geom.L / MM, "L_total_mm": caso.geom.L_total / MM,
            "D_mm": caso.geom.D / MM, "D_acostado_mm": f0.get("D_acostado_mm"),
            "D_parado_mm": f0.get("D_parado_mm"),
            "h_env_mm": f0.get("h_env_mm"), "w_env_mm": f0.get("w_env_mm"), "aletas_n": g_a.params.n,
            "aleta_h_tip_mm": g_a.h / MM, "aleta_extension_mm": g_a.params.e / MM,
            "aleta_rotacion_deg": math.degrees(g_a.params.rotacion), "r_tip_mm": g_a.r_tip / MM,
            "x_CP_mm": res.x_CP / MM, "x_CP_fuente": res.x_CP_fuente,
            "dibujo": str(ruta_dib.relative_to(base)) if dibujos and ruta_dib.is_relative_to(base) else
            (str(ruta_dib) if dibujos else ""),
        })
        if dibujos:
            from sensor_lastre.figuras import fig_dibujo
            fig_dibujo(res, rr_dib, ruta_dib)
        capas += res.filas_masas_capas()
        for rr in res.rellenos:
            filas.append(rr.fila)
            presup += rr.presupuestos
            _escribir(rr.curva, esquemas.CURVA, dir_datos / "curvas" / f"{caso.id}__{rr.relleno.nombre}.csv")
            if rr.ventana is not None:
                ventanas += [{"config_id": caso.id, "relleno": rr.relleno.nombre,
                              "desde_mm": i0 / MM, "hasta_mm": i1 / MM} for i0, i1 in rr.ventana.intervalos]
            if detalle and rr.ventana is not None:
                from sensor_lastre.figuras import fig_cg_sm
                fig_cg_sm(res, rr, dir_figs / f"cg_sm__{caso.id}__{rr.relleno.nombre}.png")
        if detalle:
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
    _escribir(pd.DataFrame(confs), esquemas.CONFIGURACIONES, dir_datos / "configuraciones.csv")
    from sensor_lastre.optimizacion import ranking
    exigir = (cfg.raw.get("optimizacion") or {}).get("exigir") or []
    opt = ranking(pd.read_csv(dir_datos / "resultados.csv"), exigir)
    _escribir(opt, esquemas.OPTIMO, dir_datos / "optimo.csv")
    print(f"\nÓptimo (máxima masa de lastre con SM ≥ {casos[0].estabilidad.SM_min:g}, exigiendo {exigir}):")
    for rel, d in opt[opt["puesto"] == 1].groupby("relleno", sort=False):
        f = d.iloc[0]
        print(f"  {rel:>18}: {f['config_id']}  m_lastre = {f['m_relleno_g']:.0f} g  "
              f"(trasero {f['m_relleno_trasero_g']:.0f} g)  SM = {f['SM_cal']:.2f}  límite: {f['limitante_masa']}")
    rutas = list(cfg.parametros_barrido)
    val = pd.DataFrame(valores).reindex(columns=["config_id"] + rutas).dropna()  # sin los casos fuera del barrido
    if figuras and len(val):
        from sensor_lastre.figuras import fig_barrido
        fig_barrido(resultados[resultados["config_id"].isin(val["config_id"])], val, rutas,
                    casos[0].estabilidad.SM_min, dir_figs)
    print(f"\n{len(casos)} configuraciones × {len(rellenos)} rellenos → {dir_datos}")
    if not con_or:
        print("CP del modelo propio (Barrowman interno). Valide los ganadores con "
              "scripts/02_validar_openrocket.py.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
