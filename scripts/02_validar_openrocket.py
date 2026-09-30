#!/usr/bin/env python
"""Script 2: valida con OpenRocket los ganadores del barrido del script 1.

    python scripts/02_validar_openrocket.py --config config/config.yaml [--n-verif N]
        [--guardar-ork] [--sin-figuras] [--dir-datos ...] [--dir-figuras ...]
    python scripts/02_validar_openrocket.py --regresion      # la base contra los valores de referencia
    python scripts/02_validar_openrocket.py --todos [--solo id1,id2]   # exporta todas (or_*.csv)

Por grupo (relleno × `validacion_or.agrupar_por`) aplica al `.ork` base los N mejores del ranking
del modelo propio, los vuelve a analizar con el CP, el perfil, la aleta y las masas de
OpenRocket, y sigue validando mientras algún candidato sin validar pueda superar al mejor
validado (hasta `max_iter`). Escribe en dir_datos:

    validacion_or.csv   interno frente a OpenRocket por configuración validada
    ganadores.csv       mejor configuración de cada grupo con los valores de OpenRocket
    resultados_or.csv   filas de resultados.csv recalculadas con OpenRocket (validadas)
    or_*.csv            exportación de OpenRocket de las configuraciones validadas

y, en dir_figuras, el dibujo y el perfil de cada ganador con el CP de OpenRocket.
"""

from __future__ import annotations

import argparse
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import pandas as pd  # noqa: E402

from sensor_lastre import esquemas  # noqa: E402
from sensor_lastre.config import ConfigError, cargar  # noqa: E402
from sensor_lastre.validacion_or import exportar_caso, validar  # noqa: E402
from sensor_lastre.verificacion import comparar  # noqa: E402

MM, G = 1e-3, 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config.yaml"))
    p.add_argument("--ork", help="sobrescribe openrocket.ork_base")
    p.add_argument("--n-verif", type=int, help="sobrescribe validacion_or.N_verif")
    p.add_argument("--todos", action="store_true",
                   help="exporta todas las configuraciones a or_*.csv (sin validar ganadores)")
    p.add_argument("--solo", help="con --todos: ids de configuración separados por comas")
    p.add_argument("--guardar-ork", action="store_true", help="guarda cada variante en dir_datos/ork/{id}.ork")
    p.add_argument("--regresion", action="store_true",
                   help="sin modificar la geometría: M de regresión y comparación con v2 §1.4")
    p.add_argument("--dir-datos", help="sobrescribe salida.dir_datos")
    p.add_argument("--dir-figuras", help="sobrescribe salida.dir_figuras")
    p.add_argument("--sin-figuras", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def _ruta(base: Path, v: str) -> Path:
    p = Path(v)
    return p if p.is_absolute() else base / p


def _escribir(filas, columnas, ruta: Path):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    df = filas if isinstance(filas, pd.DataFrame) else pd.DataFrame(filas)
    df.reindex(columns=columnas).to_csv(ruta, index=False, float_format="%.8g")


def _escribir_or(exportes, dir_datos: Path):
    _escribir([e.resumen for e in exportes], esquemas.OR_RESUMEN, dir_datos / "or_resumen.csv")
    _escribir([f for e in exportes for f in e.perfiles], esquemas.OR_PERFILES, dir_datos / "or_perfiles.csv")
    _escribir([f for e in exportes for f in e.componentes], esquemas.OR_COMPONENTES,
              dir_datos / "or_componentes.csv")
    _escribir([f for e in exportes for f in e.aletas], esquemas.OR_ALETAS, dir_datos / "or_aletas.csv")


def regresion(pr, cfg) -> int:
    """v2 §1.4 / §8.5: masa, CG, CP y estabilidad del .ork sin tocar, a M = mach_regresion_or."""
    vuelo = cfg.casos[0].vuelo
    ref = cfg.verificacion["regresion_or"]
    aero = pr.aero(vuelo.mach_regresion_or, 0.0)
    m = pr.masas()
    D = float(pr.configuracion.getReferenceLength())
    obtenido = {"masa_g": m.m_total / G, "x_CG_mm": m.x_CG_total / MM, "x_CP_mm": aero.x_CP / MM,
                "estabilidad_cal": (aero.x_CP - m.x_CG_total) / D}
    ok = True
    print(f"Regresión OpenRocket {pr.version_or} (M = {vuelo.mach_regresion_or}, AoA = 0):")
    for k, v in obtenido.items():
        r = ref[k]
        pasa = abs(v - r["valor"]) <= r["tol"]
        ok &= pasa
        print(f"  {k:16s} {v:10.4f}   ref {r['valor']} ± {r['tol']}   {'OK' if pasa else 'FALLA'}")
    print(f"  advertencias: {aero.warnings}")
    return 0 if ok else 1


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    base = Path(a.config).resolve().parents[1]
    orc = cfg.openrocket
    ork = Path(a.ork) if a.ork else _ruta(base, orc["ork_base"])
    dir_datos = Path(a.dir_datos) if a.dir_datos else _ruta(base, cfg.salida.get("dir_datos", "data"))
    dir_figs = Path(a.dir_figuras) if a.dir_figuras else _ruta(base, cfg.salida.get("dir_figuras", "figs"))
    vo = cfg.raw.get("validacion_or") or {}
    n_verif = a.n_verif or int(vo.get("N_verif", 3))
    agrupar = list(vo.get("agrupar_por") or [])
    max_iter = int(vo.get("max_iter", 3))

    resultados = None
    if not (a.regresion or a.todos):
        ruta_res = dir_datos / "resultados.csv"
        if not ruta_res.exists():
            print(f"No existe {ruta_res}: corra primero scripts/01_volumen_lastre.py", file=sys.stderr)
            return 2
        resultados = pd.read_csv(ruta_res)

    import orlab
    from sensor_lastre.or_bridge import PuenteOR, sha256

    tol_geo = cfg.verificacion.get("tol_geometria_ork_mm", 0.001) * MM
    with orlab.OpenRocketInstance(jar_path=orc.get("jar")) as inst:
        pr = PuenteOR(inst, orc["nombres_componentes"]).cargar(ork)
        if a.regresion:
            return regresion(pr, cfg)

        # la base leída por OpenRocket debe coincidir con la configuración (v2 §5.1)
        caso_base = next((c for c in cfg.casos if c.id == "base_ork"), None)
        if caso_base is not None:
            dif = comparar(pr.valores(), caso_base, tol_geo)
            if dif:
                print("La geometría del .ork no coincide con la configuración base:\n  - "
                      + "\n  - ".join(dif), file=sys.stderr)
                return 3
            print(f"Base verificada contra {ork.name} (OpenRocket {pr.version_or}).")
        else:
            print("Aviso: no hay configuración 'base_ork'; se omite la verificación de la base.", file=sys.stderr)

        meta = {"or_version": pr.version_or, "orlab_version": pr.version_orlab(), "ork_sha256": sha256(ork),
                "fecha_iso": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        dif_cfg = orc.get("cna_diferencias") or {}
        delta = math.radians(dif_cfg.get("delta_aoa_deg", 1.0)) if dif_cfg.get("activar") else None

        def exportar(caso):
            destino = dir_datos / "ork" / f"{caso.id}.ork" if a.guardar_ork else None
            return exportar_caso(pr, caso, ork, orc, meta, delta, guardar=destino, verbose=a.verbose)

        if a.todos:
            casos = cfg.casos
            if a.solo:
                ids = [s.strip() for s in a.solo.split(",")]
                casos = [c for c in casos if c.id in ids]
            exportes = []
            for caso in casos:
                ex = exportar(caso)
                exportes.append(ex)
                r = ex.resumen
                if ex.ok:
                    print(f"{caso.id:>30}  x_CP = {r['x_CP_mm']:7.2f} mm  C_Nα = {r['CN_alpha_rad']:6.3f}  "
                          f"m = {r['m_vacio_or_g']:6.1f} g  x_CG = {r['x_CG_vacio_or_mm']:6.1f} mm")
                else:
                    print(f"{caso.id:>30}  {r['estado']}: {r['mensaje']}", file=sys.stderr)
            _escribir_or(exportes, dir_datos)
            n_ok = sum(e.ok for e in exportes)
            print(f"\n{n_ok}/{len(exportes)} configuraciones exportadas → {dir_datos}")
            return 0 if n_ok == len(exportes) else 1

        print(f"Validando ganadores: {n_verif} por grupo, agrupando por relleno"
              + (f" y {agrupar}" if agrupar else "") + f", hasta {max_iter} iteraciones.")
        v = validar(cfg, resultados, exportar, n_verif, agrupar, max_iter)

    _escribir_or(list(v.exportes.values()), dir_datos)
    _escribir(v.resultados_or, esquemas.RESULTADOS, dir_datos / "resultados_or.csv")
    _escribir(v.tabla, esquemas.VALIDACION_OR, dir_datos / "validacion_or.csv")
    _escribir(v.ganadores, esquemas.GANADORES + [c for c in esquemas.RESULTADOS if c not in esquemas.GANADORES],
              dir_datos / "ganadores.csv")

    print(f"\n{len(v.exportes)} configuraciones validadas con OpenRocket.")
    t = v.tabla
    if len(t):
        dx = t["dx_CP_mm"].abs()
        print(f"|x_CP OR − x_CP interno|: mediana {dx.median():.2f} mm, máx. {dx.max():.2f} mm")
    for g in v.ganadores.itertuples():
        fila_t = t[(t["config_id"] == g.config_id) & (t["relleno"] == g.relleno)].iloc[0]
        print(f"  [{g.grupo}] {g.config_id}: m_total {g.m_total_g:.0f} g · lastre {g.m_relleno_g:.0f} g · "
              f"SM_OR {g.SM_cal:.2f} (interno {fila_t['SM_int_cal']:.2f}) · límite {g.limitante_masa} · "
              f"puesto interno {int(fila_t['puesto_interno'])}")
    if not a.sin_figuras and cfg.salida.get("figuras", True) and len(v.ganadores):
        from sensor_lastre.figuras import fig_dibujo, fig_perfil
        dir_figs.mkdir(parents=True, exist_ok=True)
        for g in v.ganadores.itertuples():
            res = v.resultados[g.config_id]
            rr = next(r for r in res.rellenos if r.relleno.nombre == g.relleno)
            fig_dibujo(res, rr, dir_figs / f"ganador_or__{g.config_id}__{g.relleno}.png")
            fig_perfil(res, rr, dir_figs / f"ganador_or_perfil__{g.config_id}__{g.relleno}.png")
        print(f"Figuras de los ganadores → {dir_figs}")
    print(f"Datos → {dir_datos}")
    return 0 if len(v.ganadores) else 1


if __name__ == "__main__":
    sys.exit(main())
