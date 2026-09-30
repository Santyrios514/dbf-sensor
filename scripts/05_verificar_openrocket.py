#!/usr/bin/env python
"""Script 5 (v3): verificación y calibración con OpenRocket del ranking del script 4.

    python scripts/05_verificar_openrocket.py --config config/optimizacion.yaml [--ranking data_opt/ranking.csv]
        [--procesos N] [--sin-figuras]

Verifica los mejores y una muestra estratificada, calibra el CP sustituto, recalcula la malla si
hace falta, elige el ganador con el CP de OpenRocket y reescribe ranking.csv. Escribe además
verificacion_or.csv, calibracion.json y los `.ork` de los mejores (con lastre y electrónica como
Mass components), y comprueba que el `.ork` del ganador reproduce su CP al reabrirlo.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_lastre.config import ConfigError  # noqa: E402
from sensor_opt import exportar  # noqa: E402
from sensor_opt.config import cargar  # noqa: E402
from sensor_opt.esquemas import VERIFICACION  # noqa: E402
from sensor_opt.verificacion import PuenteOpt, verificar_candidato, verificar_y_calibrar  # noqa: E402

MM = 1e-3
TOL_REAPERTURA_MM = 0.5  # v3 T10


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--procesos", type=int)
    p.add_argument("--sin-figuras", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    dir_d, dir_f = cfg.dir_salida(), cfg.dir_figuras()
    ruta_rank = Path(a.ranking) if a.ranking else dir_d / "ranking.csv"
    if not ruta_rank.exists():
        print(f"No existe {ruta_rank}: corra primero scripts/04_optimizar_malla.py", file=sys.stderr)
        return 2
    ranking = exportar.leer_ranking(ruta_rank)
    orc = cfg.base_raw.get("openrocket") or {}
    ork = cfg.ruta_base(orc.get("ork_base", "modelos/analisis_vol_int.ork"))

    import orlab
    with orlab.OpenRocketInstance(jar_path=orc.get("jar"), log_level="ERROR") as inst:
        pr = PuenteOpt(inst, orc["nombres_componentes"], ork)
        print(f"OpenRocket {pr.version_or} · {len(ranking):,} candidatos en {ruta_rank}")
        res = verificar_y_calibrar(cfg, pr, ranking, a.procesos)
        df, ver, cal = res.ranking, res.verificacion, res.calibracion
        df["ganador"] = df["cand_id"] == res.ganador

        # .ork de los mejores verificados (por J con el CP de OpenRocket)
        n_ork = int(cfg.salida.get("exportar_ork_top", 5))
        mejores = ver[ver["factible_or"].fillna(False).astype(bool)].sort_values(["J_or", "cand_id"]).head(n_ork)
        filas = df.set_index("cand_id")
        reapertura = math.nan
        for i, cid in enumerate(mejores["cand_id"], start=1):
            ruta = dir_d / "ork" / f"rank{i:02d}_{cid}.ork"
            fila = filas.loc[cid].copy()
            fila["cand_id"] = cid
            r = verificar_candidato(cfg, pr, fila, cal, guardar=ruta)
            if cid == res.ganador:
                pr.cargar(ruta)
                x_re = pr.aero(_mach(cfg)).x_CP
                reapertura = abs(x_re / MM - r["x_CP_or_mm"])
            print(f"  .ork → {ruta.name}")

    exportar.escribir_csv(df, dir_d / "ranking.csv", exportar.RANKING + ["ganador"])
    exportar.escribir_csv(df[df["en_pareto"].astype(bool)], dir_d / "pareto.csv", exportar.RANKING)
    exportar.escribir_csv(ver, dir_d / "verificacion_or.csv", VERIFICACION)
    exportar.escribir_json({**cal.a_dict(), "historial": res.historial, "ganador": res.ganador,
                            "reapertura_ganador_dif_CP_mm": reapertura}, dir_d / "calibracion.json")

    if res.ganador is None:
        print("\nNingún candidato verificado cumple SM en [SM_min, SM_max] con el CP de OpenRocket.")
        return 1
    g = df[df["ganador"]].iloc[0]
    print(f"\nGanador (verificado con OpenRocket): {g['cand_id']}\n"
          f"  m_total {g['m_total_or_g']:.0f} g · D_ap {g['D_ap_mm']:.1f} mm · k {g['k']:.3f} · "
          f"SM_OR {g['SM_or_cal']:.2f} (sustituto {g['SM_cal']:.2f}) · x_CP OR {g['x_CP_or_mm']:.1f} mm · "
          f"θ_eq {g['theta_eq_deg']:.1f}° · dif. masa {g['dif_masa_or_pct']:+.2f} %")
    print(f"  .ork del ganador reabierto: |ΔCP| = {reapertura:.3f} mm (tolerancia {TOL_REAPERTURA_MM} mm)")
    if not a.sin_figuras:
        dir_f.mkdir(parents=True, exist_ok=True)
        exportar.fig_calibracion(ver, cal, dir_f / "calibracion.png")
        exportar.fig_pareto(df.sort_values(["ganador"], ascending=False), dir_f / "pareto.png")
        exportar.fig_factibilidad(df, dir_f / "factibilidad.png")
        exportar.fig_ganador(cfg, g, dir_f, cal, x_CP=g["x_CP_or_mm"] * MM)
        exportar.fig_sensibilidad(cfg, g, dir_f / "sensibilidad.png", cal)
        print(f"Figuras → {dir_f}")
    print(f"Datos → {dir_d}")
    return 0 if reapertura <= TOL_REAPERTURA_MM else 4


def _mach(cfg) -> float:
    from sensor_lastre.config import cargar as cargar_lastre
    from sensor_opt.config import raw_cuerpo
    return cargar_lastre(raw_cuerpo(cfg, cfg.cuerpos()[0])).casos[0].vuelo.mach


if __name__ == "__main__":
    sys.exit(main())
