#!/usr/bin/env python
"""Script 4 (v3): barrido en malla de la geometría del sensor con el CP sustituto (sin JVM).

    python scripts/04_optimizar_malla.py --config config/optimizacion.yaml [--forzar] [--sin-refinamiento]
        [--procesos N] [--sin-figuras]

Salidas en `salida.dir` (ranking.csv, pareto.csv, tolerancias_implicitas.json) y figuras en
`salida.dir_figuras`. El script 5 verifica con OpenRocket y reescribe el ranking final.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_lastre.config import ConfigError  # noqa: E402
from sensor_opt import barrido, exportar  # noqa: E402
from sensor_opt.config import cargar  # noqa: E402
from sensor_opt.objetivo import tolerancias  # noqa: E402


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--forzar", action="store_true", help="corre aunque se supere max_evaluaciones")
    p.add_argument("--sin-refinamiento", action="store_true")
    p.add_argument("--procesos", type=int, help="sobrescribe ejecucion.procesos")
    p.add_argument("--sin-figuras", action="store_true")
    return p.parse_args(argv)


def _sugerencia(cfg) -> str:
    m = cfg.malla
    largos = sorted(((len(v), k) for k, v in m.items() if isinstance(v, list)), reverse=True)
    return ", ".join(f"{k} ({n} valores)" for n, k in largos[:3])


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    n_c, n_a = len(cfg.cuerpos()), len(cfg.aletas())
    n = n_c * n_a
    procesos = a.procesos or cfg.procesos
    lim = int(cfg.ejecucion.get("max_evaluaciones", 2_000_000))
    print(f"Malla: {n_c} cuerpos × {n_a} aletas = {n:,} evaluaciones (límite {lim:,}); {procesos} procesos.")
    print(f"Estimado: ~{n * 4e-3 / procesos / 60:.0f} min a ~4 ms por candidato (sin refinamiento).")
    if n > lim and not a.forzar:
        print(f"La malla supera max_evaluaciones. Reduzca las listas más largas ({_sugerencia(cfg)}) "
              "o use --forzar.", file=sys.stderr)
        return 3

    t0 = time.time()

    def progreso(i, total):
        if i == total or i % max(1, total // 20) == 0:
            print(f"  {i}/{total} cuerpos · {time.time() - t0:.0f} s", flush=True)

    df = barrido.optimizar(cfg, refinar=not a.sin_refinamiento, procesos=procesos, progreso=progreso)
    dir_d, dir_f = cfg.dir_salida(), cfg.dir_figuras()
    exportar.exportar_ranking(cfg, df, dir_d)
    fac = df[df["factible"]]
    print(f"\n{len(df):,} candidatos ({int(df['refinamiento'].sum()):,} del refinamiento), "
          f"{len(fac):,} factibles, {int(df['en_pareto'].sum()):,} en el frente de Pareto · {time.time() - t0:.0f} s")
    tol = tolerancias(cfg)
    print("\nTolerancias implícitas de los pesos:")
    for t in tol["lectura"]:
        print("  - " + t)
    if not fac.empty:
        g = fac.iloc[0]
        print(f"\nMejor (CP sustituto, sin verificar): {g['cand_id']}\n"
              f"  m_total {g['m_total_g']:.0f} g · D_ap {g['D_ap_mm']:.1f} mm · k {g['k']:.3f} · "
              f"SM {g['SM_cal']:.2f} · θ_eq {g['theta_eq_deg']:.1f}° · J {g['J']:.2f} · "
              f"restricción activa: {g['restriccion_activa']}")
    else:
        print("\nNingún candidato factible.")
    if not a.sin_figuras:
        rutas = exportar.figuras_barrido(cfg, df, dir_f)
        print(f"\n{len(rutas)} figuras → {dir_f}")
    print(f"Datos → {dir_d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
