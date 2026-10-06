#!/usr/bin/env python
"""Script 7: planos acotados (A3) de los ganadores de la optimización, sin JVM.

    python scripts/07_planos.py --config config/optimizacion.yaml [--ranking data_opt/ranking.csv]
        [--n N] [--cand ID ...] [--dir-planos planos] [--dpi 300]

Ganadores = candidatos verificados con OpenRocket (script 5) y factibles ahí, ordenados por J con el
CP de OpenRocket (`J_or`). Para los N primeros (y cada --cand del ranking) dibuja
`{dir_planos}/rankNN_<cand_id>.png` y `.pdf`, escribe `{dir_planos}/rankNN_perfil.csv` (x, r_e, r_i y
el polígono de la aleta, para CAD/CFD) y el comparativo de siluetas `comparativo_top.png`. Si el
script 5 no corrió, ordena por J del sustituto y el cajetín dice "sin validar".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import pandas as pd  # noqa: E402

from sensor_lastre.config import ConfigError  # noqa: E402
from sensor_opt import planos  # noqa: E402
from sensor_opt.config import cargar  # noqa: E402
from sensor_opt.exportar import escribir_csv, leer_ranking  # noqa: E402


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--n", type=int, default=5, help="número de ganadores con plano (5 por defecto)")
    p.add_argument("--cand", action="append", default=[], help="cand_id del ranking (se puede repetir)")
    p.add_argument("--dir-planos", help="por defecto, planos/<nombre de salida.dir>")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args(argv)


def _validacion(dir_d: Path) -> dict[str, dict]:
    ruta = dir_d / "verificacion_or.csv"
    if not ruta.exists():
        return {}
    v = pd.read_csv(ruta)
    v = v[v["rol"] == "verificacion"].dropna(subset=["x_CP_or_mm"])
    return {r["cand_id"]: r.to_dict() for _, r in v.drop_duplicates("cand_id", keep="last").iterrows()}


def ganadores(rk: pd.DataFrame, val: dict) -> pd.DataFrame:
    """Factibles ordenados por J_or (verificados y factibles en OR); sin verificación, por J."""
    if val and "J_or" in rk:
        ok = {c for c, v in val.items() if str(v.get("factible_or", "")).lower() in ("true", "1")}
        g = rk[rk["cand_id"].isin(ok) & rk["J_or"].notna()].sort_values(["J_or", "cand_id"], kind="mergesort")
        if len(g):
            # el ganador del script 5 primero (eligió con J_or a precisión completa; el CSV puede venir redondeado)
            if "ganador" in g:
                es = g["ganador"].astype(str).str.lower().isin(["true", "1"])
                g = pd.concat([g[es], g[~es]])
            return g
    return rk[rk["factible"]].sort_values("J", kind="mergesort")


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    dir_d = cfg.dir_salida()
    ruta_rank = Path(a.ranking) if a.ranking else dir_d / "ranking.csv"
    if not ruta_rank.exists():
        print(f"No existe {ruta_rank}: corra primero scripts/04_optimizar_malla.py", file=sys.stderr)
        return 2
    rk = leer_ranking(ruta_rank)
    val = _validacion(dir_d)
    dir_p = Path(a.dir_planos) if a.dir_planos else RAIZ / "planos" / dir_d.name
    dir_p.mkdir(parents=True, exist_ok=True)
    gan = ganadores(rk, val)
    lista = [(i, f) for i, (_, f) in enumerate(gan.head(a.n).iterrows(), start=1)]
    filas = rk.set_index("cand_id")
    for cid in a.cand:
        if cid not in filas.index:
            print(f"--cand {cid}: no está en {ruta_rank}", file=sys.stderr)
            return 2
        if cid not in [f["cand_id"] for _, f in lista]:
            pos = list(gan["cand_id"]).index(cid) + 1 if cid in set(gan["cand_id"]) else None
            lista.append((pos, rk[rk["cand_id"] == cid].iloc[0]))
    print(f"Planos de {len(lista)} candidatos ({'verificados con OpenRocket' if val else 'sin validar'}) → {dir_p}")
    hechos = []
    for puesto, f in lista:
        prefijo = f"rank{puesto:02d}" if puesto is not None else "cand"
        pl = planos.plano(cfg, f, dir_p / f"{prefijo}_{f['cand_id']}", puesto=puesto, validacion=val.get(f["cand_id"]),
                          dpi=a.dpi)
        escribir_csv(planos.perfil_csv(cfg, f), dir_p / f"{prefijo}_perfil.csv")
        hechos.append((prefijo, f))
        print(f"  {pl.rutas[0].name} (+ .pdf){' · INFACTIBLE' if pl.infactible else ''}")
    r = planos.comparativo(cfg, hechos, dir_p / "comparativo_top.png")
    if r is not None:
        print(f"  {r.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
