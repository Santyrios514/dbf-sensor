#!/usr/bin/env python
"""Script 6: convierte el ganador de la optimización (scripts 04–05) en un YAML del flujo principal
(scripts 01–02), para correr sobre esa geometría el análisis completo de lastre y su validación.

    python scripts/06_ganador_a_config.py --config config/optimizacion.yaml
        [--ranking data_opt/ranking.csv] [--cand CAND_ID] [--salida config/config_ganador.yaml]
        [--dir-datos data/ganador] [--dir-figuras figs/ganador]

El YAML hereda `base_config` (pared, electrónica, masas puntuales, lastre, remolque, vuelo) con la
masa máxima, el SM mínimo, el material y el espesor de aleta de la optimización. La aleta de
sensor_lastre (`freeform_cola`, sin extensión) termina la raíz en la base y la punta a la altura
del borde de salida de la raíz, así que solo representa exactamente las aletas contenidas con
μ = 1 y x_s + c_t = c_r (σ = 1). Si el candidato no cumple eso, el script se detiene y propone el
mejor factible que sí lo cumple.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import copy  # noqa: E402

import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from sensor_lastre.config import ConfigError  # noqa: E402
from sensor_lastre.config import cargar as cargar_lastre  # noqa: E402
from sensor_opt.config import cargar, raw_cuerpo  # noqa: E402
from sensor_opt.exportar import leer_ranking, specs_de_fila  # noqa: E402

TOL_MM = 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--cand", help="cand_id a convertir (por defecto, el ganador verificado o el primero)")
    p.add_argument("--salida", default=str(RAIZ / "config" / "config_ganador.yaml"))
    p.add_argument("--dir-datos", default="data/ganador")
    p.add_argument("--dir-figuras", default="figs/ganador")
    return p.parse_args(argv)


def representable(f) -> bool:
    return abs(f["mu_cr"] - 1.0) < 1e-9 and abs(f["x_s_mm"] + f["c_t_mm"] - f["c_r_mm"]) < TOL_MM


def elegir(df: pd.DataFrame, cand: str | None):
    if cand:
        return df[df["cand_id"] == cand].iloc[0]
    if "ganador" in df and df["ganador"].astype(str).str.lower().eq("true").any():
        return df[df["ganador"].astype(str).str.lower() == "true"].iloc[0]
    return df[df["factible"]].iloc[0]


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    ruta_rank = Path(a.ranking) if a.ranking else cfg.dir_salida() / "ranking.csv"
    df = leer_ranking(ruta_rank)
    f = elegir(df, a.cand)
    if not representable(f):
        alt = df[df["factible"]]
        alt = alt[alt.apply(representable, axis=1)]
        print(f"{f['cand_id']}: la aleta contenida (μ = {f['mu_cr']:g}, x_s + c_t = "
              f"{f['x_s_mm'] + f['c_t_mm']:.2f} mm ≠ c_r = {f['c_r_mm']:.2f} mm) no se representa en "
              "sensor_lastre.", file=sys.stderr)
        if len(alt):
            print(f"Mejor factible representable: {alt.iloc[0]['cand_id']} (use --cand).", file=sys.stderr)
        return 3
    c, _ = specs_de_fila(f)
    raw = raw_cuerpo(cfg, c)
    base = cfg.base_raw
    raw["geometria_base"]["aletas"].update({
        "c_r": {"modo": "absoluto_mm", "valor": float(f["c_r_mm"])},
        "x_tip_le": {"modo": "absoluto_mm", "valor": float(f["x_s_mm"])},
        "h_tip": {"modo": "absoluto_mm", "valor": float(f["h_mm"])},
    })
    raw["estabilidad"] = {**base["estabilidad"], "SM_min_cal": cfg.restricciones.SM_min, "SM_max_cal": None}
    raw["presupuesto_masa"] = copy.deepcopy(base.get("presupuesto_masa") or {"lastre_g": []})
    raw["validacion_or"] = {**(base.get("validacion_or") or {}), "N_verif": 1, "agrupar_por": []}
    raw["salida"] = {**base.get("salida", {}), "dir_datos": a.dir_datos, "dir_figuras": a.dir_figuras}
    raw["openrocket"] = {**base["openrocket"], "ork_base": str(cfg.ruta_base(base["openrocket"]["ork_base"]))}
    cfg_l = cargar_lastre(raw)  # valida el YAML y la aleta antes de escribirlo
    g = cfg_l.casos[0].geom.aletas
    if abs(g.r_tip - float(f["r_tip_mm"]) * 1e-3) > TOL_MM * 1e-3 * 10:
        print(f"Aviso: r_tip = {g.r_tip * 1e3:.3f} mm frente a {f['r_tip_mm']:.3f} mm del ganador.", file=sys.stderr)
    cabecera = (f"# Generado por scripts/06_ganador_a_config.py desde {ruta_rank}\n"
                f"# Candidato: {f['cand_id']}\n# Geometría del ganador de la optimización; el resto hereda "
                f"{cfg.raw['base_config']}.\n")
    Path(a.salida).write_text(cabecera + yaml.safe_dump(raw, allow_unicode=True, sort_keys=False), encoding="utf-8")
    print(f"{f['cand_id']} → {a.salida} (config {cfg_l.casos[0].id}, m_max {cfg.m_max * 1e3:.0f} g)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
