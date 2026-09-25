#!/usr/bin/env python
"""Script 3: comparación gruesa de pares (L_cuerpo, D) a partir de los resultados del script 2.

Sirve para elegir con qué L y D correr luego el barrido fino de aletas. Para cada relleno y cada
par (L, D) resume la mejor configuración factible según el mismo criterio de `optimizacion.ranking`
(máxima masa de lastre con SM ≥ SM_min que cabe en la bahía; desempate por SM, luego D acostado)
y cuántas de las variantes de aletas son factibles.

    python scripts/03_comparar_LD.py --config config/config_barrido_LD.yaml \
        [--dir-datos data/barrido_LD] [--dir-figuras figs/barrido_LD] [--relleno plomo_macizo] \
        [--sin-figuras]

Salidas (en dir_datos):
    comparacion_LD.csv           una fila por relleno × (L, D)
    comparacion_LD_pivote.csv    masa máxima [g] de cada relleno en una tabla L × D
Figura (en dir_figuras): comparacion_LD__{relleno}.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from sensor_lastre.config import ConfigError, cargar  # noqa: E402
from sensor_lastre.optimizacion import ranking  # noqa: E402

COLUMNAS = ["relleno", "L_mm", "D_mm", "L_total_mm", "n_variantes", "n_factibles", "frac_factibles",
            "m_relleno_max_g", "m_total_g", "costo_relleno_usd", "SM_cal", "tol_amarre_mm",
            "D_acostado_mm", "D_parado_mm", "cabe_largo", "limitante_masa", "config_mejor",
            "V_int_cm3", "m_casco_g", "SM_inf_max_cal", "puesto_LD"]


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config_barrido_LD.yaml"))
    p.add_argument("--dir-datos", help="sobrescribe salida.dir_datos")
    p.add_argument("--dir-figuras", help="sobrescribe salida.dir_figuras")
    p.add_argument("--relleno", help="solo figura de este relleno (por defecto, todos)")
    p.add_argument("--sin-figuras", action="store_true")
    return p.parse_args(argv)


def _ruta(base: Path, valor: str) -> Path:
    p = Path(valor)
    return p if p.is_absolute() else base / p


def comparar(res: pd.DataFrame, exigir) -> pd.DataFrame:
    """Una fila por relleno × (L, D) con la mejor configuración factible de ese par."""
    res = res[res["config_id"] != "base_ork"]           # la base ya está dentro del barrido
    top = ranking(res, exigir).drop_duplicates(["relleno", "config_id"])
    top = top.merge(res[["relleno", "config_id", "L_mm", "D_mm", "V_int_cm3", "m_casco_g"]],
                    on=["relleno", "config_id"], how="left")
    mejor = top.groupby(["relleno", "L_mm", "D_mm"], sort=False).head(1).set_index(["relleno", "L_mm", "D_mm"])

    g = res.groupby(["relleno", "L_mm", "D_mm"])
    out = pd.DataFrame({
        "L_total_mm": g["L_total_mm"].max(),
        "n_variantes": g["config_id"].nunique(),
        "V_int_cm3": g["V_int_cm3"].median(),
        "m_casco_g": g["m_casco_g"].median(),
        "SM_inf_max_cal": g["SM_inf_cal"].max(),
    })
    out["n_factibles"] = top.groupby(["relleno", "L_mm", "D_mm"])["config_id"].nunique()
    out["n_factibles"] = out["n_factibles"].fillna(0).astype(int)
    out["frac_factibles"] = out["n_factibles"] / out["n_variantes"]
    for c_out, c_in in [("m_relleno_max_g", "m_relleno_g"), ("m_total_g", "m_total_g"),
                        ("costo_relleno_usd", "costo_relleno_usd"), ("SM_cal", "SM_cal"),
                        ("tol_amarre_mm", "tol_amarre_mm"), ("D_acostado_mm", "D_acostado_mm"),
                        ("D_parado_mm", "D_parado_mm"), ("cabe_largo", "cabe_largo"),
                        ("limitante_masa", "limitante_masa"), ("config_mejor", "config_id")]:
        out[c_out] = mejor[c_in]
    out = out.reset_index()
    out = out.sort_values(["relleno", "m_relleno_max_g", "SM_cal", "D_acostado_mm"],
                          ascending=[True, False, False, True], na_position="last", kind="stable")
    out["puesto_LD"] = out.groupby("relleno").cumcount() + 1
    return out.reindex(columns=COLUMNAS)


def _figura(res: pd.DataFrame, comp: pd.DataFrame, relleno: str, exigir, alto_max, ruta: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = comp[comp["relleno"] == relleno]
    Ls, Ds = sorted(c["L_mm"].unique()), sorted(c["D_mm"].unique())
    M = c.pivot(index="L_mm", columns="D_mm", values="m_relleno_max_g").reindex(index=Ls, columns=Ds) / 1000

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1.3]})
    im = a1.imshow(M.to_numpy(), cmap="Blues", origin="lower", aspect="auto")
    a1.set_xticks(range(len(Ds)), [f"{d:g}" for d in Ds])
    a1.set_yticks(range(len(Ls)), [f"{l:g}" for l in Ls])
    a1.set_xlabel("D cuerpo [mm]")
    a1.set_ylabel("L cuerpo [mm]")
    a1.set_title(f"Masa máx. de lastre factible [kg] — {relleno}", fontsize=10)
    vmax = np.nanmax(M.to_numpy()) if np.isfinite(M.to_numpy()).any() else 1
    for i, L in enumerate(Ls):
        for j, D in enumerate(Ds):
            f = c[(c["L_mm"] == L) & (c["D_mm"] == D)].iloc[0]
            if np.isfinite(f["m_relleno_max_g"]):
                txt = (f"{f['m_relleno_max_g'] / 1000:.2f} kg\nSM {f['SM_cal']:.2f}\n"
                       f"D_ac {f['D_acostado_mm']:.0f} mm\n{f['n_factibles']}/{f['n_variantes']} fact.")
            else:
                txt = f"sin factibles\n0/{f['n_variantes']}"
            oscuro = np.isfinite(f["m_relleno_max_g"]) and f["m_relleno_max_g"] / 1000 > 0.6 * vmax
            a1.text(j, i, txt, ha="center", va="center", fontsize=8, color="white" if oscuro else "#222")
    fig.colorbar(im, ax=a1, label="kg")

    # todas las variantes factibles: masa frente a D acostado (lo que exige la bahía)
    colores = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    marcas = ["o", "s", "^", "D", "v"]
    top = ranking(res[(res["relleno"] == relleno) & (res["config_id"] != "base_ork")], exigir)
    top = top.merge(res[["relleno", "config_id", "L_mm", "D_mm"]].drop_duplicates(), on=["relleno", "config_id"])
    for i, L in enumerate(Ls):
        for j, D in enumerate(Ds):
            s = top[(top["L_mm"] == L) & (top["D_mm"] == D)]
            if s.empty:
                continue
            a2.scatter(s["D_acostado_mm"], s["m_relleno_g"] / 1000, s=40, color=colores[i % 5],
                       marker=marcas[j % 5], edgecolor="white", linewidth=1.5, alpha=0.9,
                       label=f"L {L:g} · D {D:g}")
    if alto_max:
        a2.axvline(alto_max, color="#888", lw=1, ls="--")
        a2.text(alto_max, a2.get_ylim()[0], f" alto bahía {alto_max:g} mm", color="#555", fontsize=8,
                va="bottom", rotation=90)
    a2.set_xlabel("D acostado (aletas a 45°) [mm]")
    a2.set_ylabel("Masa de lastre con SM ≥ SM_min [kg]")
    a2.set_title("Variantes de aletas factibles (color = L, marcador = D)", fontsize=10)
    a2.grid(alpha=0.25)
    a2.legend(fontsize=7, ncol=3, loc="best", frameon=False)
    fig.tight_layout()
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=130)
    plt.close(fig)


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
    exigir = (cfg.raw.get("optimizacion") or {}).get("exigir") or []
    alto_max = (cfg.raw.get("envolvente") or {}).get("alto_max_mm")

    res = pd.read_csv(dir_datos / "resultados.csv")
    comp = comparar(res, exigir)
    comp.to_csv(dir_datos / "comparacion_LD.csv", index=False, float_format="%.6g")
    piv = comp.pivot_table(index=["relleno", "L_mm"], columns="D_mm", values="m_relleno_max_g")
    piv.to_csv(dir_datos / "comparacion_LD_pivote.csv", float_format="%.0f")

    pd.set_option("display.width", 200)
    for rel, c in comp.groupby("relleno", sort=False):
        print(f"\n=== {rel}  (exigido: {', '.join(exigir) or 'nada'})")
        print(c[["puesto_LD", "L_mm", "D_mm", "m_relleno_max_g", "SM_cal", "D_acostado_mm", "D_parado_mm",
                 "n_factibles", "n_variantes", "costo_relleno_usd", "config_mejor"]]
              .to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    if not a.sin_figuras:
        for rel in ([a.relleno] if a.relleno else comp["relleno"].unique()):
            _figura(res, comp, rel, exigir, alto_max, dir_figs / f"comparacion_LD__{rel}.png")
    print(f"\n→ {dir_datos / 'comparacion_LD.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
