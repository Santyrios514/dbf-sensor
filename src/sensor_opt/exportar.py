"""Salidas de sensor_opt (v3 §7): CSV, JSON y figuras. Los `.ork` los escribe `verificacion`."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from sensor_lastre.figuras import fig_dibujo, fig_perfil  # noqa: E402

from .aletas import construir  # noqa: E402
from .barrido import evaluar_cuerpo  # noqa: E402
from .config import AletaSpec, ConfigOpt, CuerpoSpec  # noqa: E402
from .esquemas import RANKING  # noqa: E402
from .geometria import Cuerpo, construir_cuerpo  # noqa: E402
from .lastre import Llenado, llenar  # noqa: E402
from .objetivo import tolerancias  # noqa: E402
from .sustituto import IDENTIDAD, Calibracion, cp_sustituto  # noqa: E402

MM = 1e-3
AZUL, NARANJA, AQUA, AMARILLO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TINTA, TINTA2, REJILLA = "#1f1f1e", "#5f5e58", "#e4e3df"


def escribir_csv(df: pd.DataFrame, ruta: Path, columnas: list[str] | None = None):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    if columnas is not None:
        df = df.reindex(columns=list(dict.fromkeys(columnas + [c for c in df.columns if c not in columnas])))
    df.to_csv(ruta, index=False, float_format="%.6g")


def escribir_json(obj: dict, ruta: Path):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float), encoding="utf-8")


def leer_ranking(ruta: Path) -> pd.DataFrame:
    df = pd.read_csv(ruta, keep_default_na=True)
    for c in ("factible", "en_pareto", "verificado_or", "flutter_margen_bajo", "refinamiento"):
        if c in df:
            df[c] = df[c].astype(str).str.lower().isin(["true", "1"])
    df["motivos"] = df["motivos"].fillna("")
    return df


def exportar_ranking(cfg: ConfigOpt, df: pd.DataFrame, dir_datos: Path):
    escribir_csv(df, dir_datos / "ranking.csv", RANKING)
    escribir_csv(df[df["en_pareto"].astype(bool)], dir_datos / "pareto.csv", RANKING)
    escribir_json(tolerancias(cfg), dir_datos / "tolerancias_implicitas.json")


# --------------------------------------------------------------------------- detalle de un candidato


@dataclass
class Detalle:
    cu: Cuerpo
    llenado: Llenado | None
    motivos: list[str]


def specs_de_fila(r) -> tuple[CuerpoSpec, AletaSpec]:
    par = None if pd.isna(r["cola_parametro"]) else float(r["cola_parametro"])
    c = CuerpoSpec(float(r["D_mm"]), r["cola_forma"], par, float(r["L_t_rel_D"]), float(r["k"]))
    a = AletaSpec(float(r["lambda_LE"]), float(r["mu_cr"]), float(r["r_tip_rel_R"]), float(r["gamma_ct"]),
                  float(r["sigma_flecha"]))
    return c, a


def detalle(cfg: ConfigOpt, c: CuerpoSpec, a: AletaSpec, cal: Calibracion = IDENTIDAD,
            x_CP: float | None = None) -> Detalle:
    """Reconstruye el candidato con todo el detalle (res, rr) para figuras y verificación.
    x_CP fuerza el CP (p. ej. el de OpenRocket) en lugar del sustituto o calibrado."""
    cu = construir_cuerpo(cfg, c)
    if not cu.ok:
        return Detalle(cu, None, cu.motivos)
    g, d, mot = construir(cfg, cu, a)
    if g is None:
        return Detalle(cu, None, mot)
    cp, m = cp_sustituto(cu, g, cu.caso.vuelo.mach, cu.caso.numerico.n_franjas_aleta, cfg.restricciones.eps_CN)
    if cp is None:
        return Detalle(cu, None, [m])
    x = float(cal.aplicar(cp.x_CP)) if x_CP is None else x_CP
    fuente = "openrocket" if x_CP is not None else ("sustituto" if cal.identidad else "calibrado")
    cp_uso = type(cp)(x_CP=x, CNa=cp.CNa, partes=cp.partes)
    Ll = llenar(cu, g, cp_uso, cfg.restricciones.SM_max, fuente)
    return Detalle(cu, Ll, Ll.motivos)


def _rr_dibujo(Ll: Llenado):
    """Las figuras de sensor_lastre dibujan el tapón delantero hasta ventana.ell_SM_max; con el
    tope de masa el tapón real es más corto (ℓ del llenado), así que se pasa una copia."""
    rr = Ll.rr
    if rr.ventana is None or not math.isfinite(Ll.ell):
        return rr
    return replace(rr, ventana=replace(rr.ventana, ell_SM_max=Ll.ell))


# --------------------------------------------------------------------------- figuras


def _estilo(ax):
    ax.grid(color=REJILLA, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_ganador(cfg: ConfigOpt, fila: pd.Series, dir_fig: Path, cal: Calibracion = IDENTIDAD,
                x_CP: float | None = None, sufijo: str = "") -> list[Path]:
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a, cal, x_CP)
    if det.llenado is None:
        return []
    dir_fig.mkdir(parents=True, exist_ok=True)
    res, rr = det.llenado.res, _rr_dibujo(det.llenado)
    res.caso = replace(res.caso, id=str(fila["cand_id"]))
    rutas = [dir_fig / f"ganador_perfil{sufijo}.png", dir_fig / f"ganador_dibujo{sufijo}.png"]
    fig_perfil(res, rr, rutas[0])
    fig_dibujo(res, rr, rutas[1])
    return rutas


def fig_pareto(df: pd.DataFrame, ruta: Path):
    p = df[df["en_pareto"].astype(bool)]
    if p.empty:
        return None
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    vmin, vmax = float(p["SM_cal"].min()), float(p["SM_cal"].max())
    sc = axs[0].scatter(p["D_ap_mm"], p["k"], c=p["SM_cal"], cmap="Blues", vmin=vmin, vmax=vmax, s=36,
                        edgecolor=TINTA2, linewidth=0.5)
    axs[0].set_xlabel("D aparente [mm]")
    axs[0].set_ylabel("k = d_a / D")
    axs[0].set_title("Frente de Pareto (color = SM)", loc="left", fontsize=10)
    fig.colorbar(sc, ax=axs[0], label="SM [cal]")
    axs[1].scatter(p["D_ap_mm"], p["SM_cal"], s=30, color=AZUL, edgecolor="white", linewidth=1)
    axs[1].set_xlabel("D aparente [mm]")
    axs[1].set_ylabel("SM [cal]")
    axs[2].scatter(p["k"], p["SM_cal"], s=30, color=NARANJA, edgecolor="white", linewidth=1)
    axs[2].set_xlabel("k = d_a / D")
    axs[2].set_ylabel("SM [cal]")
    for ax in axs:
        _estilo(ax)
    g = df[df["factible"].astype(bool)].head(1)
    if not g.empty:
        axs[0].plot(g["D_ap_mm"], g["k"], marker="*", ms=16, color=AMARILLO, mec=TINTA, label="ganador")
        axs[0].legend(frameon=False, loc="upper right")
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def fig_factibilidad(df: pd.DataFrame, ruta: Path):
    m = df.loc[~df["factible"].astype(bool), "motivos"].fillna("").str.split(";").explode()
    m = m[m != ""].value_counts().sort_values()
    if m.empty:
        return None
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(m) + 1.4), constrained_layout=True)
    ax.barh(m.index, m.values, color=AZUL, height=0.6)
    for y, v in enumerate(m.values):
        ax.text(v, y, f" {v:,}", va="center", fontsize=8, color=TINTA)
    ax.set_xlabel("candidatos descartados (un candidato puede tener varios motivos)")
    ax.set_title(f"Motivos de descarte · {int((~df['factible'].astype(bool)).sum()):,} de {len(df):,} candidatos",
                 loc="left", fontsize=10)
    _estilo(ax)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


VARIABLES = [("D_mm", "D [mm]"), ("L_t_rel_D", "L_t / D"), ("k", "k"), ("lambda_LE", "λ_LE"),
             ("mu_cr", "μ"), ("r_tip_rel_R", "r_tip / R"), ("gamma_ct", "γ"), ("sigma_flecha", "σ")]


def sensibilidad(cfg: ConfigOpt, fila: pd.Series, cal: Calibracion = IDENTIDAD) -> pd.DataFrame:
    """J y SM variando una variable a la vez (valores de la malla configurada) alrededor de `fila`."""
    from .barrido import completar
    c0, a0 = specs_de_fila(fila)
    filas = []
    for var, _ in VARIABLES:
        for v in sorted(set(float(x) for x in cfg.malla[var]) | {float(fila[var])}):
            c = replace(c0, **{var: v}) if var in ("D_mm", "L_t_rel_D", "k") else c0
            a = replace(a0, **{var: v}) if var not in ("D_mm", "L_t_rel_D", "k") else a0
            for f in evaluar_cuerpo(cfg, c, [a], cal):
                f["variable"], f["valor"] = var, v
                filas.append(f)
    return completar(cfg, pd.DataFrame(filas))


def fig_sensibilidad(cfg: ConfigOpt, fila: pd.Series, ruta: Path, cal: Calibracion = IDENTIDAD):
    s = sensibilidad(cfg, fila, cal)
    fig, axs = plt.subplots(2, len(VARIABLES), figsize=(2.1 * len(VARIABLES), 5.2), sharey="row",
                            constrained_layout=True)
    J0 = float(fila["J"])
    for i, (var, etq) in enumerate(VARIABLES):
        d = s[s["variable"] == var].sort_values("valor")
        ok = d["factible"].astype(bool)
        for ax, col, color in ((axs[0, i], "J", AZUL), (axs[1, i], "SM_cal", NARANJA)):
            y = d[col] - J0 if col == "J" else d[col]
            ax.plot(d["valor"], y, color=color, lw=1.5, zorder=1)
            ax.scatter(d.loc[ok, "valor"], y[ok], s=22, color=color, zorder=2)
            ax.scatter(d.loc[~ok, "valor"], y[~ok], s=22, facecolor="white", edgecolor=color, zorder=2)
            ax.axvline(float(fila[var]), color=TINTA2, lw=0.8, ls=":")
            _estilo(ax)
        axs[1, i].set_xlabel(etq)
    axs[0, 0].set_ylabel("J − J_ganador")
    axs[0, 0].set_yscale("symlog", linthresh=1.0)
    axs[1, 0].set_ylabel("SM [cal]")
    for ax in axs[1]:
        ax.axhspan(cfg.restricciones.SM_min, cfg.restricciones.SM_max, color=REJILLA, alpha=0.5, lw=0, zorder=0)
    fig.suptitle("Sensibilidad alrededor del ganador (punto hueco = infactible)", x=0.01, ha="left", fontsize=10)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def fig_calibracion(ver: pd.DataFrame, cal: Calibracion, ruta: Path):
    if ver.empty:
        return None
    fig, ax = plt.subplots(figsize=(5.6, 5), constrained_layout=True)
    for rol, color in (("verificacion", AZUL), ("calibracion", NARANJA)):
        d = ver[ver["rol"].str.contains(rol)]
        ax.scatter(d["x_CP_sust_mm"], d["x_CP_or_mm"], s=30, color=color, edgecolor="white", linewidth=1,
                   label={"verificacion": "mejores", "calibracion": "muestra estratificada"}[rol])
    xs = np.linspace(ver["x_CP_sust_mm"].min(), ver["x_CP_sust_mm"].max(), 50)
    ax.plot(xs, cal.alpha / MM + cal.beta * xs, color=TINTA, lw=1.2,
            label=f"x_OR = {cal.alpha / MM:.2f} + {cal.beta:.4f} x_sust  (R² = {cal.R2:.4f})")
    ax.plot(xs, xs, color=TINTA2, lw=0.8, ls="--", label="identidad")
    ax.set_xlabel("x_CP sustituto [mm]")
    ax.set_ylabel("x_CP OpenRocket [mm]")
    ax.set_title(f"Calibración del CP · residuo máx. {cal.res_max / MM:.2f} mm", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _estilo(ax)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def figuras_barrido(cfg: ConfigOpt, df: pd.DataFrame, dir_fig: Path, cal: Calibracion = IDENTIDAD) -> list[Path]:
    dir_fig.mkdir(parents=True, exist_ok=True)
    rutas = [fig_pareto(df, dir_fig / "pareto.png"), fig_factibilidad(df, dir_fig / "factibilidad.png")]
    fac = df[df["factible"].astype(bool)]
    if not fac.empty:
        g = fac.iloc[0]
        rutas += fig_ganador(cfg, g, dir_fig, cal)
        rutas.append(fig_sensibilidad(cfg, g, dir_fig / "sensibilidad.png", cal))
    return [r for r in rutas if r is not None]
