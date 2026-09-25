"""Figuras por configuración y mapas de calor L×D (sección 6.2, paso 3)."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from .analisis import ResultadoCaso, ResultadoRelleno  # noqa: E402
from .config import ESTACIONES  # noqa: E402

MM = 1e-3

# Paleta de referencia (categórica en orden fijo + rampa secuencial azul)
AZUL, NARANJA, AQUA, AMARILLO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TINTA, TINTA2, REJILLA, SUPERFICIE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SECUENCIAL = LinearSegmentedColormap.from_list(
    "azul", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
CAPAS = ["#eda100", "#f5b89c", "#a8dfc7", "#e87ba4"]

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE, "axes.facecolor": SUPERFICIE, "axes.edgecolor": TINTA2,
    "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2, "text.color": TINTA,
    "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 9,
})


def _sombrear_ventana(ax, rr: ResultadoRelleno):
    if rr.ventana is None:
        return
    for i, (a, b) in enumerate(rr.ventana.intervalos):
        ax.axvspan(a / MM, b / MM, color=AQUA, alpha=0.15, lw=0,
                   label="cumple SM" if i == 0 else None)


def fig_cg_sm(res: ResultadoCaso, rr: ResultadoRelleno, ruta: Path):
    c = rr.curva
    caso = res.caso
    D = caso.geom.D
    x_max = res.x_CP - caso.estabilidad.SM_min * D
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)

    a1.plot(c["ell_mm"], c["x_CG_mm"], color=AZUL)
    a1.axhline(res.x_CP / MM, color=NARANJA, lw=1.2, ls="--")
    a1.axhline(x_max / MM, color=TINTA2, lw=1.2, ls=":")
    a1.annotate(f"x_CP = {res.x_CP / MM:.1f} mm", (0, res.x_CP / MM), textcoords="offset points",
                xytext=(4, 4), color=TINTA2)
    a1.annotate(f"x_max (SM = {caso.estabilidad.SM_min:g}) = {x_max / MM:.1f} mm", (0, x_max / MM),
                textcoords="offset points", xytext=(4, -12), color=TINTA2)
    _sombrear_ventana(a1, rr)
    a1.set_xlabel("ℓ lastre [mm]")
    a1.set_ylabel("x_CG [mm]")
    a1.set_title("Centro de gravedad", loc="left")

    a2.plot(c["ell_mm"], c["SM_cal"], color=AZUL)
    a2.axhline(caso.estabilidad.SM_min, color=TINTA2, lw=1.2, ls=":")
    if caso.estabilidad.SM_max is not None:
        a2.axhline(caso.estabilidad.SM_max, color=TINTA2, lw=1.2, ls=":")
    _sombrear_ventana(a2, rr)
    a2.set_xlabel("ℓ lastre [mm]")
    a2.set_ylabel("SM [cal]")
    a2.set_title("Margen estático", loc="left")
    if rr.ventana is not None and rr.ventana.intervalos:
        a2.legend(frameon=False, loc="lower center")
    fig.suptitle(f"{caso.id} · {rr.relleno.nombre} (ρ_b = {rr.relleno.rho_b:.0f} kg/m³)",
                 x=0.01, ha="left", fontsize=10)
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_perfil(res: ResultadoCaso, rr: ResultadoRelleno | None, ruta: Path):
    cav, caso = res.cav, res.caso
    x = cav.x / MM
    fig, ax = plt.subplots(figsize=(10, 2.8), constrained_layout=True)
    capas_dib = set()
    for i, s in enumerate(ESTACIONES):
        m = cav.estacion == i
        fr = cav.fronteras[s]
        for k, capa in enumerate(caso.pared[s], start=1):
            lab = capa.material if capa.material not in capas_dib else None
            capas_dib.add(capa.material)
            ax.fill_between(x[m], fr[k][m] / MM, fr[k - 1][m] / MM, color=CAPAS[(k - 1) % len(CAPAS)],
                            lw=0, label=lab)
    ax.plot(x, cav.r_e / MM, color=TINTA, lw=1.0)
    ax.plot(x, cav.r_i / MM, color=TINTA2, lw=0.8)
    if rr is not None and rr.ventana is not None and math.isfinite(rr.ventana.ell_SM_max):
        a, b = res.x_b0, res.x_b0 + rr.ventana.ell_SM_max
        m = (cav.x >= a) & (cav.x <= b)
        ax.fill_between(x[m], 0, cav.r_i[m] / MM, color=AZUL, alpha=0.55, lw=0,
                        label=f"lastre útil ({rr.relleno.nombre})")
    b_geo = res.x_b0 + res.lim.ell_geo
    ax.axvline(b_geo / MM, color=NARANJA, lw=1.2, ls="--")
    ax.annotate("ℓ_geo", (b_geo / MM, caso.geom.R / MM), textcoords="offset points",
                xytext=(3, -12), color=TINTA2)
    for mp in caso.mamparos:
        ax.axvspan(mp.x / MM, (mp.x + mp.e) / MM, color=TINTA2, alpha=0.4, lw=0)
    ax.set_aspect("equal")
    ax.set_ylim(0, caso.geom.R / MM * 1.15)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(f"{caso.id} · perfil exterior, capas y región de lastre", loc="left")
    ax.legend(frameon=False, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.25))
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_mapas(resultados: pd.DataFrame, dir_fig: Path):
    """Mapa de calor L × D de la masa de relleno útil y del SM, por relleno."""
    barr = resultados[resultados["config_id"].str.match(r"^L[\d.]+_D[\d.]+$")]
    for rel, df in barr.groupby("relleno", sort=False):
        Ls = sorted(df["L_mm"].unique())
        Ds = sorted(df["D_mm"].unique())
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
        for ax, col, titulo, fmt in ((axs[0], "m_relleno_g", "Masa de relleno útil [g]", "{:.0f}"),
                                     (axs[1], "SM_cal", "SM en ℓ_SM,max [cal]", "{:.2f}")):
            Z = np.full((len(Ls), len(Ds)), np.nan)
            for _, f in df.iterrows():
                Z[Ls.index(f["L_mm"]), Ds.index(f["D_mm"])] = f[col]
            im = ax.imshow(Z, cmap=SECUENCIAL, origin="lower", aspect="auto")
            fin = Z[np.isfinite(Z)]
            medio = (np.nanmin(Z) + np.nanmax(Z)) / 2 if fin.size else 0
            for i in range(len(Ls)):
                for j in range(len(Ds)):
                    v = Z[i, j]
                    txt = "—" if not np.isfinite(v) else fmt.format(v)
                    ax.text(j, i, txt, ha="center", va="center",
                            color="white" if np.isfinite(v) and v > medio else TINTA)
            ax.set_xticks(range(len(Ds)), [f"{d:g}" for d in Ds])
            ax.set_yticks(range(len(Ls)), [f"{l:g}" for l in Ls])
            ax.set_xlabel("D [mm]")
            ax.set_ylabel("L [mm]")
            ax.set_title(titulo, loc="left")
            ax.grid(False)
            fig.colorbar(im, ax=ax, shrink=0.85)
        fig.suptitle(f"Barrido L × D · {rel}", x=0.01, ha="left", fontsize=10)
        fig.savefig(dir_fig / f"mapa__{rel}.png", dpi=150)
        plt.close(fig)
