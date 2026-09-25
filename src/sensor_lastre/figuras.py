"""Figuras por configuración y mapas de calor L×D (sección 6.2, paso 3)."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .analisis import ResultadoCaso, ResultadoRelleno  # noqa: E402
from .config import ESTACIONES  # noqa: E402

MM = 1e-3

# Paleta de referencia (categórica en orden fijo + rampa secuencial azul)
AZUL, NARANJA, AQUA, AMARILLO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TINTA, TINTA2, REJILLA, SUPERFICIE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
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
    """Perfil exterior, capas, región de lastre útil y polígono de la aleta superpuesto."""
    cav, caso = res.cav, res.caso
    a = caso.geom.aletas
    x = cav.x / MM
    fig, ax = plt.subplots(figsize=(10, 3.4), constrained_layout=True)
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
    if a is not None:
        P = np.vstack([a.poligono, a.poligono[:1]]) / MM
        ax.fill(P[:, 0], P[:, 1], color=AMARILLO, alpha=0.35, lw=0, label="aleta")
        ax.plot(P[:, 0], P[:, 1], color=TINTA2, lw=0.8)
    if rr is not None and rr.ventana is not None and math.isfinite(rr.ventana.ell_SM_max):
        a0, b0 = res.x_b0, res.x_b0 + rr.ventana.ell_SM_max
        m = (cav.x >= a0) & (cav.x <= b0)
        ax.fill_between(x[m], 0, cav.r_i[m] / MM, color=AZUL, alpha=0.55, lw=0,
                        label=f"lastre útil ({rr.relleno.nombre})")
    b_geo = res.x_b0 + res.lim.ell_geo
    ax.axvline(b_geo / MM, color=NARANJA, lw=1.2, ls="--")
    ax.annotate("ℓ_geo", (b_geo / MM, caso.geom.R / MM), textcoords="offset points",
                xytext=(3, 4), color=TINTA2)
    ax.axvline(res.x_CP / MM, color=NARANJA, lw=1.2, ls=":")
    ax.annotate("x_CP", (res.x_CP / MM, caso.geom.R / MM), textcoords="offset points",
                xytext=(3, 4), color=TINTA2)
    for mp in caso.mamparos:
        ax.axvspan(mp.x / MM, (mp.x + mp.e) / MM, color=TINTA2, alpha=0.4, lw=0)
    ax.set_aspect("equal")
    r_max = max(caso.geom.R, a.r_tip if a is not None else 0.0)
    ax.set_ylim(0, r_max / MM * 1.1)
    ax.set_xlim(0, caso.geom.L_total / MM * 1.01)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(f"{caso.id} · perfil, capas, aleta y región de lastre", loc="left")
    ax.legend(frameon=False, loc="upper center", ncol=5, bbox_to_anchor=(0.5, -0.22))
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


SERIES = [AZUL, NARANJA, AQUA, AMARILLO]


def _eje_x(rutas: list[str]) -> tuple[str, str | None]:
    """Parámetro del eje x (el que contiene 'h_tip', si existe) y el de las series."""
    x = next((r for r in rutas if "h_tip" in r), rutas[0])
    resto = [r for r in rutas if r != x]
    return x, (resto[0] if resto else None)


def _etiqueta(ruta: str) -> str:
    return ".".join(p for p in ruta.split(".") if p != "valor")


def fig_barrido(resultados: pd.DataFrame, valores: pd.DataFrame, rutas: list[str], SM_min: float,
                dir_fig: Path):
    """Curvas de diseño del barrido (v2 §6.6): SM_∞ y masa de lastre necesaria para SM_min
    frente a h_tip, una serie por valor del segundo parámetro. Marcador hueco = no cabe en alto."""
    if not rutas:
        return
    px, ps = _eje_x(rutas)
    df = resultados.merge(valores, on="config_id")
    for rel, d in df.groupby("relleno", sort=False):
        for col, archivo, titulo, ref in (
            ("SM_inf_cal", "sm_inf", "SM_∞ (techo por geometría) [cal]", SM_min),
            ("SM_max_alcanzable_cal", "sm_max", "SM máximo alcanzable con este relleno [cal]", SM_min),
            ("m_lastre_SM_min_g", "masa_SM_min", f"Masa de lastre para SM = {SM_min:g} cal [g]", None),
        ):
            if col not in d:
                continue
            fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
            grupos = d.groupby(ps, sort=True) if ps else [(None, d)]
            for k, (val, g) in enumerate(grupos):
                g = g.sort_values(px)
                c = SERIES[k % len(SERIES)]
                ax.plot(g[px], g[col], color=c, lw=2)
                cabe = g["cabe_alto"].fillna(True).astype(bool).to_numpy()
                ax.scatter(g[px][cabe], g[col][cabe], s=36, color=c, zorder=3)
                ax.scatter(g[px][~cabe], g[col][~cabe], s=36, facecolor=SUPERFICIE, edgecolor=c,
                           lw=1.5, zorder=3)
                if ps and len(g):
                    ult = g.iloc[-1]
                    if np.isfinite(ult[col]):
                        # desfase vertical por serie: evita que etiquetas de puntos iguales se encimen
                        ax.annotate(f"{_etiqueta(ps)} = {val:g}", (ult[px], ult[col]),
                                    textcoords="offset points", xytext=(6, 11 * k - 11), va="center",
                                    color=TINTA2)
            if ref is not None:
                ax.axhline(ref, color=TINTA2, lw=1.2, ls=":")
                ax.annotate(f"SM_min = {ref:g}", (d[px].min(), ref), textcoords="offset points",
                            xytext=(0, 4), color=TINTA2)
            ax.set_xlabel(f"{_etiqueta(px)} [mm]")
            ax.set_title(titulo, loc="left")
            ax.margins(x=0.18)
            fig.suptitle(f"{rel} · marcador hueco: no cabe en la bahía (alto)", x=0.01, ha="left",
                         fontsize=9, color=TINTA2)
            fig.savefig(dir_fig / f"{archivo}__{rel}.png", dpi=150)
            plt.close(fig)
