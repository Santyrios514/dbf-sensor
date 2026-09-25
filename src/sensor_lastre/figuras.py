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
        ell2 = rr.fila.get("ell_trasero_mm", 0.0) * MM
        if ell2 > 0:
            a2 = rr.modelo.x_r0(rr.ventana.ell_SM_max)
            m2 = (cav.x >= a2) & (cav.x <= a2 + ell2)
            ax.fill_between(x[m2], 0, cav.r_i[m2] / MM, color=AZUL, alpha=0.55, lw=0)
        el = caso.electronica
        x_e = float(rr.modelo.x_e(rr.ventana.ell_SM_max))
        ax.axvspan(x_e / MM, (x_e + el.Le) / MM, ymax=0.3, color=AQUA, alpha=0.35, lw=0,
                   label="electrónica")
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
    resto = [r for r in rutas if r not in (px, ps)]  # parámetros extra: una figura por combinación
    df = resultados.merge(valores, on="config_id")
    grupos_fig = df.groupby(resto, sort=True) if resto else [((), df)]
    for clave, dfr in grupos_fig:
        clave = clave if isinstance(clave, tuple) else (clave,)
        sufijo = "".join(f"__{_etiqueta(r).split('.')[-1]}{v}" for r, v in zip(resto, clave))
        _fig_barrido_grupo(dfr, px, ps, SM_min, dir_fig, sufijo)


def _fig_barrido_grupo(df, px, ps, SM_min, dir_fig, sufijo):
    for rel, d in df.groupby("relleno", sort=False):
        for col, archivo, titulo, ref in (
            ("SM_inf_cal", "sm_inf", "SM_∞ (techo por geometría) [cal]", SM_min),
            ("SM_max_alcanzable_cal", "sm_max", "SM máximo alcanzable con este relleno [cal]", SM_min),
            ("m_lastre_SM_min_g", "masa_SM_min", f"Masa de lastre mínima para SM = {SM_min:g} cal [g]", None),
            ("m_relleno_g", "masa_max", f"Masa MÁXIMA de lastre con SM ≥ {SM_min:g} cal [g]", None),
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
            fig.suptitle(f"{rel}{sufijo.replace('__', ' · ')} · marcador hueco: no cabe en la bahía (alto)",
                         x=0.01, ha="left",
                         fontsize=9, color=TINTA2)
            fig.savefig(dir_fig / f"{archivo}__{rel}{sufijo}.png", dpi=150)
            plt.close(fig)


# --------------------------------------------------------------------------- dibujo de la configuración


def _cota_h(ax, x0, x1, y, texto, color=TINTA):
    ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<->", color=color, lw=0.9,
                                                      shrinkA=0, shrinkB=0))
    ax.text((x0 + x1) / 2, y, texto, ha="center", va="bottom", fontsize=8, color=color,
            bbox=dict(fc=SUPERFICIE, ec="none", pad=0.5))


def _cota_v(ax, x, y0, y1, texto, color=TINTA, lado="left"):
    ax.annotate("", (x, y0), (x, y1), arrowprops=dict(arrowstyle="<->", color=color, lw=0.9,
                                                      shrinkA=0, shrinkB=0))
    ax.text(x + (-2 if lado == "left" else 2), (y0 + y1) / 2, texto, ha="right" if lado == "left" else "left",
            va="center", fontsize=8, color=color, rotation=90)


def fig_dibujo(res: ResultadoCaso, rr: ResultadoRelleno | None, ruta: Path):
    """Dibujo acotado de la configuración: vista lateral (cuerpo, aletas proyectadas según su
    rotación, lastre delantero y trasero, electrónica, CG y CP) y vista frontal (desde la nariz)
    con el diámetro máximo con aletas y el rectángulo de bahía que exige."""
    caso, cav = res.caso, res.cav
    g, a = caso.geom, caso.geom.aletas
    pa = a.params
    x = cav.x / MM
    R, r_tip, Lt = g.R / MM, a.r_tip / MM, g.L_total / MM
    phis = pa.rotacion + 2 * np.pi * np.arange(pa.n) / pa.n  # desde la vertical
    D_max = 2 * max(g.R, a.r_tip) / MM
    from .aletas import envolvente
    h_env, w_env = (v / MM for v in envolvente(a.r_tip, g.R, pa.n, pa.rotacion))

    fig, (al, af) = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [3.3, 1]},
                                 constrained_layout=True)

    # --- vista lateral: aletas (proyección r·cos φ sobre el plano vertical), detrás del cuerpo
    P = a.poligono / MM
    for phi in phis:
        c = np.cos(phi)
        if abs(c) < 1e-6:  # aleta de canto: solo asoma detrás de la base
            al.plot([P[:, 0].min(), P[:, 0].max()], [0, 0], color=TINTA2, lw=1.2, zorder=0.5)
            continue
        al.fill(P[:, 0], P[:, 1] * c, color=AMARILLO, alpha=0.45, lw=0, zorder=0.5)
        al.plot(np.r_[P[:, 0], P[0, 0]], np.r_[P[:, 1], P[0, 1]] * c, color=TINTA2, lw=0.8, zorder=0.5)
    al.fill_between(x, -cav.r_e / MM, cav.r_e / MM, color="#e4e3df", lw=0)
    al.plot(x, cav.r_e / MM, color=TINTA, lw=1.1)
    al.plot(x, -cav.r_e / MM, color=TINTA, lw=1.1)
    al.plot(x, cav.r_i / MM, color=TINTA2, lw=0.5)
    al.plot(x, -cav.r_i / MM, color=TINTA2, lw=0.5)

    titulo = f"{caso.id}"
    if rr is not None and rr.ventana is not None and np.isfinite(rr.ventana.ell_SM_max):
        mod, ell1 = rr.modelo, rr.ventana.ell_SM_max
        ell2 = rr.fila.get("ell_trasero_mm", 0.0) * MM
        tramos = [(res.x_b0, res.x_b0 + ell1)]
        if ell2 > 0:
            tramos.append((mod.x_r0(ell1), mod.x_r0(ell1) + ell2))
        for k, (t0, t1) in enumerate(tramos):
            m = (cav.x >= t0) & (cav.x <= t1)
            al.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AZUL, alpha=0.7, lw=0,
                            label=f"lastre ({rr.relleno.nombre})" if k == 0 else None)
        x_e = float(mod.x_e(ell1))
        m = (cav.x >= x_e) & (cav.x <= x_e + caso.electronica.Le)
        al.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AQUA, alpha=0.55, lw=0,
                        label="electrónica")
        f = rr.fila
        xcg = f["x_CG_mm"]
        al.plot(xcg, 0, "o", ms=9, mfc=SUPERFICIE, mec=TINTA, mew=1.5, zorder=5)
        al.plot(xcg, 0, "o", ms=3, color=TINTA, zorder=6)
        al.annotate(f"CG {xcg:.0f}", (xcg, 0), textcoords="offset points", xytext=(0, -16), ha="center",
                    fontsize=8, color=TINTA)
        costo = f.get("costo_relleno_usd", np.nan)
        titulo += (f" · lastre {f['m_relleno_g'] / 1000:.2f} kg de {rr.relleno.nombre}"
                   f" (≈ US$ {costo:,.0f}) · SM {f['SM_cal']:.2f}")
    xcp = res.x_CP / MM
    al.plot(xcp, 0, "D", ms=7, mfc=NARANJA, mec=TINTA, mew=0.8, zorder=5)
    al.annotate(f"CP {xcp:.0f}", (xcp, 0), textcoords="offset points", xytext=(0, 9), ha="center",
                fontsize=8, color=TINTA)

    ymax = max(R, r_tip) * 1.25
    _cota_h(al, 0, Lt, -ymax * 0.95, f"L total = {Lt:.1f} mm (cuerpo {g.L / MM:.0f})")
    _cota_v(al, g.nariz.Ln / MM + 20, -R, R, f"D = {g.D / MM:.0f}")
    al.set_xlim(-8, Lt + 8)
    al.set_ylim(-ymax * 1.08, ymax)
    al.set_aspect("equal")
    al.grid(False)
    al.set_xlabel("x [mm] (desde la punta)")
    al.set_title(f"Vista lateral · aletas a {np.degrees(pa.rotacion):.0f}°", loc="left", fontsize=9)
    al.legend(frameon=False, loc="upper left", fontsize=8, ncol=2)

    # --- vista frontal (desde la nariz): aletas como placas de espesor t, cuerpo encima
    t = pa.t / MM
    r_in = pa.r_in / MM
    for phi in phis:
        u = np.array([np.sin(phi), np.cos(phi)])  # dirección radial (x horizontal, y vertical)
        n = np.array([u[1], -u[0]])
        c = [r_in * u + t / 2 * n, r_tip * u + t / 2 * n, r_tip * u - t / 2 * n, r_in * u - t / 2 * n]
        c = np.array(c)
        af.fill(c[:, 0], c[:, 1], color=AMARILLO, alpha=0.8, ec=TINTA2, lw=0.8)
    af.add_patch(plt.Circle((0, 0), R, fc="#e4e3df", ec=TINTA, lw=1.1))
    af.add_patch(plt.Circle((0, 0), g.cola.Ra / MM, fc="none", ec=TINTA2, lw=0.6, ls=":"))
    af.add_patch(plt.Circle((0, 0), D_max / 2, fc="none", ec=NARANJA, lw=1.0, ls="--"))
    af.add_patch(plt.Rectangle((-w_env / 2, -h_env / 2), w_env, h_env, fc="none", ec=AZUL, lw=1.2, ls="--"))
    lim = max(D_max, h_env, w_env) / 2 * 1.32
    _cota_h(af, -w_env / 2, w_env / 2, -h_env / 2 - lim * 0.12, f"ancho {w_env:.1f}", AZUL)
    _cota_v(af, -w_env / 2 - lim * 0.1, -h_env / 2, h_env / 2, f"alto {h_env:.1f}", AZUL)
    af.text(0, lim * 0.93, f"D máx. con aletas = {D_max:.1f} mm", ha="center", va="top", fontsize=8,
            color=NARANJA, bbox=dict(fc=SUPERFICIE, ec="none", pad=0.5))
    af.set_xlim(-lim, lim)
    af.set_ylim(-lim, lim)
    af.set_aspect("equal")
    af.grid(False)
    af.set_xticks([])
    af.set_yticks([])
    for s in af.spines.values():
        s.set_visible(False)
    af.set_title("Vista frontal · bahía (azul)", loc="left", fontsize=9)

    fig.suptitle(titulo, x=0.01, ha="left", fontsize=10)
    fig.savefig(ruta, dpi=160)
    plt.close(fig)
