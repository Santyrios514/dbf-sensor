"""Planos de los ganadores de la optimización (hoja A3 acotada), sin JVM.

Cada plano es una hoja A3 apaisada a escala normalizada (1:x, la mayor que cabe), con:

- vista lateral en corte: perfil, capas de pared, tapones de plomo, cavidad de electrónica, herraje
  y aletas en verdadera magnitud; CG y CP con el brazo SM·D acotado;
- vista posterior (desde popa): cuerpo, popa, aletas con su rotación y el círculo del D aparente;
- cotas en mm (L, L_n, L_c, L_t, D, d_popa, D_ap, Δx_LE, x_s, c_t, c_r, h, espesor de aleta, tapones
  de plomo y cavidad de electrónica), barra de escala, tabla de cotas, cajetín y notas.

El candidato se reconstruye desde su fila del ranking (`exportar.detalle`) con el CP de OpenRocket
cuando está verificado, así que el lastre dibujado es el del llenado validado. `plano` devuelve las
cotas tal como se dibujaron (valor y texto) para contrastarlas con el ranking. Los PNG y PDF no
llevan fecha ni software en los metadatos: la misma entrada da los mismos bytes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, Polygon, Rectangle  # noqa: E402

from sensor_lastre.config import ESTACIONES  # noqa: E402
from sensor_lastre.figuras import CAPAS, SUPERFICIE  # noqa: E402

from .config import ConfigOpt  # noqa: E402
from .exportar import AMARILLO, AQUA, AZUL, NARANJA, REJILLA, TINTA, TINTA2, Detalle, detalle, specs_de_fila  # noqa: E402

MM, G = 1e-3, 1e-3
HOJA_MM = (420.0, 297.0)  # A3 apaisado
MARCO_MM = 10.0
ESCALAS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0)  # 1:x normalizadas
ROJO = "#c8322b"
METADATOS_PNG = {"Software": None}
METADATOS_PDF = {"CreationDate": None, "Producer": None, "Creator": None}
FUENTE = 6.5  # pt de las cotas
TITULO = "Sensor remolcado · DBF 2026-27 (UPB) · plano de optimización"

# cota del plano → columna del ranking con el mismo valor (geometría; el lastre depende del CP usado)
COL_RANKING = {"L_t": "L_t_mm", "D": "D_mm", "D_ap": "D_ap_mm", "c_r": "c_r_mm", "c_t": "c_t_mm",
               "x_s": "x_s_mm", "h": "h_mm", "L_c": "L_c_mm"}

NOMBRES = {"L": "L total", "L_n": "L_n nariz", "L_c": "L_c cuerpo cilíndrico", "L_t": "L_t cola",
           "D": "D cuerpo", "d_popa": "d_popa diámetro de popa", "D_ap": "D_ap aparente",
           "H_ap": "H_ap altura aparente (caja mínima)", "W_caja": "ancho de la caja con H_ap",
           "x_LE": "x_LE borde de ataque (desde la punta)", "dx_LE": "Δx_LE desde el inicio de la cola",
           "c_r": "c_r cuerda de raíz", "c_t": "c_t cuerda de punta", "x_s": "x_s flecha", "h": "h envergadura",
           "t_aleta": "t espesor de aleta", "plomo_delantero": "tapón de plomo delantero",
           "x_electronica": "x inicio de la cavidad de electrónica", "L_electronica": "largo de la cavidad",
           "plomo_trasero": "tapón de plomo trasero", "SM_D": "brazo SM·D (CG → CP)"}


@dataclass
class Plano:
    rutas: list[Path]
    cotas: dict[str, tuple[float, str]] = field(default_factory=dict)  # nombre → (valor mm, texto dibujado)
    infactible: bool = False


# --------------------------------------------------------------------------- cotas


def _txt(v: float) -> str:
    """Cota en mm con hasta dos decimales (15.75, 46.13, 90)."""
    t = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if t in ("-0", "") else t


class _Cotas:
    """Dibuja cotas y guarda (valor, texto) de cada una."""

    def __init__(self, ax):
        self.ax = ax
        self.d: dict[str, tuple[float, str]] = {}

    def _reg(self, nombre, v, etiqueta):
        self.d[nombre] = (float(v), _txt(v))
        return f"{etiqueta} {_txt(v)}" if etiqueta else _txt(v)

    def h(self, nombre, x0, x1, y, etiqueta="", y_ext=(None, None)):
        """Cota horizontal entre x0 y x1 a la altura y; y_ext: desde dónde salen las referencias."""
        ax = self.ax
        t = self._reg(nombre, x1 - x0, etiqueta)
        for x, ye in zip((x0, x1), y_ext):
            if ye is not None:
                ax.plot([x, x], [ye, y + (1.5 if y > ye else -1.5)], color=TINTA2, lw=0.35)
        ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<|-|>", color=TINTA, lw=0.5,
                                                          mutation_scale=5, shrinkA=0, shrinkB=0))
        return ax.text((x0 + x1) / 2, y + 0.8, t, ha="center", va="bottom", fontsize=FUENTE, color=TINTA,
                       bbox=dict(fc="white", ec="none", pad=0.3))

    def v(self, nombre, x, y0, y1, etiqueta="", x_ext=(None, None), izquierda=True):
        ax = self.ax
        t = self._reg(nombre, y1 - y0, etiqueta)
        for y, xe in zip((y0, y1), x_ext):
            if xe is not None:
                ax.plot([xe, x + (-1.5 if x < xe else 1.5)], [y, y], color=TINTA2, lw=0.35)
        ax.annotate("", (x, y0), (x, y1), arrowprops=dict(arrowstyle="<|-|>", color=TINTA, lw=0.5,
                                                          mutation_scale=5, shrinkA=0, shrinkB=0))
        ax.text(x + (-0.8 if izquierda else 0.8), (y0 + y1) / 2, t, ha="right" if izquierda else "left",
                va="center", rotation=90, fontsize=FUENTE, color=TINTA, bbox=dict(fc="white", ec="none", pad=0.3))

    def valor(self, nombre, v):
        """Cota anotada como texto (no como línea)."""
        self.d[nombre] = (float(v), _txt(v))
        return _txt(v)


# --------------------------------------------------------------------------- estado del llenado


@dataclass
class Estado:
    """Lo que se dibuja del lastre: tapones, electrónica, CG y CP."""

    ell: float
    ell2: float
    x_r0: float
    x_e: float
    x_CG: float
    x_CP: float
    m_total: float
    factible: bool
    fila: dict


def estado(det: Detalle) -> Estado | None:
    Ll = det.llenado
    if Ll is None or not math.isfinite(Ll.ell):
        return None
    f, mod = Ll.fila, Ll.rr.modelo
    ell, ell2 = Ll.ell, float(f.get("ell_trasero_mm", 0.0) or 0.0) * MM
    return Estado(ell, ell2, float(mod.x_r0(ell)), float(mod.x_e(ell)), f["x_CG_mm"] * MM, Ll.res.x_CP,
                  f["m_total_g"] * G, Ll.ok, f)


# --------------------------------------------------------------------------- vistas


def _vista_lateral(ax, det: Detalle, est: Estado | None, cot: _Cotas, Rmax: float):
    caso = det.llenado.res.caso if det.llenado is not None else det.cu.caso
    g, cav, a = caso.geom, det.cu.cav, caso.geom.aletas
    x = cav.x / MM
    R, Ra, L, Ln, x_t0 = g.R / MM, g.cola.Ra / MM, g.L / MM, g.nariz.Ln / MM, g.x_cola / MM
    re = lambda xx: float(np.interp(xx, x, cav.r_e / MM))  # noqa: E731
    # aletas en verdadera magnitud (arriba y abajo)
    if a is not None:
        P = a.poligono / MM
        for s in (1, -1):
            ax.add_patch(Polygon(np.c_[P[:, 0], s * P[:, 1]], closed=True, fc=AMARILLO, alpha=0.55, ec=TINTA,
                                 lw=0.6, zorder=1))
    # pared por capas
    for i, s in enumerate(ESTACIONES):
        m = cav.estacion == i
        fr = cav.fronteras[s]
        for k in range(1, len(caso.pared[s]) + 1):
            for sg in (1, -1):
                ax.fill_between(x[m], sg * fr[k][m] / MM, sg * fr[k - 1][m] / MM, color=CAPAS[(k - 1) % len(CAPAS)],
                                lw=0, zorder=2)
    for sg in (1, -1):
        ax.plot(x, sg * cav.r_e / MM, color=TINTA, lw=0.8, zorder=3)
        ax.plot(x, sg * cav.r_i / MM, color=TINTA2, lw=0.3, zorder=3)
    for xu in (Ln, x_t0):  # uniones nariz–cuerpo y cuerpo–cola
        ax.plot([xu, xu], [-re(xu), re(xu)], color=TINTA2, lw=0.4, ls=(0, (4, 2)), zorder=3)
    ax.plot([-6, L + 6], [0, 0], color=TINTA2, lw=0.4, ls=(0, (10, 2, 2, 2)), zorder=3)  # eje
    # lastre y cavidad de electrónica
    Le = caso.electronica.Le / MM
    if est is not None:
        tramos = [(det.cu.x_b0, det.cu.x_b0 + est.ell)] + ([(est.x_r0, est.x_r0 + est.ell2)] if est.ell2 > 0 else [])
        for t0, t1 in tramos:
            m = (cav.x >= t0) & (cav.x <= t1)
            ax.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AZUL, alpha=0.55, lw=0, hatch="////",
                            ec="white", zorder=2.5)
        m = (cav.x >= est.x_e) & (cav.x <= est.x_e + caso.electronica.Le)
        ax.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AQUA, alpha=0.6, lw=0, zorder=2.5)
        ax.text(est.x_e / MM + Le / 2, R * 0.5, "electrónica", rotation=90, ha="center", va="center", fontsize=5.5,
                color=TINTA, zorder=6)  # en la mitad superior: el CP suele caer dentro de la cavidad
    for pm in caso.puntuales:
        ax.add_patch(Rectangle((pm.x / MM - 2, -2), 4, 4, fc=TINTA, ec="none", zorder=6))
        ax.annotate(pm.nombre.replace("_", " "), (pm.x / MM, -2), xytext=(pm.x / MM, -R * 0.55),
                    fontsize=5.5, ha="center", color=TINTA, arrowprops=dict(arrowstyle="-", lw=0.4, color=TINTA2),
                    zorder=6, bbox=dict(fc="white", ec="none", pad=0.2))
    # CG, CP y brazo SM·D
    if est is not None:
        xcg, xcp = est.x_CG / MM, est.x_CP / MM
        ax.plot(xcg, 0, "o", ms=6, mfc="white", mec=TINTA, mew=1.0, zorder=7)
        ax.plot(xcg, 0, marker=(2, 0, 45), ms=6, color=TINTA, mew=0.8, zorder=7)
        ax.text(xcg, -3, f"CG {xcg:.1f}", ha="center", va="top", fontsize=FUENTE, zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.3))
        ax.plot(xcp, 0, "D", ms=5, mfc=NARANJA, mec=TINTA, mew=0.6, zorder=7)
        ax.text(xcp, -3, f"CP {xcp:.1f}", ha="center", va="top", fontsize=FUENTE, zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.3))
        t = cot.h("SM_D", xcg, xcp, 6.0, "SM·D =")
        t.set_text(f"SM·D = {_txt(xcp - xcg)} ({(est.x_CP - est.x_CG) / g.D:.2f} cal)")
    # cotas del cuerpo (abajo)
    y1, y2 = -(Rmax + 10), -(Rmax + 20)
    cot.h("L_n", 0, Ln, y1, "L_n", (0, -re(Ln)))
    cot.h("L_c", Ln, x_t0, y1, "L_c", (None, -re(x_t0)))
    cot.h("L_t", x_t0, L, y1, "L_t", (None, -re(L)))
    cot.h("L", 0, L, y2, "L", (y1, y1))
    cot.v("D", -14, -R, R, "D", (Ln, Ln))
    cot.v("d_popa", L + 24, -Ra, Ra, "d_popa", (L, L), izquierda=False)
    # aleta (arriba)
    if a is not None:
        P = a.puntos / MM
        xLE, rLE = a.x_LE / MM, a.r_LE / MM
        cr, xs, ct = a.params.c_r / MM, P[1, 0], P[2, 0] - P[1, 0]
        r_tip = a.r_tip / MM
        ya, yb = r_tip + 7, r_tip + 15
        cot.valor("x_LE", xLE)
        if xLE - x_t0 > 1e-6:
            cot.h("dx_LE", x_t0, xLE, ya, "Δx_LE", (re(x_t0), rLE))
        else:
            cot.valor("dx_LE", 0.0)
        if xs > 1e-6:
            cot.h("x_s", xLE, xLE + xs, ya, "x_s", (None if xLE - x_t0 > 1e-6 else rLE, r_tip))
        else:
            cot.valor("x_s", 0.0)
        cot.h("c_t", xLE + xs, xLE + xs + ct, ya, "c_t", (None if xs > 1e-6 else r_tip, r_tip))
        cot.h("c_r", xLE, xLE + cr, yb, "c_r", (rLE, re(xLE + cr)))
        cot.v("h", L + 10, rLE, r_tip, "h", (xLE, xLE + xs + ct), izquierda=False)
    # lastre y cavidad (arriba del cuerpo)
    if est is not None:
        yp = R + 7
        b0 = det.cu.x_b0 / MM
        cot.h("plomo_delantero", b0, b0 + est.ell / MM, yp, "plomo", (re(b0), re(b0 + est.ell / MM)))
        xe = est.x_e / MM
        t = cot.h("electronica", xe, xe + Le, R + 15, "", (re(xe), re(xe + Le)))
        t.set_text(f"electrónica {_txt(Le)} @ x = {cot.valor('x_electronica', xe)}")
        if est.ell2 > 0:
            a0 = est.x_r0 / MM
            t = cot.h("plomo_trasero", a0, a0 + est.ell2 / MM, yp, "", (re(a0), re(a0 + est.ell2 / MM)))
            t.set_text(f"plomo {_txt(est.ell2 / MM)}")
            t.set_position((a0 + est.ell2 / MM + 2, yp + 0.8))  # tapón corto: el texto va a su derecha
            t.set_ha("left")
        else:
            cot.valor("plomo_trasero", 0.0)
    # barra de escala
    yb_ = -(Rmax + 31)
    for i in range(5):
        ax.add_patch(Rectangle((i * 10, yb_), 10, 2.2, fc=TINTA if i % 2 == 0 else "white", ec=TINTA, lw=0.4))
    for v in (0, 50):
        ax.text(v, yb_ - 1, f"{v}", ha="center", va="top", fontsize=FUENTE)
    ax.text(52, yb_ + 1.1, "mm", ha="left", va="center", fontsize=FUENTE)


def _puntas(r_tip: float, t: float, n: int, phi0: float) -> np.ndarray:
    """Esquinas de las puntas de las n aletas (rectángulos de espesor t) en la sección transversal,
    con φ medido desde la vertical: (x, y) = r_tip·(sin φ, cos φ) ± (t/2)·(cos φ, −sin φ)."""
    phi = phi0 + 2 * np.pi * np.arange(n) / n
    u = np.c_[np.sin(phi), np.cos(phi)]
    nrm = np.c_[u[:, 1], -u[:, 0]]
    return np.r_[r_tip * u + t / 2 * nrm, r_tip * u - t / 2 * nrm]


def _alto_ancho(R: float, r_tip: float, t: float, n: int, phi0: float) -> tuple[float, float]:
    P = _puntas(r_tip, t, n, phi0)
    alto = max(R, P[:, 1].max()) + max(R, -P[:, 1].min())
    ancho = max(R, P[:, 0].max()) + max(R, -P[:, 0].min())
    return float(alto), float(ancho)


def altura_aparente(R: float, r_tip: float, t: float, n: int, n_ang: int = 3601) -> tuple[float, float, float]:
    """(H_ap, ancho, φ*): mínima altura de la caja que contiene la sección (cuerpo + aletas con su
    espesor) con el sensor acostado, girándolo sobre su eje. Basta barrer φ en [0, 2π/n): el
    conjunto de aletas se repite con ese periodo. Con 4 aletas φ* = 45° y H = √2 (r_tip + t/2)."""
    from scipy.optimize import minimize_scalar
    per = 2 * np.pi / n
    phis = np.linspace(0.0, per, n_ang)
    alts = np.array([_alto_ancho(R, r_tip, t, n, p)[0] for p in phis])
    i = int(np.argmin(alts))
    lo, hi = phis[max(i - 1, 0)], phis[min(i + 1, n_ang - 1)]
    r = minimize_scalar(lambda p: _alto_ancho(R, r_tip, t, n, p)[0], bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-10})
    phi = float(r.x) if r.fun <= alts[i] else float(phis[i])
    H, W = _alto_ancho(R, r_tip, t, n, phi)
    return H, W, phi % per


def _vista_posterior(ax, det: Detalle, cot: _Cotas, lim: float):
    caso = det.llenado.res.caso if det.llenado is not None else det.cu.caso
    g, a = caso.geom, caso.geom.aletas
    R, Ra = g.R / MM, g.cola.Ra / MM
    ax.add_patch(Circle((0, 0), R, fc=REJILLA, ec=TINTA, lw=0.8))
    ax.add_patch(Circle((0, 0), Ra, fc=SUPERFICIE, ec=TINTA, lw=0.6))
    ax.plot([-lim * 0.9, lim * 0.9], [0, 0], color=TINTA2, lw=0.3, ls=(0, (10, 2, 2, 2)))
    ax.plot([0, 0], [-lim * 0.9, lim * 0.9], color=TINTA2, lw=0.3, ls=(0, (10, 2, 2, 2)))
    if a is None:
        return
    pa = a.params
    t, r_tip = pa.t / MM, a.r_tip / MM
    r_raiz = float(a.poligono[:, 1].min()) / MM  # la raíz sigue la cola: desde su radio más bajo
    for phi in pa.rotacion + 2 * np.pi * np.arange(pa.n) / pa.n:
        u = np.array([np.sin(phi), np.cos(phi)])
        nrm = np.array([u[1], -u[0]])
        cc = np.array([r_raiz * u + t / 2 * nrm, r_tip * u + t / 2 * nrm, r_tip * u - t / 2 * nrm,
                       r_raiz * u - t / 2 * nrm])
        ax.add_patch(Polygon(cc, closed=True, fc=AMARILLO, ec=TINTA, lw=0.5, zorder=3))
    D_ap = 2 * max(R, r_tip)
    ax.add_patch(Circle((0, 0), D_ap / 2, fc="none", ec=NARANJA, lw=0.8, ls=(0, (5, 3))))
    ax.text(0, D_ap / 2 + 2, f"D_ap = {cot.valor('D_ap', D_ap)}", ha="center", va="bottom", fontsize=FUENTE,
            color=NARANJA)
    ax.text(0, -lim + 2, f"{pa.n} aletas a {math.degrees(pa.rotacion):g}° · t = {cot.valor('t_aleta', t)} · "
            f"r_tip = {r_tip:.1f}", ha="center", va="bottom", fontsize=FUENTE)
    # altura aparente: caja mínima con el sensor acostado (aletas giradas a φ*)
    H, W, phi = altura_aparente(R, r_tip, t, pa.n)
    per = 2 * np.pi / pa.n
    d = (phi - pa.rotacion) % per
    if min(d, per - d) > 1e-6:  # la orientación de vuelo no es la de la caja: aletas fantasma a φ*
        for q in _puntas(r_tip, t, pa.n, phi)[:pa.n]:
            ax.plot([0, q[0]], [0, q[1]], color=AZUL, lw=0.5, ls=(0, (2, 2)))
    ax.add_patch(Rectangle((-W / 2, -H / 2), W, H, fc="none", ec=AZUL, lw=0.8, ls=(0, (6, 2)), zorder=4))
    cot.v("H_ap", -(lim - 5), -H / 2, H / 2, "H_ap", (-W / 2, -W / 2))
    cot.valor("W_caja", W)
    ax.text(0, H / 2 + 1.5, f"caja mínima (aletas a {math.degrees(phi):.0f}°)", ha="center", va="bottom",
            fontsize=5.5, color=AZUL)


# --------------------------------------------------------------------------- cajetín


def _filas_cajetin(fila: pd.Series, det: Detalle, est: Estado | None, puesto, val: dict | None,
                   escala: float) -> list[tuple[str, str]]:
    f = est.fila if est is not None else {}
    caso = det.llenado.res.caso if det.llenado is not None else det.cu.caso
    g = caso.geom
    m_rel = f.get("m_relleno_g", math.nan)
    m_tras = f.get("m_relleno_trasero_g", 0.0) or 0.0
    SM = (est.x_CP - est.x_CG) / g.D if est is not None else math.nan
    banderas = [b for b in (f.get("banderas") or "").split(";") if b]
    if bool(fila.get("flutter_margen_bajo", False)):
        banderas.append("flutter_margen_bajo")
    if val is None:
        linea_val = "sin validar (CP del sustituto)"
    else:
        ok = str(val.get("factible_or", "")).lower() in ("true", "1")
        linea_val = (f"OpenRocket: x_CP {val['x_CP_or_mm']:.2f} mm (Δ vs calibrado {val['dx_or_cal_mm']:+.2f}), "
                     f"SM {val['SM_or_cal']:.3f}, masa Δ {val['dif_masa_or_pct']:+.2f} % · "
                     + ("VALIDADO" if ok else "NO VALIDADO"))
    estado_txt = "factible" if est is not None and est.factible else "INFACTIBLE"
    return [
        ("Candidato", str(fila["cand_id"])),
        ("Puesto (verificado con OR)", "—" if puesto is None else f"{int(puesto)} · {estado_txt}"),
        ("Masa total", f"{(est.m_total / G if est else math.nan):.1f} g"),
        ("  casco / aletas", f"{f.get('m_casco_g', math.nan):.1f} g / {f.get('m_aletas_g', math.nan):.1f} g"),
        ("  plomo delantero / trasero", f"{m_rel - m_tras:.1f} g / {m_tras:.1f} g"),
        ("  electrónica / herraje", f"{f.get('m_electronica_g', math.nan):.1f} g / "
                                    f"{f.get('m_puntuales_g', math.nan):.1f} g"),
        ("x_CG / x_CP (desde la punta)", f"{(est.x_CG / MM if est else math.nan):.1f} / "
                                         f"{(est.x_CP / MM if est else math.nan):.1f} mm"),
        ("SM · C_Nα", f"{SM:.3f} cal · {f.get('CN_alpha_rad', math.nan):.3f}"),
        ("Tolerancia de amarre", f"{f.get('tol_amarre_mm', math.nan):.2f} mm"),
        ("θ_eq cola · fineza · f_b", f"{fila['theta_eq_deg']:.1f}° · {fila['fineza_cola']:.2f} · "
                                     f"{fila['f_base_roma']:.2f}"),
        ("V_flutter", f"{fila.get('V_flutter_m_s', math.nan):.0f} m/s"),
        ("Restricción activa", str(f.get("limitante_masa", "")) or "—"),
        ("Banderas", ", ".join(banderas) or "ninguna"),
        ("Validación", linea_val),
        ("Escala · unidades", f"1:{escala:g} en A3 · cotas en mm"),
    ]


def _cajetin(ax, filas: list[tuple[str, str]], titulo: str):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, fc="white", ec=TINTA, lw=0.8))
    n = len(filas) + 1
    alto = 1 / n
    ax.add_patch(Rectangle((0, 1 - alto), 1, alto, fc=REJILLA, ec=TINTA, lw=0.6))
    ax.text(0.015, 1 - alto / 2, titulo, va="center", fontsize=7.5, weight="bold")
    for i, (k, v) in enumerate(filas, start=1):
        y = 1 - (i + 0.5) * alto
        ax.plot([0, 1], [1 - i * alto] * 2, color=TINTA2, lw=0.3)
        ax.text(0.015, y, k, va="center", fontsize=6.2, color=TINTA2)
        ax.text(0.33, y, v, va="center", fontsize=6.2, color=TINTA)
    ax.plot([0.32, 0.32], [0, 1 - alto], color=TINTA2, lw=0.3)


# --------------------------------------------------------------------------- plano


def _escala(ancho_mm: float, alto_mm: float, ancho_disp: float, alto_disp: float) -> float:
    for e in ESCALAS:
        if ancho_mm / e <= ancho_disp and alto_mm / e <= alto_disp:
            return e
    return ESCALAS[-1]


def plano(cfg: ConfigOpt, fila: pd.Series, ruta_base: Path, puesto=None, validacion: dict | None = None,
          dpi: int = 300) -> Plano:
    """Dibuja el plano de `fila` (una fila del ranking) en ruta_base.png y ruta_base.pdf. Con
    `validacion` (fila de verificacion_or.csv) se usa el CP de OpenRocket para el llenado."""
    c, a_spec = specs_de_fila(fila)
    x_CP = None if validacion is None else float(validacion["x_CP_or_mm"]) * MM
    det = detalle(cfg, c, a_spec, x_CP=x_CP)
    if not det.cu.ok or det.llenado is None:
        raise ValueError(f"{fila['cand_id']}: no se puede reconstruir ({';'.join(det.motivos)})")
    est = estado(det)
    caso = det.llenado.res.caso
    g, a = caso.geom, caso.geom.aletas
    L, R = g.L / MM, g.R / MM
    Rmax = max(R, a.r_tip / MM if a is not None else R)
    lat = (-34.0, L + 40.0, -(Rmax + 37.0), Rmax + 24.0)
    lim = Rmax + 14.0
    W, H = HOJA_MM
    esc = _escala((lat[1] - lat[0]) + 2 * lim, max(lat[3] - lat[2], 2 * lim), W - 2 * MARCO_MM - 26.0,
                  H - 2 * MARCO_MM - 120.0)
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    fig.patch.set_facecolor("white")

    def ejes(x0, y0, w, h):
        return fig.add_axes([x0 / W, y0 / H, w / W, h / H])

    marco = ejes(0, 0, W, H)
    marco.set_xlim(0, W)
    marco.set_ylim(0, H)
    marco.axis("off")
    marco.add_patch(Rectangle((MARCO_MM, MARCO_MM), W - 2 * MARCO_MM, H - 2 * MARCO_MM, fc="none", ec=TINTA, lw=1.0))
    w_lat, h_lat = (lat[1] - lat[0]) / esc, (lat[3] - lat[2]) / esc
    y_vistas = H - MARCO_MM - 14.0 - max(h_lat, 2 * lim / esc)
    al = ejes(MARCO_MM + 8.0, y_vistas, w_lat, h_lat)
    al.set_xlim(lat[0], lat[1])
    al.set_ylim(lat[2], lat[3])
    al.set_aspect("equal")
    al.axis("off")
    cot = _Cotas(al)
    cot.valor("L_electronica", caso.electronica.Le / MM)
    _vista_lateral(al, det, est, cot, Rmax)
    x_post = MARCO_MM + 8.0 + w_lat + 10.0
    ap = ejes(x_post, y_vistas + (h_lat - 2 * lim / esc) / 2, 2 * lim / esc, 2 * lim / esc)
    ap.set_xlim(-lim, lim)
    ap.set_ylim(-lim, lim)
    ap.set_aspect("equal")
    ap.axis("off")
    cot.ax = ap
    _vista_posterior(ap, det, cot, lim)
    y_tit = y_vistas + max(h_lat, 2 * lim / esc) + 3.0
    marco.text(MARCO_MM + 8.0, y_tit, "VISTA LATERAL EN CORTE (aletas en verdadera magnitud)", fontsize=8, weight="bold")
    marco.text(x_post, y_tit, "VISTA POSTERIOR (desde popa)", fontsize=8, weight="bold")
    infactible = est is None or not est.factible
    if infactible:
        marco.text(MARCO_MM + 8.0 + w_lat / 2, y_vistas + h_lat / 2, "INFACTIBLE", fontsize=46, color=ROJO,
                   alpha=0.35, ha="center", va="center", rotation=12, weight="bold", zorder=20)
        marco.text(MARCO_MM + 8.0, y_vistas - 4.0, f"INFACTIBLE · {';'.join(det.motivos) or 'restricción'}",
                   fontsize=9, color=ROJO, weight="bold", va="top")
    # tabla de cotas (los mismos textos que en las vistas)
    items = [(NOMBRES[k], cot.d[k][1]) for k in NOMBRES if k in cot.d]
    y_tab = y_vistas - 12.0
    marco.text(MARCO_MM + 6.0, y_tab, "COTAS [mm]", fontsize=7.5, weight="bold", va="top")
    por_col = math.ceil(len(items) / 3)
    for i, (k, v) in enumerate(items):
        xc = MARCO_MM + 6.0 + (i // por_col) * 66.0
        yc = y_tab - 6.0 - (i % por_col) * 4.6
        marco.text(xc, yc, k, fontsize=6.4, color=TINTA2, va="top")
        marco.text(xc + 60.0, yc, v, fontsize=6.4, color=TINTA, va="top", ha="right")
    # cajetín y notas
    filas = _filas_cajetin(fila, det, est, puesto, validacion, esc)
    ancho_caj, alto_caj = 205.0, 92.0
    caj = ejes(W - MARCO_MM - ancho_caj, MARCO_MM, ancho_caj, alto_caj)
    _cajetin(caj, filas, TITULO)
    pared = " + ".join(f"{cp.material} {cp.t / MM:.1f}" for cp in caso.pared["cuerpo"])
    el = caso.electronica
    cav_txt = (f", cono D {el.D_cav[0] / MM:g} → {el.D_cav[1] / MM:g} recortado por la pared" if el.D_cav else "")
    rho_b = det.llenado.rr.relleno.rho_b if hasattr(det.llenado.rr.relleno, "rho_b") else math.nan
    notas = [
        "NOTAS",
        "1. Cotas en mm; x desde la punta de la nariz. Escala válida al imprimir el PDF en A3 al 100 %.",
        f"2. Pared: {pared} mm (afuera → adentro). Plomo macizo ({rho_b:.0f} kg/m³), rayado.",
        f"3. Aletas: {a.params.n} × {a.params.material}, espesor {a.params.t / MM:.1f} mm, raíz sobre la superficie "
        f"de la cola." if a is not None else "3. Sin aletas.",
        f"4. Cavidad de electrónica de {el.Le / MM:g} mm en la junta cuerpo–cola{cav_txt}; sin plomo en ese tramo.",
        "5. CG y CP del modelo propio; con validación, CP de OpenRocket (ver cajetín).",
        "6. H_ap: mínima altura de caja con el sensor acostado, girado a su mejor ángulo; incluye el espesor "
        "de las aletas, sin holgura.",
    ]
    for i, t in enumerate(notas):
        marco.text(MARCO_MM + 6.0, MARCO_MM + alto_caj - 4.0 - i * 5.0, t, fontsize=6.8 if i else 7.5,
                   weight="bold" if i == 0 else "normal", va="top")
    ruta_base.parent.mkdir(parents=True, exist_ok=True)
    # sin with_suffix: el cand_id tiene puntos (r2.2) que Path tomaría como extensión
    rutas = [ruta_base.parent / f"{ruta_base.name}.png", ruta_base.parent / f"{ruta_base.name}.pdf"]
    fig.savefig(rutas[0], dpi=dpi, metadata=METADATOS_PNG)
    fig.savefig(rutas[1], metadata=METADATOS_PDF)
    plt.close(fig)
    return Plano(rutas=rutas, cotas=cot.d, infactible=infactible)


# --------------------------------------------------------------------------- comparativo


def comparativo(cfg: ConfigOpt, filas: list[tuple[str, pd.Series]], ruta: Path, dpi: int = 200) -> Path | None:
    """Siluetas laterales superpuestas a la misma escala, con masa y D aparente en la leyenda."""
    if not filas:
        return None
    fig, ax = plt.subplots(figsize=(14, 4.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for i, (etq, f) in enumerate(filas):
        c, a = specs_de_fila(f)
        det = detalle(cfg, c, a)
        cav = det.cu.cav
        col = cmap(i / max(len(filas) - 1, 1))
        x, r = cav.x / MM, cav.r_e / MM
        m = f.get("m_total_or_g", f["m_total_g"])
        lab = f"{etq} · {m / 1000:.2f} kg · D_ap {f['D_ap_mm']:.1f} mm · {f['cola_forma']}"
        ax.plot(np.r_[x, x[::-1]], np.r_[r, -r[::-1]], color=col, lw=1.2, label=lab)
        if det.llenado is not None and det.llenado.res.caso.geom.aletas is not None:
            P = det.llenado.res.caso.geom.aletas.poligono / MM
            for s in (1, -1):
                ax.plot(np.r_[P[:, 0], P[0, 0]], s * np.r_[P[:, 1], P[0, 1]], color=col, lw=0.9)
    ax.axhline(0, color=TINTA2, lw=0.4, ls=(0, (10, 2, 2, 2)))
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm] (desde la punta)")
    ax.set_ylabel("r [mm]")
    ax.set_title("Siluetas de los ganadores (misma escala; aletas en verdadera magnitud)", loc="left", fontsize=10)
    ax.grid(color=REJILLA, lw=0.5)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=dpi, metadata=METADATOS_PNG)
    plt.close(fig)
    return ruta


def perfil_csv(cfg: ConfigOpt, fila: pd.Series) -> pd.DataFrame:
    """x, r_e, r_i cada 0.5 mm y el polígono de la aleta (x, r), para CAD/CFD."""
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a)
    cav = det.cu.cav
    xs = np.arange(0.0, cav.x[-1] + 1e-12, 0.5 * MM)
    df = pd.DataFrame({"tipo": "perfil", "x_mm": xs / MM, "r_e_mm": np.interp(xs, cav.x, cav.r_e) / MM,
                       "r_i_mm": np.interp(xs, cav.x, cav.r_i) / MM})
    if det.llenado is not None and det.llenado.res.caso.geom.aletas is not None:
        P = det.llenado.res.caso.geom.aletas.poligono / MM
        df = pd.concat([df, pd.DataFrame({"tipo": "aleta", "x_mm": P[:, 0], "r_e_mm": P[:, 1]})], ignore_index=True)
    return df
