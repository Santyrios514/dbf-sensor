"""Función objetivo ponderada, tolerancias implícitas y frente de Pareto (v3 §3.8, §4.4).

    J = Σ w_j f_j,  f_j ∈ [0, 1] normalizados con cotas a priori de la malla configurada:

    f1 = (m_max − m_total)/m_max
    f2 = (D_ap − D_lo)/(D_hi − D_lo),       D_lo = min D, D_hi = 2 max r_tip
    f3 = (k − k_lo)/(k_hi − k_lo)
    f4 = |SM − SM_c| / ((SM_max − SM_min)/2)     (modo 'centro')
         (SM_max − SM)/(SM_max − SM_min)          (modo 'max')

Tolerancia implícita: como cada f_i ∈ [0, 1], los criterios inferiores juntos pueden compensar
en J a lo sumo Σ_{i>j} w_i. Una diferencia en f_j mayor que

    ε_j = Σ_{i>j} w_i / w_j = (w_{j+1}/w_j)(1 + w_{j+2}/w_{j+1} + …)

decide el orden sin importar los demás criterios (es la cota exacta de la aproximación de la spec).
"""

from __future__ import annotations

import numpy as np

from .config import ConfigOpt

NOMBRES = ("f1_masa", "f2_diametro", "f3_taper", "f4_SM")


def f_valores(cfg: ConfigOpt, m_total, D_ap, k, SM) -> np.ndarray:
    """Matriz (n, 4) de objetivos normalizados; cada columna recortada a [0, 1]."""
    c, r, ob = cfg.cotas, cfg.restricciones, cfg.objetivo
    m_total, D_ap, k, SM = (np.asarray(v, dtype=float) for v in (m_total, D_ap, k, SM))
    f1 = (cfg.m_max - m_total) / cfg.m_max
    f2 = (D_ap - c["D_ap_lo"]) / (c["D_ap_hi"] - c["D_ap_lo"])
    f3 = (k - c["k_lo"]) / (c["k_hi"] - c["k_lo"]) if c["k_hi"] > c["k_lo"] else np.zeros_like(k)
    if ob.f4_modo == "centro":
        f4 = np.abs(SM - ob.SM_centro) / ((r.SM_max - r.SM_min) / 2)
    else:
        f4 = (r.SM_max - SM) / (r.SM_max - r.SM_min)
    return np.clip(np.column_stack([f1, f2, f3, f4]), 0.0, 1.0)


def J(cfg: ConfigOpt, F: np.ndarray) -> np.ndarray:
    return np.asarray(F, dtype=float) @ np.asarray(cfg.objetivo.pesos, dtype=float)


def tolerancias(cfg: ConfigOpt) -> dict:
    """ε_j adimensionales y en unidades físicas (g, mm, k, cal)."""
    w = np.asarray(cfg.objetivo.pesos, dtype=float)
    eps = [float(w[j + 1:].sum() / w[j]) for j in range(4)]
    c, r, ob = cfg.cotas, cfg.restricciones, cfg.objetivo
    escala_SM = (r.SM_max - r.SM_min) / 2 if ob.f4_modo == "centro" else (r.SM_max - r.SM_min)
    fis = {
        "masa_g": eps[0] * cfg.m_max * 1e3,
        "D_ap_mm": eps[1] * (c["D_ap_hi"] - c["D_ap_lo"]) * 1e3,
        "k": eps[2] * (c["k_hi"] - c["k_lo"]),
        "SM_cal": eps[3] * escala_SM,
    }
    textos = [
        f"Un candidato con hasta {fis['masa_g']:.1f} g menos de masa total puede ganar si es mejor en "
        "diámetro aparente, taper y SM.",
        f"Con la misma masa, uno con hasta {fis['D_ap_mm']:.2f} mm más de diámetro aparente puede ganar "
        "si tiene mejor taper y SM.",
        f"Con la misma masa y diámetro, uno con hasta {fis['k']:.4f} más de k puede ganar si su SM "
        "está más cerca del objetivo.",
        "El SM es el último criterio: no compensa nada.",
    ]
    return {"pesos": w.tolist(), "eps": dict(zip(NOMBRES, eps)), "fisicas": fis,
            "cotas": {"D_ap_lo_mm": c["D_ap_lo"] * 1e3, "D_ap_hi_mm": c["D_ap_hi"] * 1e3,
                      "k_lo": c["k_lo"], "k_hi": c["k_hi"], "m_max_g": cfg.m_max * 1e3},
            "lectura": textos}


def pareto(F: np.ndarray) -> np.ndarray:
    """Máscara de los no dominados (minimización en todas las columnas)."""
    F = np.asarray(F, dtype=float)
    n = len(F)
    if n == 0:
        return np.zeros(0, dtype=bool)
    orden = np.lexsort(F.T[::-1])  # por la primera columna, luego las demás
    nd = np.zeros(n, dtype=bool)
    frente = np.empty((0, F.shape[1]))
    for i in orden:
        f = F[i]
        if frente.size and np.any(np.all(frente <= f, axis=1) & np.any(frente < f, axis=1)):
            continue
        if frente.size and np.any(np.all(frente == f, axis=1)):
            nd[i] = True  # duplicado de un punto del frente: también no dominado
            continue
        nd[i] = True
        frente = np.vstack([frente, f])
    return nd
