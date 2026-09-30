"""Fase O2: T4 (guarda del CP), T5 (llenado), T6 (objetivo lexicográfico), T8 (calibración)."""

import math

import numpy as np
import pytest

from sensor_lastre.barrowman import Contribucion
from sensor_opt.aletas import construir
from sensor_opt.config import AletaSpec, cargar
from sensor_opt.geometria import construir_cuerpo
from sensor_opt.lastre import llenar
from sensor_opt.objetivo import J, f_valores, pareto, tolerancias
from sensor_opt.sustituto import MOTIVO_CN, IDENTIDAD, ajustar, combinar, cp_sustituto

G = 1e-3


def _cuerpo(cfg, D=65, forma="parabolica", param=1.0, lt=1.0, k=0.3):
    spec = next(c for c in cfg.cuerpos() if c.D_mm == D and c.forma == forma and c.parametro == param
                and c.L_t_rel_D == lt and c.k == k)
    return construir_cuerpo(cfg, spec)


def _candidato(cfg, cu, a):
    g, d, mot = construir(cfg, cu, a)
    assert not mot, mot
    cp, m = cp_sustituto(cu, g, cu.caso.vuelo.mach, 48, cfg.restricciones.eps_CN)
    assert m is None
    return g, cp


# --------------------------------------------------------------------------- T4


@pytest.mark.parametrize("k, CNa_f", [(0.0, 0.0), (1e-4, 1e-6), (0.05, 1e-3), (0.0, 0.3)])
def test_T4_guarda_cola_cerrada_y_aletas_minimas(k, CNa_f):
    """Nariz (+2) + cola (2(k² − 1)) + aletas ≈ 0: se descarta, nunca inf ni NaN."""
    partes = (Contribucion("nariz", 2.0, 0.03), Contribucion("cola", 2 * (k**2 - 1), 0.37),
              Contribucion("aletas", CNa_f, 0.38))
    res, motivo = combinar(partes, eps_CN=0.5)
    assert res is None and motivo == MOTIVO_CN


def test_T4_guarda_no_descarta_con_aletas():
    partes = (Contribucion("nariz", 2.0, 0.03), Contribucion("cola", -2.0, 0.37), Contribucion("aletas", 3.0, 0.38))
    res, motivo = combinar(partes, eps_CN=0.5)
    assert motivo is None and math.isfinite(res.x_CP) and res.CNa == pytest.approx(3.0)
    assert res.x_CP == pytest.approx((2 * 0.03 - 2 * 0.37 + 3 * 0.38) / 3)  # finito con k = 0


def test_T4_candidato_real_con_aleta_diminuta_se_descarta(cfg_default):
    """Con la guarda muy alta, un candidato real se descarta por CN_total_pequeno sin NaN."""
    cu = _cuerpo(cfg_default, k=0.15)
    g, d, mot = construir(cfg_default, cu, AletaSpec(0.6, 0.5, 1.2, 0.5, 0.0))
    assert not mot
    cp, motivo = cp_sustituto(cu, g, cu.caso.vuelo.mach, 48, eps_CN=1e3)
    assert cp is None and motivo == MOTIVO_CN


def test_nariz_mas_cola_suman_2k2(cfg_default):
    cu = _cuerpo(cfg_default, k=0.3)
    assert sum(p.CNa for p in cu.partes_cuerpo) == pytest.approx(2 * 0.3**2, rel=1e-3)


# --------------------------------------------------------------------------- T5


def test_T5_masa_activa_lastre_contiguo_y_m_max(raw):
    raw["masa"]["m_max_g"] = 2000
    cfg = cargar(raw)
    cu = _cuerpo(cfg, D=70, lt=1.0, k=0.3)
    g, cp = _candidato(cfg, cu, AletaSpec(0.0, 1.0, 2.0, 0.5, 1.0))
    L = llenar(cu, g, cp, cfg.restricciones.SM_max)
    f = L.fila
    assert L.restriccion_activa == "masa"
    assert f["m_total_g"] == pytest.approx(2000, abs=0.1)
    assert f["ell_trasero_mm"] == 0 and f["m_relleno_trasero_g"] == 0  # contiguo desde la nariz
    assert L.ell > 0 and L.ell < cu.lim.ell_geo
    assert L.m_delantero == pytest.approx(f["m_relleno_g"] * G)
    V = L.rr.modelo.V_b(L.ell)  # el tapón empieza en x_b0 y ocupa [x_b0, x_b0 + ℓ]
    assert L.rr.modelo.rho_b * V == pytest.approx(L.m_delantero, abs=0.1 * G)  # ℓ con tol_raiz_mm


def test_T5_SM_activo_da_SM_min(raw):
    raw["masa"]["m_max_g"] = 1e6
    cfg = cargar(raw)
    cfg.base_raw["remolque"]["tol_amarre_min_mm"] = None  # aísla la restricción de SM
    cu = _cuerpo(cfg, D=65, lt=1.0, k=0.3)
    encontrado = False
    for a in cfg.aletas():
        g, d, mot = construir(cfg, cu, a)
        if mot:
            continue
        cp, m = cp_sustituto(cu, g, cu.caso.vuelo.mach, 48, cfg.restricciones.eps_CN)
        if m:
            continue
        L = llenar(cu, g, cp, cfg.restricciones.SM_max)
        if L.ok and L.restriccion_activa == "SM_min":
            assert L.fila["SM_cal"] == pytest.approx(cfg.restricciones.SM_min, abs=0.005)
            encontrado = True
    assert encontrado


def test_sobreestabilidad_es_infactible_sin_retrasar_el_lastre(raw):
    raw["restricciones"]["SM_max_cal"] = 1.2
    cfg = cargar(raw)
    cu = _cuerpo(cfg, D=70, lt=1.0, k=0.3)
    g, cp = _candidato(cfg, cu, AletaSpec(0.0, 1.0, 2.4, 0.7, 1.0))
    L = llenar(cu, g, cp, cfg.restricciones.SM_max)
    assert "SM_sobre_max" in L.motivos and L.fila["SM_cal"] > 1.2
    assert L.fila["ell_trasero_mm"] == 0 or L.restriccion_activa != "SM_min"


# --------------------------------------------------------------------------- T6 y objetivo


def test_T6_orden_lexicografico(cfg_default):
    rng = np.random.default_rng(0)
    tol = tolerancias(cfg_default)["eps"]
    eps = list(tol.values())
    for j in range(3):
        FA = rng.uniform(0, 1, (4000, 4))
        FB = rng.uniform(0, 1, (4000, 4))
        FB[:, :j] = FA[:, :j]  # empatan en los criterios superiores
        delta = rng.uniform(1.0001, 1.5, 4000) * eps[j]
        FA[:, j] = rng.uniform(0, 1 - delta)
        FB[:, j] = FA[:, j] + delta
        FA[:, j + 1:] = 1.0  # el peor caso para A en los inferiores
        FB[:, j + 1:] = 0.0
        assert np.all(J(cfg_default, FA) < J(cfg_default, FB))


def test_tolerancias_implicitas_en_unidades(cfg_default):
    t = tolerancias(cfg_default)
    assert t["eps"]["f1_masa"] == pytest.approx(0.010101, rel=1e-4)
    assert t["fisicas"]["masa_g"] == pytest.approx(0.010101 * 6200, rel=1e-4)
    assert t["fisicas"]["D_ap_mm"] == pytest.approx(0.010101 * (168 - 60), rel=1e-4)
    assert t["fisicas"]["SM_cal"] == 0.0


def test_f_valores(cfg_default):
    F = f_valores(cfg_default, [6.2, 3.1], [0.060, 0.168], [0.15, 1.0], [1.5, 2.0])  # k efectivo
    assert np.allclose(F[0], [0, 0, 0, 0]) and np.allclose(F[1], [0.5, 1, 1, 1])
    F = f_valores(cfg_default, [6.2], [0.060], [0.6], [1.5])
    assert F[0, 2] == pytest.approx((0.6 - 0.15) / 0.85)


def test_pareto():
    F = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0], [1, 1, 1], [0.5, 0.5, 0.5], [0, 1, 1]])
    assert pareto(F).tolist() == [True, True, True, False, True, True]


# --------------------------------------------------------------------------- T8


def test_T8_calibracion_recupera_alpha_beta():
    rng = np.random.default_rng(12345)
    x = rng.uniform(0.15, 0.30, 30)
    y = 0.0123 + 0.94 * x + rng.normal(0, 1e-5, 30)
    c = ajustar(x, y)
    assert c.alpha == pytest.approx(0.0123, rel=0.01) and c.beta == pytest.approx(0.94, rel=0.01)
    assert c.R2 > 0.999 and c.res_max < 1e-4
    assert IDENTIDAD.aplicar(0.2) == 0.2
