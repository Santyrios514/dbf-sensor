"""Fase O1: T1 (paquete viejo intacto), T2 (aletas contenidas en toda la malla), T3 (misma
cavidad, masas y CG que sensor_lastre)."""

import math
import os
import subprocess
import sys

import numpy as np
import pytest

from sensor_lastre.config import cargar as cargar_lastre
from sensor_lastre.geometria import construir_cavidad
from sensor_lastre.masas import masa_vacia
from sensor_opt.aletas import construir, dimensiones, velocidad_flutter
from sensor_opt.config import ConfigError, cargar, raw_cuerpo
from sensor_opt.geometria import construir_cuerpo, theta_eq

from .conftest import RAIZ

MM = 1e-3


@pytest.mark.skipif(os.environ.get("SKIP_T1") == "1", reason="SKIP_T1=1")
def test_T1_tests_viejos_en_verde():
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "tests/"], cwd=RAIZ,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:]


def test_malla_por_defecto(cfg_default):
    assert len(cfg_default.cuerpos()) == 3 * 4 * 5 * 7 == 420
    assert len(cfg_default.aletas()) == 4 * 3 * 7 * 3 * 3 == 756
    assert cfg_default.n_evaluaciones == 420 * 756
    c = cfg_default.cotas
    assert c["D_ap_lo"] == pytest.approx(0.060) and c["D_ap_hi"] == pytest.approx(2.4 * 0.070)
    assert (c["k_lo"], c["k_hi"]) == (0.15, 1.0)  # f3 usa k efectivo, que llega a 1 con cola roma


def test_T2_aleta_contenida_en_toda_la_malla(cfg_default):
    """x_LE ≥ x_t0, la aleta no pasa de la base y P3 cae sobre la cola, en todos los candidatos.
    Las dimensiones se revisan en toda la malla; el polígono completo, en un cuerpo por forma."""
    tol = 0.01 * MM
    aletas = cfg_default.aletas()
    formas_vistas = set()
    for spec in cfg_default.cuerpos():
        cu = construir_cuerpo(cfg_default, spec)
        if cu.caso is None:
            continue
        L, x_t0 = cfg_default.L, cu.x_t0
        r = cu.r_sup()
        for a in aletas:
            d = dimensiones(cu, a, L)
            assert d.x_LE >= x_t0 - 1e-12
            assert d.x_LE + d.c_r <= L + 1e-12 and d.x_LE + d.x_s + d.c_t <= L + 1e-12
            assert abs(float(r(np.array([d.x_LE + d.c_r]))[0]) - d.r_TE_raiz) < tol
        clave = (spec.forma, spec.parametro)
        if clave in formas_vistas:
            continue
        formas_vistas.add(clave)
        for a in aletas:
            g, d, mot = construir(cfg_default, cu, a)
            if g is None:
                continue
            assert g.poligono[:, 0].min() >= x_t0 - 1e-12
            assert g.poligono[:, 0].max() <= L + 1e-12
            x3, r3 = d.x_LE + g.puntos[3, 0], d.r_LE + g.puntos[3, 1]
            assert abs(r3 - float(r(np.array([x3]))[0])) < tol
    assert len(formas_vistas) == 4


def test_T3_misma_cavidad_masas_y_cg_que_sensor_lastre(cfg_default):
    """Cuerpo de sensor_opt frente al mismo cuerpo escrito a mano en el YAML de sensor_lastre."""
    spec = next(c for c in cfg_default.cuerpos() if c.D_mm == 65 and c.forma == "parabolica"
                and c.parametro == 1.0 and c.L_t_rel_D == 1.0 and c.k == 0.6)
    cu = construir_cuerpo(cfg_default, spec)
    raw = dict(cfg_default.base_raw)
    raw["geometria_base"] = {**raw["geometria_base"],
                             "nariz": {"forma": "elipsoide", "parametro": None,
                                       "longitud": {"modo": "relativo_D", "valor": 1.0}},
                             "cola": {**raw["geometria_base"]["cola"], "forma": "parabolica", "parametro": 1.0,
                                      "longitud": {"modo": "relativo_D", "valor": 1.0},
                                      "diametro_popa": {"modo": "relativo_D", "valor": 0.6}}}
    raw["geometria_base"]["aletas"] = raw_cuerpo(cfg_default, spec)["geometria_base"]["aletas"]
    raw["barrido"] = {"parametros": {"L_mm": [400], "D_mm": [65]}}
    raw["configuraciones"] = {}
    caso = cargar_lastre(raw).casos[0]
    cav = construir_cavidad(caso)
    mv = masa_vacia(caso, cav)

    def rel(a, b):
        return abs(a - b) / abs(b)

    assert rel(cu.cav.V_int, cav.V_int) < 1e-9 and rel(cu.cav.V_ext, cav.V_ext) < 1e-9
    assert rel(cu.mv.m_casco, mv.m_casco) < 1e-9
    assert rel(cu.mv.M_casco / cu.mv.m_casco, mv.M_casco / mv.m_casco) < 1e-9
    assert np.allclose(cu.cav.r_i, cav.r_i, rtol=0, atol=1e-15)


def test_cuerpo_y_aleta_de_referencia(cfg_default):
    spec = next(c for c in cfg_default.cuerpos() if c.L_t_rel_D == 1.5)
    cu = construir_cuerpo(cfg_default, spec)
    assert cu.ok and cu.caso.geom.L == pytest.approx(0.400)
    assert cu.caso.geom.nariz.Ln == pytest.approx(spec.D) and cu.caso.geom.cola.Ra == pytest.approx(spec.k * spec.D / 2)
    assert cu.theta_eq == pytest.approx(math.atan(spec.D / 2 * (1 - spec.k) / spec.L_t))
    assert theta_eq(0.03, 1.0, 0.05) == 0.0


def test_validaciones_de_aleta_descartan_sin_error(cfg_default):
    cu = construir_cuerpo(cfg_default, cfg_default.cuerpos()[0])
    from sensor_opt.config import AletaSpec
    g, d, mot = construir(cfg_default, cu, AletaSpec(0.0, 1.0, 1.01, 0.5, 0.0))  # h ≈ 0.3 mm
    assert g is None and "h_menor_minimo" in mot
    g, d, mot = construir(cfg_default, cu, AletaSpec(0.9, 0.1, 2.0, 0.3, 0.0))  # cuerdas diminutas
    assert g is None and "c_r_menor_minimo" in mot and "c_t_menor_minimo" in mot


def test_cuerpo_con_cola_demasiado_larga_se_descarta(raw):
    raw["malla"]["L_t_rel_D"] = [6.0]
    cfg = cargar(raw)
    cu = construir_cuerpo(cfg, cfg.cuerpos()[0])
    assert not cu.ok and "L_c_no_positivo" in cu.motivos


def test_config_invalida(raw):
    raw["restricciones"]["k_min"] = 0
    raw["malla"]["lambda_LE"] = [1.0]
    with pytest.raises(ConfigError) as e:
        cargar(raw)
    assert any("k_min" in m for m in e.value.errores) and any("lambda_LE" in m for m in e.value.errores)


def test_flutter_referencia_spec():
    A = 0.5 * (0.060 + 0.030) * 0.040
    v2 = velocidad_flutter(0.040, A, 0.060, 0.030, 0.002, 2.4e9, 0.38, 0.0)
    v1 = velocidad_flutter(0.040, A, 0.060, 0.030, 0.001, 2.4e9, 0.38, 0.0)
    assert 350 < v2 < 450 and 120 < v1 < 170  # spec v3 §3.2.1: ≈ 425 y ≈ 150 m/s
    assert v2 / v1 == pytest.approx(2 ** 1.5, rel=1e-9)  # V_f ∝ t^{3/2}


def test_arrastre_de_cola_como_openrocket():
    """f_b = (3 − L_t/ΔD)/2 entre las finezas 3 y 1 (SymmetricComponentCalc, OpenRocket 24.12)."""
    from sensor_opt.geometria import fineza_cola, fraccion_base_roma, k_efectivo
    assert fineza_cola(0.070, 0.15, 0.056) == pytest.approx(0.056 / (0.070 * 0.85))
    assert [fraccion_base_roma(f) for f in (4.0, 3.0, 2.0, 1.0, 0.5)] == [0.0, 0.0, 0.5, 1.0, 1.0]
    assert k_efectivo(0.3, 0.0) == pytest.approx(0.3) and k_efectivo(0.3, 1.0) == pytest.approx(1.0)
    assert k_efectivo(0.3, 0.5) == pytest.approx(math.sqrt(0.09 + 0.5 * 0.91))


def test_cola_base_roma_se_descarta(cfg_default):
    """El ganador de 10 kg anterior (parabólica, L_t = 0.8 D, k = 0.15, θ ≈ 28°) queda fuera."""
    spec = next(c for c in cfg_default.cuerpos() if c.D_mm == 70 and c.forma == "parabolica"
                and c.parametro == 1.0 and c.L_t_rel_D == 0.8 and c.k == 0.15)
    cu = construir_cuerpo(cfg_default, spec)
    assert "cola_base_roma" in cu.motivos and cu.fineza < 1 and cu.k_ef == pytest.approx(1.0)
    suave = next(c for c in cfg_default.cuerpos() if c.D_mm == 70 and c.L_t_rel_D == 2.0 and c.k == 0.6)
    cu = construir_cuerpo(cfg_default, suave)
    assert cu.ok and cu.f_base_roma == 0.0 and cu.k_ef == pytest.approx(0.6)


def test_config_base_roma(raw):
    raw["restricciones"]["cola_base_roma"] = {"angulo_inicio_deg": 30, "angulo_base_roma_deg": 20}
    with pytest.raises(ConfigError):
        cargar(raw)
