"""Herraje de remolque en el CG (masa puntual `x: {modo: en_CG}`) y su alojamiento en el plomo."""

from __future__ import annotations

import math

import numpy as np
import pytest

from conftest import solo
from sensor_lastre.analisis import analizar_caso
from sensor_lastre.config import cargar
from sensor_lastre.estabilidad import ModeloLastre

MM, G = 1e-3, 1e-3
RHO_PB = 11340.0


def _caso(raw, herraje: dict):
    raw["masas_puntuales"] = [{"nombre": "herraje_remolque", **herraje}]
    return cargar(solo(raw, 350, 65)).casos[0]


def _modelo(caso):
    res = analizar_caso(caso, rellenos=[caso.rellenos[0]])
    return res, ModeloLastre(caso, res.cav, res.mv, res.x_b0, rho_b=RHO_PB)


def test_masa_en_el_CG_no_mueve_el_CG(raw):
    sin = _caso(raw, {"masa_g": 0, "x": {"modo": "en_CG"}})
    con = _caso(raw, {"masa_g": 50, "x": {"modo": "en_CG"}})
    (_, m0), (_, m1) = _modelo(sin), _modelo(con)
    ells = np.linspace(0.02, 0.2, 7)
    np.testing.assert_allclose(m1.x_CG(ells), m0.x_CG(ells), rtol=0, atol=1e-12)
    np.testing.assert_allclose(m1.m(ells) - m0.m(ells), 0.050, rtol=0, atol=1e-12)
    np.testing.assert_allclose(m1.x_CG_con_trasero(0.1, 0.01), m0.x_CG_con_trasero(0.1, 0.01), atol=1e-12)


def test_equivale_a_fijar_la_pieza_en_el_CG_final(raw):
    con = _caso(raw, {"masa_g": 50, "x": {"modo": "en_CG"}})
    _, m1 = _modelo(con)
    ell = 0.15
    x_cg = float(m1.x_CG(ell))
    fija = _caso(raw, {"masa_g": 50, "x": {"modo": "absoluto_mm", "valor": x_cg / MM}})
    _, m2 = _modelo(fija)
    assert float(m2.x_CG(ell)) == pytest.approx(x_cg, abs=1e-9)
    assert float(m2.m(ell)) == pytest.approx(float(m1.m(ell)), abs=1e-12)


def test_alojamiento_descuenta_plomo_sin_mover_el_CG(raw):
    sin = _caso(raw, {"masa_g": 15, "x": {"modo": "en_CG"}})
    con = _caso(raw, {"masa_g": 15, "x": {"modo": "en_CG"}, "alojamiento": {"diametro_mm": 10, "largo_mm": 30}})
    V = math.pi * 0.010**2 / 4 * 0.030
    assert con.puntuales[0].V_aloj == pytest.approx(V)
    (_, m0), (_, m1) = _modelo(sin), _modelo(con)
    assert m1.m_aloj == pytest.approx(RHO_PB * V)
    ells = np.linspace(0.02, 0.2, 7)
    np.testing.assert_allclose(m1.x_CG(ells), m0.x_CG(ells), atol=1e-12)
    np.testing.assert_allclose(m0.m(ells) - m1.m(ells), RHO_PB * V, atol=1e-12)


def test_largo_auto_es_el_radio_interior_del_cilindro(raw):
    caso = _caso(raw, {"masa_g": 15, "x": {"modo": "en_CG"}, "alojamiento": {"diametro_mm": 10, "largo_mm": "auto"}})
    r_i = 0.065 / 2 - sum(c.t for c in caso.pared["cuerpo"])
    assert caso.puntuales[0].V_aloj == pytest.approx(math.pi * 0.010**2 / 4 * r_i)


def test_fila_reporta_plomo_neto(raw):
    caso = _caso(raw, {"masa_g": 15, "x": {"modo": "en_CG"}, "alojamiento": {"diametro_mm": 10, "largo_mm": 30}})
    res = analizar_caso(caso, rellenos=[caso.rellenos[0]])
    f = res.rellenos[0].fila
    assert f["m_lastre_alojamiento_g"] == pytest.approx(RHO_PB * math.pi * 0.010**2 / 4 * 0.030 / G)
    m_bruto = RHO_PB * f["V_util_cm3"] * 1e-6 / G
    assert f["m_relleno_g"] == pytest.approx(m_bruto - f["m_lastre_alojamiento_g"])
    m_vacio = f["m_casco_g"] + f["m_aletas_g"] + f["m_mamparos_g"] + f["m_electronica_g"] + f["m_puntuales_g"]
    assert f["m_total_g"] == pytest.approx(m_vacio + f["m_relleno_g"], abs=1e-6)
