"""Electrónica en la junta cuerpo–cola (modo junta_cola, cavidad cónica recortada por la pared)."""

from __future__ import annotations

import math

import yaml

from conftest import CONFIG_REPO, solo
from sensor_lastre.analisis import analizar_caso
from sensor_lastre.config import cargar
from sensor_lastre.estabilidad import ModeloLastre, cabe_cavidad_electronica

MM = 1e-3


def _junta(raw: dict, **cambios) -> dict:
    with open(CONFIG_REPO, encoding="utf-8") as f:
        raw["electronica"] = {**yaml.safe_load(f)["electronica"], **cambios}
    raw["lastre"].update({"fraccion_max_L": None, "masa_max_sensor_g": None})
    return raw


def test_posicion_relativa_a_la_junta(raw):
    caso = cargar(solo(_junta(raw), 350, 65)).casos[0]
    el, g = caso.electronica, caso.geom
    assert el.modo == "fija" and el.junta_cola
    assert math.isclose(el.x_inicio, g.x_cola - 10 * MM)
    assert math.isclose(el.x_inicio + el.Le, g.x_cola + 10 * MM)
    assert math.isclose(el.x_cg(el.x_inicio), g.x_cola)


def test_la_posicion_sigue_a_la_junta(raw):
    c1 = cargar(solo(_junta(raw), 350, 65)).casos[0]
    c2 = cargar(solo(_junta(raw, desfase_junta_mm=-15), 350, 65)).casos[0]
    assert math.isclose(c1.electronica.x_inicio - c2.electronica.x_inicio, 5 * MM)


def test_tapones_a_ambos_lados_de_la_cavidad(raw):
    caso = cargar(solo(_junta(raw), 350, 65)).casos[0]
    res = analizar_caso(caso, rellenos=[caso.rellenos[0]])
    el = caso.electronica
    assert res.lim.limitante == "electronica_fija"
    assert math.isclose(res.x_b0 + res.lim.ell_geo, el.x_inicio)
    # el tapón trasero arranca al final de la cavidad, ya dentro de la cola
    mod = ModeloLastre(caso, res.cav, res.mv, res.x_b0, rho_b=11340.0)
    assert math.isclose(mod.x_r0(res.lim.ell_geo), el.x_inicio + el.Le)
    assert mod.x_r0(0.0) > caso.geom.x_cola and mod.ell2_max(res.lim.ell_geo) > 0


def test_cono_recortado_es_informativo(raw):
    caso = cargar(solo(_junta(raw), 350, 65)).casos[0]
    res = analizar_caso(caso, rellenos=[caso.rellenos[0]])
    assert cabe_cavidad_electronica(caso, res.cav) is False  # D_fin = 550 mm > interior
    assert "cavidad_electronica_recortada" in res.banderas
    assert "inviable_geo" not in res.banderas


def test_cono_dentro_de_la_pared(raw):
    caso = cargar(solo(_junta(raw, cavidad={"D_inicio_mm": 40, "D_fin_mm": 40}), 350, 65)).casos[0]
    res = analizar_caso(caso, rellenos=[caso.rellenos[0]])
    assert cabe_cavidad_electronica(caso, res.cav) is True
    assert "cavidad_electronica_recortada" not in res.banderas
