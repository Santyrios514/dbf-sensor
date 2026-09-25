import math

import pytest

from sensor_lastre.config import ConfigError, cargar, deep_merge
from sensor_lastre.atmosfera import isa

from conftest import CONFIG, solo


def test_carga_barrido_completo():
    cfg = cargar(CONFIG)
    ids = [c.id for c in cfg.casos]
    assert len(ids) == 9
    assert ids[0] == "L270_D60" and ids[-1] == "L400_D70"
    c = cfg.caso("L350_D65")
    assert c.geom.L == pytest.approx(0.350)
    assert c.geom.nariz.Ln == pytest.approx(0.065)  # relativo_D = 1.0
    assert c.geom.cola.Ra == pytest.approx(0.6 * 0.065 / 2)
    assert c.geom.aletas.xs == pytest.approx(c.geom.aletas.cr - c.geom.aletas.ct)
    assert c.geom.aletas.x_r0 == pytest.approx(0.350 - 0.25 * 0.350)


def test_rellenos_granular():
    cfg = cargar(CONFIG)
    rho = {r.nombre: r.rho_b for r in cfg.rellenos}
    assert rho["plomo_macizo"] == 11340
    assert rho["perdigon_pb_epoxy"] == pytest.approx(0.62 * 11340 + 0.38 * 1150)


def test_override_deep_merge(raw):
    raw["configuraciones"] = {"C350_65_ogiva": {"L_mm": 350, "D_mm": 65,
                                                 "nariz": {"forma": "ogiva_tangente"}}}
    cfg = cargar(raw)
    assert len(cfg.casos) == 10
    c = cfg.caso("C350_65_ogiva")
    assert c.geom.nariz.forma == "ogiva_tangente"
    assert c.geom.nariz.Ln == pytest.approx(0.065)  # conserva la longitud de la base


def test_deep_merge_no_muta():
    a = {"x": {"y": 1, "z": 2}}
    b = deep_merge(a, {"x": {"y": 5}})
    assert b == {"x": {"y": 5, "z": 2}} and a == {"x": {"y": 1, "z": 2}}


@pytest.mark.parametrize("mutar, fragmento", [
    (lambda r: r["geometria_base"]["nariz"]["longitud"].update(valor=4.0), "L_n + L_t"),
    (lambda r: r["geometria_base"]["cola"]["diametro_popa"].update(valor=1.2), "R_a"),
    (lambda r: r["pared"]["por_defecto"][0].update(material="unobtainium"), "unobtainium"),
    (lambda r: r["pared"]["por_defecto"][1].update(espesor_mm=0), "espesor_mm"),
    (lambda r: r["lastre"].update(factor_llenado=0), "factor_llenado"),
    (lambda r: r["lastre"].update(factor_llenado=1.2), "factor_llenado"),
    (lambda r: r["rellenos"][4]["granular"].update(empaquetamiento=1.0), "empaquetamiento"),
    (lambda r: r["estabilidad"].update(SM_max_cal=1.0), "SM_max_cal"),
    (lambda r: r["geometria_base"]["cola"].update(popa_cerrada=False) or
               r["geometria_base"]["cola"]["diametro_popa"].update(valor=0.01), "popa abierta"),
])
def test_validacion(raw, mutar, fragmento):
    solo(raw, 270, 60)
    mutar(raw)
    with pytest.raises(ConfigError) as e:
        cargar(raw)
    assert fragmento in str(e.value)


def test_isa_nivel_del_mar_y_medellin():
    s0 = isa(0.0)
    assert s0.rho == pytest.approx(1.225, rel=1e-3)
    assert s0.a == pytest.approx(340.29, rel=1e-3)
    s = isa(1495.0)
    assert s.rho == pytest.approx(1.0586, rel=1e-3)  # tabla ISA: 1.0581 a 1500 m
    assert math.isfinite(s.a)
