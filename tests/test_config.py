import math

import pytest

from sensor_lastre.atmosfera import isa
from sensor_lastre.config import ConfigError, cargar, deep_merge

from conftest import CONFIG, barrido, solo


def test_carga_barrido_v2():
    cfg = cargar(CONFIG)
    ids = [c.id for c in cfg.casos]
    assert len(ids) == 4 * 3 + 1
    assert ids[0] == "L350_D65_h_tip20_extension0" and ids[-1] == "base_ork"
    c = cfg.caso("L350_D65_h_tip30_extension20")
    assert c.barrido == {"aletas.h_tip.valor": 30, "aletas.extension.valor": 20}
    assert c.params_aleta.h == pytest.approx(0.030)
    assert c.params_aleta.e == pytest.approx(0.020)
    assert c.geom.L == pytest.approx(0.350)
    assert c.geom.L_total == pytest.approx(0.350 + 0.020 + 0.02e-3)
    b = cfg.caso("base_ork")
    assert b.geom.nariz.Ln == pytest.approx(0.065)
    assert b.geom.cola.forma == "parabolica" and b.geom.cola.parametro == 1.0
    assert b.geom.cola.Ra == pytest.approx(0.0195)
    assert b.geom.cola.popa == "abierta"


def test_mach_auto():
    c = cargar(CONFIG).casos[0]
    assert c.vuelo.mach == pytest.approx(30.0 / isa(1495).a)
    assert c.vuelo.mach == pytest.approx(0.089, abs=1e-3)
    assert c.vuelo.mach_regresion_or == 0.3


def test_rellenos_granular():
    cfg = cargar(CONFIG)
    rho = {r.nombre: r.rho_b for r in cfg.rellenos}
    assert rho["plomo_macizo"] == 11340
    assert rho["perdigon_pb_epoxy"] == pytest.approx(0.62 * 11340 + 0.38 * 1150)


def test_override_deep_merge(raw):
    solo(raw, 350, 65)
    raw["configuraciones"] = {"ogiva": {"L_mm": 350, "D_mm": 65, "nariz": {"forma": "ogiva_tangente"}}}
    cfg = cargar(raw)
    assert len(cfg.casos) == 2
    c = cfg.caso("ogiva")
    assert c.geom.nariz.forma == "ogiva_tangente"
    assert c.geom.nariz.Ln == pytest.approx(0.065)  # conserva la longitud de la base


def test_deep_merge_no_muta():
    a = {"x": {"y": 1, "z": 2}}
    b = deep_merge(a, {"x": {"y": 5}})
    assert b == {"x": {"y": 5, "z": 2}} and a == {"x": {"y": 1, "z": 2}}


def test_barrido_zip_y_rutas_de_secciones(raw):
    raw["barrido"] = {"modo": "zip", "parametros": {
        "L_mm": [300, 350], "D_mm": [60, 65], "electronica.masa_g": [100, 200]}}
    raw["configuraciones"] = {}
    cfg = cargar(raw)
    assert [c.id for c in cfg.casos] == ["L300_D60_masa_g100", "L350_D65_masa_g200"]
    assert cfg.casos[1].electronica.me == pytest.approx(0.200)


def test_barrido_ruta_inexistente(raw):
    barrido(raw, 350, 65, aletas__h_punta__valor=[10])
    with pytest.raises(ConfigError) as e:
        cargar(raw)
    assert "aletas.h_punta.valor" in str(e.value)


def test_longitud_referencia_total(raw):
    raw["geometria_base"]["longitud_referencia"] = "total"
    cfg = cargar(solo(raw, 390, 65))
    g = cfg.casos[0].geom
    assert g.L == pytest.approx(0.390 - 0.03998)
    assert g.L_total == pytest.approx(0.390 + 0.02e-3)  # el offset bottom asoma 0.02 mm


def test_h_tip_por_radio_de_punta(raw):
    raw["geometria_base"]["aletas"]["h_tip"] = {"modo": "r_tip_absoluto_mm", "valor": 50.74}
    g = cargar(solo(raw, 350, 65)).casos[0].geom.aletas
    assert g.r_tip == pytest.approx(0.05074)
    assert g.h == pytest.approx(0.05074 - g.r_LE)


def test_popa_cerrada_agrega_tapa(raw):
    raw["geometria_base"]["cola"]["popa"] = "cerrada"
    c = cargar(solo(raw, 350, 65)).casos[0]
    tapa = [m for m in c.mamparos if m.nombre == "tapa_popa"]
    assert len(tapa) == 1
    assert tapa[0].x + tapa[0].e == pytest.approx(c.geom.L)


@pytest.mark.parametrize("mutar, fragmento", [
    (lambda r: r["geometria_base"]["nariz"]["longitud"].update(valor=5.0), "L_n + L_t"),
    (lambda r: r["geometria_base"]["cola"]["diametro_popa"].update(valor=1.2), "R_a"),
    (lambda r: r["pared"]["por_defecto"][0].update(material="unobtainium"), "unobtainium"),
    (lambda r: r["pared"]["por_defecto"][1].update(espesor_mm=0), "espesor_mm"),
    (lambda r: r["lastre"].update(factor_llenado=0), "factor_llenado"),
    (lambda r: r["lastre"].update(factor_llenado=1.2), "factor_llenado"),
    (lambda r: r["rellenos"][4]["granular"].update(empaquetamiento=1.0), "empaquetamiento"),
    (lambda r: r["estabilidad"].update(SM_max_cal=1.0), "SM_max_cal"),
    (lambda r: r["geometria_base"]["cola"].update(forma="potencia", parametro=None), "parametro"),
    (lambda r: r["geometria_base"]["cola"].update(parametro=1.5), "K de 'parabolica'"),
    (lambda r: r["geometria_base"]["aletas"]["c_r"].update(valor=70), "c_r"),
    (lambda r: r["geometria_base"]["aletas"]["h_tip"].update(valor=-1), "h_tip"),
    (lambda r: r["geometria_base"]["aletas"].update(r_interior_mm=25), "r_interior"),
    (lambda r: r["geometria_base"]["aletas"]["x_tip_le"].update(valor=200), "x_tip_le"),
    (lambda r: r["geometria_base"]["aletas"]["x_tip_le"].update(valor=80) or
               r["geometria_base"]["aletas"]["extension"].update(valor=0), "x_tip_le"),
    (lambda r: r["geometria_base"]["cola"].update(popa="cerrada", mamparo_popa={"espesor_mm": 0,
                                                                               "material": "PLA"}),
     "mamparo_popa"),
])
def test_validacion(raw, mutar, fragmento):
    solo(raw, 350, 65)
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
