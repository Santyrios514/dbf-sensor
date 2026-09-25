"""Fase A2 (v2 §8, puntos 1–4, 6 y 8): .ork sin JVM, perfiles de transición, aleta freeform,
fórmula del polígono, CP de cuerpos y envolvente."""

import math

import numpy as np
import pytest

from sensor_lastre import ork_xml, verificacion
from sensor_lastre.aletas import area_centroide, envolvente, es_simple
from sensor_lastre.barrowman import cuerpo_revolucion
from sensor_lastre.config import cargar
from sensor_lastre.geometria import malla
from sensor_lastre.perfiles import forma_or, radio_cola, radio_exterior

from conftest import CONFIG, ORK, barrido, solo

MM = 1e-3


@pytest.fixture(scope="module")
def cfg():
    return cargar(CONFIG)


@pytest.fixture(scope="module")
def ork():
    return ork_xml.leer(ORK)


# --------------------------------------------------------------------------- 1. lectura del .ork


def test_ork_metadatos(ork):
    assert ork.creador == "OpenRocket 24.12" and ork.version_formato == "1.10"
    assert ork.referencia == "maximum"
    assert [c.tipo for c in ork.componentes] == ["NoseCone", "BodyTube", "Transition", "FreeformFinSet"]
    assert ork.por_nombre("Aletas de forma libre").padre == "Transición"


def test_ork_valores_seccion_1_1(ork, cfg):
    v = verificacion.valores_desde_xml(ork, cfg.openrocket["nombres_componentes"])
    tol = 1e-6
    assert v.Ln == pytest.approx(0.065, abs=tol) and v.R_nariz == pytest.approx(0.0325, abs=tol)
    assert v.forma_nariz == "elipsoide"
    assert v.Lc == pytest.approx(0.220, abs=tol) and v.R_cuerpo == pytest.approx(0.0325, abs=tol)
    assert v.Lt == pytest.approx(0.065, abs=tol) and v.Ra == pytest.approx(0.0195, abs=tol)
    assert v.forma_cola == "parabolica" and v.param_cola == 1.0 and v.popa_abierta
    for t in (v.t_nariz, v.t_cuerpo, v.t_cola, v.fin_t):
        assert t == pytest.approx(0.0018, abs=tol)
    assert v.fin_n == 4 and v.fin_rot == 0.0
    assert v.fin_offset_bottom == pytest.approx(0.02e-3, abs=tol)
    assert all(r == 1342.0 for r in v.rho.values())


def test_ork_puntos_seccion_1_2(ork):
    pts = np.array(ork.por_nombre("Aletas de forma libre").puntos) / MM
    ref = np.array([[0, 0], [40, 20], [81.08, 20], [81.08, -30.11], [41.18, -30.10], [41.10, -11.24]])
    assert np.allclose(pts, ref, atol=0.005)


def test_config_base_coincide_con_ork(ork, cfg):
    v = verificacion.valores_desde_xml(ork, cfg.openrocket["nombres_componentes"])
    assert verificacion.comparar(v, cfg.caso("base_ork"), tol=cfg.verificacion["tol_geometria_ork_mm"] * MM) == []


def test_comparar_detecta_diferencias(ork, raw):
    raw["geometria_base"]["cola"]["forma"] = "conica"
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 25
    cfg = cargar(solo(raw, 350, 65))
    v = verificacion.valores_desde_xml(ork, cfg.openrocket["nombres_componentes"])
    dif = verificacion.comparar(v, cfg.casos[0], tol=1e-6)
    assert any("cola.forma" in d for d in dif) and any("aletas.P1" in d for d in dif)


# --------------------------------------------------------------------------- 2. perfiles de transición


def test_transicion_parabolica_seccion_1_3(cfg):
    g = cfg.caso("base_ork").geom
    x_LE = 350.02e-3 - 41.10e-3
    assert radio_exterior(np.array([x_LE]), g)[0] / MM == pytest.approx(30.74, abs=0.01)
    assert radio_exterior(np.array([g.L]), g)[0] / MM == pytest.approx(19.50, abs=1e-9)
    u = np.linspace(0, 1, 11)
    ana = 19.5 + 13 * (1 - u**2)
    assert np.allclose(radio_cola(u * 0.065, "parabolica", 0.0325, 0.0195, 0.065, 1.0) / MM, ana, atol=1e-9)


@pytest.mark.parametrize("forma, p", [("conica", None), ("ogiva", 1.0), ("ogiva", 0.5), ("elipsoide", None),
                                      ("potencia", 0.5), ("parabolica", 0.0), ("parabolica", 0.75),
                                      ("haack", 0.0), ("haack", 1 / 3)])
def test_formas_cierran(forma, p):
    """f(0) = 0 y f(1) = radio: la cola empalma con el cuerpo y termina en R_a."""
    r = forma_or(np.array([0.0, 0.065]), forma, 0.013, 0.065, p)
    assert r[0] == pytest.approx(0.0, abs=1e-12) and r[1] == pytest.approx(0.013, rel=1e-9)


def test_parabolica_formula_general():
    xi = np.linspace(0, 1, 7)
    for K in (0.0, 0.5, 1.0):
        assert np.allclose(forma_or(xi, "parabolica", 1.0, 1.0, K), (2 * xi - K * xi**2) / (2 - K))


def test_ogiva_parametro_1_es_tangente():
    R, L = 0.03, 0.06
    x = np.linspace(0, L, 50)
    rho_o = (R**2 + L**2) / (2 * R)
    tang = np.sqrt(rho_o**2 - (L - x) ** 2) + R - rho_o
    assert np.allclose(forma_or(x, "ogiva", R, L, 1.0), tang, atol=1e-12)


# --------------------------------------------------------------------------- 3–4. aleta y polígono


def test_reconstruccion_aleta_reproduce_xml(ork, cfg):
    g = cfg.caso("base_ork").geom.aletas
    xml = np.array(ork.por_nombre("Aletas de forma libre").puntos)
    assert np.max(np.abs(g.puntos - xml)) < 0.1 * MM
    assert g.x_LE / MM == pytest.approx(308.92, abs=1e-6)
    assert g.r_LE / MM == pytest.approx(30.74, abs=0.01)
    assert g.r_tip / MM == pytest.approx(50.74, abs=0.01)
    assert g.x_TE / MM == pytest.approx(390.0, abs=0.01)
    assert g.puntos[4, 1] + g.r_LE == pytest.approx(0.63 * MM, abs=1e-9)  # borde interior en r_in
    assert es_simple(g.poligono)
    assert set(g.origen) == {"punto", "superficie"}


def test_aleta_sin_extension(raw):
    raw["geometria_base"]["aletas"]["extension"]["valor"] = 0
    g = cargar(solo(raw, 350, 65)).casos[0].geom.aletas
    assert len(g.puntos) == 4 and es_simple(g.poligono)
    assert g.x_TE == pytest.approx(g.x_LE + g.params.c_r)


@pytest.mark.parametrize("P, A, xc", [
    (np.array([[0, 0], [3, 0], [3, 2], [0, 2]], float), 6.0, 1.5),
    (np.array([[0, 0], [4, 0], [0, 3]], float), 6.0, 4 / 3),
    (np.array([[0, 0], [6, 0], [4, 2], [1, 2]], float), 9.0, (6**2 + 6 * 3 + 3**2 + 1 * (6 + 2 * 3)) / (3 * 9)),
])
def test_shoelace_analitico(P, A, xc):
    a, x = area_centroide(P)
    assert a == pytest.approx(A, abs=1e-9) and x == pytest.approx(xc, abs=1e-9)
    a2, x2 = area_centroide(P[::-1])  # la orientación no importa
    assert a2 == pytest.approx(A, abs=1e-9) and x2 == pytest.approx(xc, abs=1e-9)


def test_aleta_base_area_y_masa(cfg):
    g = cfg.caso("base_ork").geom.aletas
    assert g.area / MM**2 == pytest.approx(2.62e3, rel=0.03)
    assert g.masa * 1e3 == pytest.approx(25, rel=0.03)
    assert g.masa == pytest.approx(4 * 1342 * 1.8e-3 * g.area, rel=1e-12)


def test_solape_eje_opcional(raw):
    raw["geometria_base"]["aletas"].update(r_interior_mm=0.0, descontar_solape_eje=True)
    g = cargar(solo(raw, 350, 65)).casos[0].geom.aletas
    Lx = g.puntos[3, 0] - g.puntos[4, 0]
    assert g.volumen_solape() == pytest.approx(g.params.t**2 * Lx, rel=1e-12)
    # La spec dice "< 0.1 g"; con r_in = 0 son t²·L·ρ ≈ 0.17 g. Menor al 1 % de las aletas.
    assert g.volumen_solape() * g.params.rho / g.masa < 0.01


# --------------------------------------------------------------------------- 6. CP de cuerpos


@pytest.mark.parametrize("k", [0.0, 0.3, 0.6, 0.9])
def test_cp_tronco_de_cono(k):
    R, Lt = 0.0325, 0.065
    x = malla(Lt, 0.01 * MM)
    r = R + (k * R - R) * x / Lt
    c = cuerpo_revolucion(x, r, 0.0, Lt, np.pi * R**2, "cola")
    assert c.x == pytest.approx(Lt / 3 * (1 + 2 * k) / (1 + k), abs=1e-6)
    assert c.CNa == pytest.approx(2 * (k**2 - 1), rel=1e-9)
    if k == 0:
        assert c.x == pytest.approx(Lt / 3, abs=1e-6)  # cola cerrada: el error de la v1


# --------------------------------------------------------------------------- 8. envolvente


def test_envolvente_rotacion(cfg):
    g = cfg.caso("base_ork").geom
    r_tip = g.aletas.r_tip
    h, w = envolvente(r_tip, g.R, 4, 0.0)
    assert h / MM == pytest.approx(101.5, abs=0.05) and w == pytest.approx(h)
    h45, w45 = envolvente(r_tip, g.R, 4, math.radians(45))
    assert h45 == pytest.approx(math.sqrt(2) * r_tip) and w45 == pytest.approx(h45)


def test_envolvente_n_impar_y_cuerpo():
    h, w = envolvente(0.05, 0.03, 3, 0.0)
    assert h == pytest.approx(0.05 + 0.03)  # una aleta arriba, abajo manda el cuerpo (0.5 r_tip < R)
    h, _ = envolvente(0.02, 0.03, 4, 0.0)
    assert h == pytest.approx(0.06)  # aletas más cortas que el radio: manda el cuerpo


def test_barrido_de_rotacion(raw):
    barrido(raw, 350, 65, aletas__rotacion_deg=[0, 45])
    cfg = cargar(raw)
    assert [round(c.params_aleta.rotacion, 6) for c in cfg.casos] == [0.0, round(math.pi / 4, 6)]
