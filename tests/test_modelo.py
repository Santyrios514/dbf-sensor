"""Modelo completo (sección 8.2) y piezas de Fase B: Barrowman, ventana, trim."""

import math

import numpy as np
import pytest

from sensor_lastre.analisis import analizar_caso
from sensor_lastre.barrowman import aletas, cp_interno
from sensor_lastre.config import cargar
from sensor_lastre.estabilidad import alpha_trim
from sensor_lastre.geometria import malla
from sensor_lastre.perfiles import radio_exterior, radio_nariz

from conftest import solo

MM = 1e-3


def _caso_cilindro(raw, modo_e="detras_del_lastre"):
    """Nariz cónica corta y lastre que arranca en el cuerpo: ℓ* cae en el tramo cilíndrico."""
    raw["geometria_base"]["nariz"].update(forma="conica", longitud={"modo": "absoluto_mm", "valor": 20})
    raw["lastre"].update(x_inicio={"modo": "absoluto_mm", "valor": 30}, fraccion_max_L=None)
    raw["electronica"]["modo"] = modo_e
    if modo_e == "fija":
        raw["electronica"]["x_inicio_mm"] = 250
    cfg = cargar(solo(raw, 350, 65))
    return cfg.casos[0], cfg.rellenos[0]


@pytest.mark.parametrize("modo_e", ["detras_del_lastre", "fija"])
def test_condicion_del_minimo(raw, modo_e):
    caso, plomo = _caso_cilindro(raw, modo_e)
    rr = analizar_caso(caso, rellenos=[plomo]).rellenos[0]
    mod, v = rr.modelo, rr.ventana
    ls = v.ell_star
    assert 0 < ls < rr.fila["ell_geo_mm"] * MM
    xb = mod.x_b0 + ls
    assert caso.geom.nariz.Ln < xb < caso.geom.x_cola  # en el cilindro
    rhs = float(mod.x_CG(ls)) - caso.electronica.me * mod.dxe_dell() / (
        plomo.rho_b * caso.lastre.fll * float(mod.cav.area(xb)))
    assert abs(xb - rhs) < 0.5 * MM


def test_derivada_analitica(raw):
    caso, plomo = _caso_cilindro(raw)
    mod = analizar_caso(caso, rellenos=[plomo]).rellenos[0].modelo
    h = 0.05 * MM
    for ell in (10 * MM, 40 * MM, 90 * MM):
        num = (mod.x_CG(ell + h) - mod.x_CG(ell - h)) / (2 * h)
        assert float(mod.dxCG_dell(ell)) == pytest.approx(float(num), rel=2e-3, abs=1e-5)


def test_monotonia_e_inversion(raw):
    cfg = cargar(raw)
    for caso in cfg.casos:
        rr = analizar_caso(caso, rellenos=[cfg.rellenos[0]]).rellenos[0]
        mod = rr.modelo
        ells = np.linspace(0, mod.ell_cavidad, 2000)
        V = mod.V_b(ells)
        assert np.all(np.diff(V) > 0)
        for ell in np.linspace(0.05, 0.95, 7) * mod.ell_cavidad:
            ell_inv = mod.ell_de_masa(float(mod.m_b(ell)), caso.numerico.tol)
            assert abs(ell_inv - ell) < 0.01 * MM


def test_ventana_toca_SM_min(raw):
    cfg = cargar(solo(raw, 350, 65))
    caso = cfg.casos[0]
    rr = analizar_caso(caso, rellenos=[cfg.rellenos[0]]).rellenos[0]
    v, mod = rr.ventana, rr.modelo
    D = caso.geom.D
    x_CP = rr.fila["x_CP_mm"] * MM
    ell_geo = rr.fila["ell_geo_mm"] * MM
    for ell in (v.ell_SM_min, v.ell_SM_max):
        if 0 < ell < ell_geo:
            assert (x_CP - float(mod.x_CG(ell))) / D == pytest.approx(caso.estabilidad.SM_min, abs=1e-4)
    assert v.ell_SM_min <= v.ell_star <= v.ell_SM_max


def test_SM_max_acota_la_ventana(raw):
    raw["estabilidad"]["SM_max_cal"] = 1.7
    cfg = cargar(solo(raw, 350, 65))
    rr = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0]
    sm = rr.curva["SM_cal"].to_numpy()
    assert sm.max() > 1.7  # la restricción está activa
    for a, b in rr.ventana.intervalos:
        for ell in np.linspace(a, b, 20):
            s = (rr.fila["x_CP_mm"] * MM - float(rr.modelo.x_CG(ell))) / cfg.casos[0].geom.D
            assert 1.5 - 1e-4 <= s <= 1.7 + 1e-4
    assert "no_unimodal" in rr.fila["banderas"]  # la ventana queda partida en dos


def test_SM_inalcanzable(raw):
    raw["estabilidad"]["SM_min_cal"] = 6.0
    cfg = cargar(solo(raw, 270, 60))
    rr = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert "SM_inalcanzable" in rr.fila["banderas"]
    assert math.isnan(rr.fila.get("V_util_cm3", math.nan))


def test_inviable_geo(raw):
    raw["lastre"]["zonas_prohibidas"] = [{"desde_mm": 0, "hasta_mm": 300}]
    cfg = cargar(solo(raw, 270, 60))
    rr = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert "inviable_geo" in rr.fila["banderas"]


@pytest.mark.parametrize("hol", [0.0, 2.0])
@pytest.mark.parametrize("cid, ref", [("L270_D60", 221.7), ("L350_D65", 477.3), ("L400_D70", 707.0)])
def test_regresion_chat(raw, cid, ref, hol):
    raw["pared"]["por_defecto"] = [{"material": "PLA", "espesor_mm": 2.0}]
    raw["lastre"].update(x_inicio={"modo": "absoluto_mm", "valor": 0}, fraccion_max_L=None,
                         margen_cola_mm=0)
    raw["electronica"].update(longitud_mm=100, masa_g=150, modo="detras_del_lastre", holgura_mm=hol)
    cfg = cargar(raw)
    fila = analizar_caso(cfg.caso(cid), rellenos=[cfg.rellenos[0]]).rellenos[0].fila
    assert fila["V_util_geo_cm3"] == pytest.approx(ref, rel=0.05)


# --------------------------------------------------------------------------- Barrowman


@pytest.mark.parametrize("forma, p, kn", [
    ("conica", None, 2 / 3), ("elipsoide", None, 1 / 3), ("potencia", 0.5, 0.5), ("potencia", 0.75, 0.6),
])
def test_kn_nariz(raw, forma, p, kn):
    raw["geometria_base"]["nariz"].update(forma=forma, parametro=p)
    g = cargar(solo(raw, 350, 65)).casos[0].geom
    x = malla(g.L, 0.1 * MM)
    nariz = cp_interno(g, x, radio_exterior(x, g)).partes[0]
    assert nariz.CNa == pytest.approx(2.0, rel=1e-9)
    assert nariz.x / g.nariz.Ln == pytest.approx(kn, rel=1e-3)


def test_ogiva_esbelta_kn():
    R, Ln = 10 * MM, 200 * MM  # fineza 10
    x = malla(Ln, 0.05 * MM)
    r = radio_nariz(x, "ogiva_tangente", None, R, Ln)
    V = np.trapezoid(np.pi * r**2, x)
    kn = 1 - V / (np.pi * R**2 * Ln)
    assert kn == pytest.approx(0.466, abs=2e-3)


def test_cola_conica_barrowman(raw):
    g = cargar(solo(raw, 350, 65)).casos[0].geom
    x = malla(g.L, 0.1 * MM)
    cola = cp_interno(g, x, radio_exterior(x, g)).partes[1]
    k = (2 * g.cola.Ra) / g.D
    assert cola.CNa == pytest.approx(2 * (k**2 - 1), rel=1e-9)
    dF_dR = 1 / k  # Barrowman usa delantero/trasero
    x_ref = g.x_cola + g.cola.Lt / 3 * (1 + (1 - dF_dR) / (1 - dF_dR**2))
    assert cola.x == pytest.approx(x_ref, abs=0.01 * MM)


def test_aletas_barrowman_referencia():
    """Caso de mano: 4 aletas, cr = 60, ct = 30, s = 50, x_s = 30, D = 60 (mm)."""
    from dataclasses import replace
    from sensor_lastre.config import Aletas, Cola, Geometria, Nariz
    a = Aletas(n=4, cr=60 * MM, ct=30 * MM, s=50 * MM, xs=30 * MM, tf=2 * MM, material="PLA",
               rho=1240, phi=1, x_r0=100 * MM)
    g = Geometria(L=200 * MM, D=60 * MM, nariz=Nariz("conica", None, 60 * MM),
                  cola=Cola("conica", 0.0, 30 * MM, True), aletas=a)
    c = aletas(g)
    lm = math.hypot(50, 30 + 15 - 30)
    CNa = (1 + 30 / 80) * (4 * 4 * (50 / 60) ** 2) / (1 + math.sqrt(1 + (2 * lm / 90) ** 2))
    xf = 100 + 30 * (60 + 60) / (3 * 90) + (90 - 1800 / 90) / 6
    assert c.CNa == pytest.approx(CNa, rel=1e-12)
    assert c.x / MM == pytest.approx(xf, rel=1e-12)


# --------------------------------------------------------------------------- trim


def test_trim_en_CG_es_cero(raw):
    raw["remolque"]["x_T"] = {"modo": "en_CG"}
    cfg = cargar(solo(raw, 270, 60))
    rr = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert np.allclose(rr.curva["alpha_trim_deg"], 0.0)
    assert rr.fila["x_T_mm"] == pytest.approx(rr.fila["x_CG_mm"])
    assert rr.fila["x_T_trim_cero_mm"] == pytest.approx(rr.fila["x_CG_mm"])


def test_trim_signo_y_remolque_detras_del_CP(raw):
    assert alpha_trim(1.0, 0.10, 0.05, 0.20, 500.0, 3e-3, 8.0) > 0
    assert alpha_trim(1.0, 0.04, 0.05, 0.20, 500.0, 3e-3, 8.0) < 0
    assert np.isnan(alpha_trim(1.0, 0.10, 0.25, 0.20, 500.0, 3e-3, 8.0))
    raw["remolque"]["x_T"] = {"modo": "relativo_L", "valor": 0.9}
    cfg = cargar(solo(raw, 270, 60))
    rr = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert "inestable_respecto_remolque" in rr.fila["banderas"]


def test_envolvente(raw):
    cfg = cargar(raw)
    rr = analizar_caso(cfg.caso("L270_D60"), rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert rr.fila["cabe_largo"] is True
    assert rr.fila["cabe_diametro"] is False  # 60 + 2·48 = 156 > 122
    raw["geometria_base"]["aletas"]["semienvergadura"] = {"modo": "absoluto_mm", "valor": 30}
    cfg = cargar(raw)
    rr = analizar_caso(cfg.caso("L350_D60"), rellenos=[cfg.rellenos[0]]).rellenos[0]
    assert rr.fila["cabe_largo"] is False and rr.fila["cabe_diametro"] is True
