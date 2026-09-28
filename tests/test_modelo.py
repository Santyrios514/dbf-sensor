"""Modelo completo (sección 8.2) y piezas de Fase B: Barrowman, ventana, trim."""

import copy
import math

import numpy as np
import pytest

from sensor_lastre import aletas as mod_aletas
from sensor_lastre.analisis import Tolerancias, analizar_caso
from sensor_lastre.barrowman import cuerpo_revolucion
from sensor_lastre.config import cargar
from sensor_lastre.estabilidad import ModeloLastre, alpha_trim
from sensor_lastre.geometria import malla
from sensor_lastre.perfiles import radio_exterior, radio_nariz

from conftest import PERDIGON, barrido, con_rellenos, solo

MM = 1e-3


def _aletas_efectivas(raw):
    """Aletas del barrido que sí permiten SM = 1.5 (h_tip = 30, extensión = 40)."""
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 30
    raw["geometria_base"]["aletas"]["extension"]["valor"] = 40
    return raw


def _caso_cilindro(raw, modo_e="detras_del_lastre"):
    """Nariz cónica corta y lastre que arranca en el cuerpo: ℓ* cae en el tramo cilíndrico."""
    raw["geometria_base"]["nariz"].update(forma="conica", longitud={"modo": "absoluto_mm", "valor": 20})
    raw["lastre"].update(x_inicio={"modo": "absoluto_mm", "valor": 30}, fraccion_max_L=None)
    raw["electronica"]["modo"] = modo_e
    raw["lastre"]["lastre_trasero"] = modo_e == "detras_del_lastre"  # el tapón trasero exige electrónica detrás
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
    cfg = cargar(solo(_aletas_efectivas(raw), 350, 65))
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
    raw["estabilidad"].update(SM_min_cal=1.5, SM_max_cal=1.7)
    cfg = cargar(solo(_aletas_efectivas(raw), 350, 65))
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
    cfg = cargar(barrido(raw, [270, 350, 400], [60, 65, 70]))
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
    nariz = cuerpo_revolucion(x, radio_exterior(x, g), 0.0, g.nariz.Ln, np.pi * g.R**2, "nariz")
    assert nariz.CNa == pytest.approx(2.0, rel=1e-9)
    assert nariz.x / g.nariz.Ln == pytest.approx(kn, rel=1e-3)


def test_ogiva_esbelta_kn():
    R, Ln = 10 * MM, 200 * MM  # fineza 10
    x = malla(Ln, 0.05 * MM)
    r = radio_nariz(x, "ogiva_tangente", None, R, Ln)
    V = np.trapezoid(np.pi * r**2, x)
    kn = 1 - V / (np.pi * R**2 * Ln)
    assert kn == pytest.approx(0.466, abs=2e-3)


def test_freeform_trapezoidal_igual_a_barrowman_clasico():
    """Una aleta trapezoidal sobre un cilindro, pasada por el método de franjas, reproduce las
    fórmulas cerradas de Barrowman (C_Nα y x_CP en el cuarto de la cuerda media)."""
    n, cr, ct, s, xs, R, x0 = 4, 60 * MM, 30 * MM, 50 * MM, 30 * MM, 30 * MM, 100 * MM
    P = np.array([[x0, R], [x0 + xs, R + s], [x0 + xs + ct, R + s], [x0 + cr, R]])
    p = mod_aletas.ParamsAleta(n=n, rotacion=0, c_r=cr, x_s=xs, h=s, r_tip_obj=None, e=0, r_in=0,
                               delta_b=0, eps=0, t=2 * MM, material="PLA", rho=1240, phi=1)
    g = mod_aletas.GeomAleta(params=p, x_LE=x0, r_LE=R, R_a=R, h=s, puntos=P - [x0, R], poligono=P,
                             origen=("punto",) * 4)
    c = mod_aletas.barrowman_freeform(g, np.pi * R**2, 0.0, 481)
    D = 2 * R
    lm = np.hypot(s, xs + ct / 2 - cr / 2)
    CNa = (1 + R / (s + R)) * (4 * n * (s / D) ** 2) / (1 + np.sqrt(1 + (2 * lm / (cr + ct)) ** 2))
    xf = x0 + xs * (cr + 2 * ct) / (3 * (cr + ct)) + (cr + ct - cr * ct / (cr + ct)) / 6
    assert c.CNa == pytest.approx(CNa, rel=1e-4)
    assert c.x_CP == pytest.approx(xf, abs=0.01 * MM)


def test_barrowman_interno_vs_openrocket_base(raw):
    """Aletas del .ork: OpenRocket 24.12 da C_Nα = 2.638 y x_CP = 342.41 mm (M = 0.3, medido con
    el jar); el método interno queda dentro de 3 % y 2 mm."""
    raw["condiciones_vuelo"]["mach"] = 0.3
    caso = cargar(raw).caso("base_ork")
    c = mod_aletas.barrowman_freeform(caso.geom.aletas, np.pi * caso.geom.R**2, 0.3, 48)
    assert c.CNa == pytest.approx(2.638, rel=0.03)
    assert c.x_CP / MM == pytest.approx(342.41, abs=2.0)
    assert c.span / MM == pytest.approx(31.24, abs=0.01)  # como FinSet.getSpan()


# --------------------------------------------------------------------------- trim


def test_trim_en_CG_es_cero(raw):
    raw["remolque"]["x_T"] = {"modo": "en_CG"}
    cfg = cargar(solo(_aletas_efectivas(raw), 270, 60))
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
    rr = analizar_caso(cargar(raw).caso("base_ork"), rellenos=[cargar(raw).rellenos[0]]).rellenos[0]
    f = rr.fila
    assert f["L_total_mm"] == pytest.approx(390.0, abs=0.01)
    assert f["cabe_largo"] is False  # 390 > 340 mm
    assert f["h_env_mm"] == pytest.approx(101.48, abs=0.01) and f["cabe_alto"] is True
    assert f["cabe_ancho"] is None  # ancho de la bahía pendiente


# --------------------------------------------------------------------------- SM_∞ y banderas v2


def test_SM_inf_limite_masa_infinita(raw):
    caso = cargar(solo(_aletas_efectivas(raw), 350, 65)).casos[0]
    rr = analizar_caso(caso, rellenos=[cargar(raw).rellenos[0]]).rellenos[0]
    mod, x_CP, D = rr.modelo, rr.fila["x_CP_mm"] * MM, caso.geom.D
    for ell in (20 * MM, rr.ventana.ell_star, 120 * MM):
        m0 = float(mod.m(ell) - mod.m_b(ell))
        rho = 1e4 * m0 / float(mod.V_b(ell))
        pesado = ModeloLastre(caso=caso, cav=mod.cav, mv=mod.mv, x_b0=mod.x_b0, rho_b=rho)
        SM = (x_CP - float(pesado.x_CG(ell))) / D
        assert abs(SM - float(pesado.SM_inf(ell, x_CP, D))) < 0.01


def test_SM_descomposicion_exacta(raw):
    """SM(ℓ) = SM_∞(ℓ) − m_0/(m_0 + m_b)·(x_0 − x̄_b)/D con masa finita (v2 §3.6)."""
    caso = cargar(solo(_aletas_efectivas(raw), 350, 65)).casos[0]
    mod = analizar_caso(caso, rellenos=[cargar(raw).rellenos[0]]).rellenos[0].modelo
    x_CP, D, ell = 0.2, caso.geom.D, 60 * MM
    mb = float(mod.m_b(ell))
    m0 = float(mod.m(ell)) - mb
    x0 = (float(mod.x_CG(ell)) * float(mod.m(ell)) - mb * float(mod.xbar_b(ell))) / m0
    SM = (x_CP - float(mod.x_CG(ell))) / D
    rhs = float(mod.SM_inf(ell, x_CP, D)) - m0 / (m0 + mb) * (x0 - float(mod.xbar_b(ell))) / D
    assert SM == pytest.approx(rhs, abs=1e-12)


def test_banderas_por_geometria_y_margen(raw):
    raw["estabilidad"]["SM_min_cal"] = 1.5
    raw["barrido"]["parametros"].pop("aletas.rotacion_deg")
    cfg = cargar(raw)
    plomo = cfg.rellenos[0]
    f = analizar_caso(cfg.caso("L350_D65_h_tip20_extension0"), rellenos=[plomo]).rellenos[0].fila
    assert "SM_inalcanzable_por_geometria" in f["banderas"]
    assert math.isnan(f["ell_SM_min_mm"]) and math.isnan(f["ell_SM_max_mm"])
    assert f["V_util_geo_cm3"] > 0  # las columnas geométricas sí se escriben
    f = analizar_caso(cfg.caso("base_ork"), rellenos=[plomo]).rellenos[0].fila
    assert "margen_SM_bajo" in f["banderas"] and "SM_inalcanzable_por_geometria" not in f["banderas"]
    assert f["SM_inf_cal"] == pytest.approx(1.78, abs=0.05)
    f = analizar_caso(cfg.caso("L350_D65_h_tip30_extension40"), rellenos=[plomo]).rellenos[0].fila
    assert not ({"margen_SM_bajo", "SM_inalcanzable"} & set(f["banderas"].split(";")))


def test_dif_masa_alta(raw):
    cfg = cargar(raw)
    caso = cfg.caso("base_ork")
    ok = {"nariz": 26.0144e-3, "cuerpo": 105.5151e-3, "cola": 27.4883e-3, "aletas": 25.317e-3}  # OR 24.12
    res = analizar_caso(caso, rellenos=[cfg.rellenos[0]], masas_or=ok, tol=Tolerancias.desde(cfg.verificacion))
    assert "dif_masa_alta" not in res.banderas
    assert abs(res.dif_masa_total) < 0.01
    malo = dict(ok, cola=2 * ok["cola"])  # p. ej. transición marcada Filled en el .ork
    res = analizar_caso(caso, rellenos=[cfg.rellenos[0]], masas_or=malo, tol=Tolerancias.desde(cfg.verificacion))
    assert "dif_masa_alta" in res.banderas


# --------------------------------------------------------------------------- criterio: máxima masa


def _caso_opt(raw, trasero, h=30, e=40, SM_min=1.0):
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = h
    raw["geometria_base"]["aletas"]["extension"]["valor"] = e
    raw["lastre"]["lastre_trasero"] = trasero
    raw["estabilidad"]["SM_min_cal"] = SM_min
    raw["lastre"]["masa_max_sensor_g"] = None  # óptimo sin tope de masa
    cfg = cargar(solo(raw, 350, 65))
    return cfg.casos[0], cfg.rellenos[0]


def test_tapon_trasero_llena_hasta_SM_min(raw):
    caso, plomo = _caso_opt(raw, True)
    f = analizar_caso(caso, rellenos=[plomo]).rellenos[0].fila
    assert f["limitante_masa"] == "SM"
    assert f["SM_cal"] == pytest.approx(1.0, abs=1e-4)
    assert f["ell_SM_max_mm"] == pytest.approx(f["ell_geo_mm"])  # el delantero se llenó primero
    assert 0 < f["ell_trasero_mm"] and f["m_relleno_trasero_g"] > 0


def test_tapon_trasero_nunca_reduce_masa(raw):
    c0, plomo = _caso_opt(copy.deepcopy(raw), False)
    c1, _ = _caso_opt(copy.deepcopy(raw), True)
    f0 = analizar_caso(c0, rellenos=[plomo]).rellenos[0].fila
    f1 = analizar_caso(c1, rellenos=[plomo]).rellenos[0].fila
    assert f0["limitante_masa"] == "geometria" and f0["SM_cal"] > 1.0
    assert f1["m_relleno_g"] > f0["m_relleno_g"]
    assert f1["m_relleno_g"] - f1["m_relleno_trasero_g"] == pytest.approx(f0["m_relleno_g"], rel=1e-9)


def test_tapon_trasero_sin_margen_de_SM(raw):
    """Aletas débiles: el SM limita al delantero antes de ℓ_geo y no se agrega lastre trasero."""
    caso, plomo = _caso_opt(raw, True, h=20, e=20)
    f = analizar_caso(caso, rellenos=[plomo]).rellenos[0].fila
    assert f["limitante_masa"] == "SM" and f["ell_trasero_mm"] == 0
    assert f["ell_SM_max_mm"] < f["ell_geo_mm"]


def test_tapon_trasero_llena_la_cola_si_sobra_SM(raw):
    caso, plomo = _caso_opt(raw, True, h=50, e=40)
    rr = analizar_caso(caso, rellenos=[plomo]).rellenos[0]
    f, mod = rr.fila, rr.modelo
    assert f["limitante_masa"] == "geometria" and f["SM_cal"] > 1.0
    assert f["ell_trasero_mm"] * MM == pytest.approx(mod.ell2_max(rr.ventana.ell_SM_max))


def test_optimo_con_SM_min_maximiza_masa(raw):
    """Barrer ℓ_2 a mano: ninguna ℓ_2 mayor que la elegida cumple SM ≥ 1."""
    caso, plomo = _caso_opt(raw, True)
    rr = analizar_caso(caso, rellenos=[plomo]).rellenos[0]
    mod, ell1 = rr.modelo, rr.ventana.ell_SM_max
    x_CP, D = rr.fila["x_CP_mm"] * MM, caso.geom.D
    ell2 = rr.fila["ell_trasero_mm"] * MM
    for l2 in np.linspace(ell2 + 0.1 * MM, mod.ell2_max(ell1), 20):
        assert (x_CP - float(mod.x_CG_con_trasero(ell1, l2))) / D < 1.0


def test_presupuesto_usa_tapon_trasero(raw):
    raw["presupuesto_masa"]["lastre_g"] = [1000, 6000, 20000]
    caso, plomo = _caso_opt(raw, True)
    pres = analizar_caso(caso, rellenos=[plomo]).rellenos[0].presupuestos
    p1, p6, p20 = pres
    assert p1["ell_trasero_mm"] == 0 and not p1["no_cabe"]
    assert p6["ell_trasero_mm"] > 0 and not p6["no_cabe"]  # 6 kg > 5.2 kg que caben adelante
    assert p6["m_total_g"] - p1["m_total_g"] == pytest.approx(5000, abs=0.5)  # la masa pedida se respeta
    assert p20["no_cabe"]


def test_ranking_filtra_y_ordena():
    import pandas as pd
    from sensor_lastre.optimizacion import ranking
    df = pd.DataFrame({
        "config_id": ["a", "b", "c", "d"], "relleno": ["pb"] * 4,
        "m_relleno_g": [5000, 7000, np.nan, 9000],
        "cabe_alto": [True, True, True, False], "cabe_ancho": [None, None, None, None],
    })
    r = ranking(df, ["cabe_alto", "cabe_ancho"])
    assert list(r["config_id"]) == ["b", "a"] and list(r["puesto"]) == [1, 2]


def test_ranking_desempata_por_SM_y_bahia():
    import pandas as pd
    from sensor_lastre.optimizacion import ranking
    df = pd.DataFrame({
        "config_id": ["a", "b", "c"], "relleno": ["pb"] * 3, "m_relleno_g": [6944.2, 6943.9, 6944.0],
        "SM_cal": [1.2, 1.5, 1.5], "D_acostado_mm": [100, 114, 100], "D_parado_mm": [141, 161, 141],
        "cabe_alto": [True] * 3, "cabe_ancho": [None] * 3,
    })
    assert list(ranking(df, ["cabe_alto"])["config_id"]) == ["c", "b", "a"]


def test_ranking_con_tope_ordena_por_masa_total():
    """Con el tope, la masa total empata y gana el SM aunque la otra lleve más plomo."""
    import pandas as pd
    from sensor_lastre.optimizacion import ranking
    df = pd.DataFrame({
        "config_id": ["liviana", "grande"], "relleno": ["pb"] * 2, "m_relleno_g": [5864.0, 5830.0],
        "m_total_g": [6200.0, 6200.0], "SM_cal": [1.26, 1.80], "D_acostado_mm": [114, 114],
        "D_parado_mm": [161, 161], "cabe_alto": [True] * 2, "cabe_ancho": [None] * 2,
    })
    assert list(ranking(df, [])["config_id"]) == ["grande", "liviana"]


# --------------------------------------------------------------------------- costo, D máx. y dibujo


def test_costo_relleno(raw):
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 40
    cfg = cargar(con_rellenos(solo(raw, 350, 65), PERDIGON))
    for rel in cfg.rellenos:
        f = analizar_caso(cfg.casos[0], rellenos=[rel]).rellenos[0].fila
        assert f["costo_usd_kg"] == pytest.approx(rel.costo_usd_kg)
        assert f["costo_relleno_usd"] == pytest.approx(f["m_relleno_g"] * 1e-3 * rel.costo_usd_kg, rel=1e-9)
    c = {r.nombre: r.costo_usd_kg for r in cfg.rellenos}
    assert c["plomo_macizo"] == raw["costos_usd_kg"]["plomo"]
    # perdigón + epoxy: promedio ponderado por fracción de masa
    w = 0.62 * 11340 / (0.62 * 11340 + 0.38 * 1150)
    assert c["perdigon_pb_epoxy"] == pytest.approx(w * 6 + (1 - w) * raw["costos_usd_kg"]["epoxy"])


def test_costo_sin_precio_es_nan(raw):
    raw["costos_usd_kg"].pop("acero")
    cfg = cargar(con_rellenos(solo(raw, 350, 65), {"nombre": "acero_macizo", "material": "acero"}))
    assert math.isnan({r.nombre: r.costo_usd_kg for r in cfg.rellenos}["acero_macizo"])


@pytest.mark.parametrize("rot", [0, 45])
def test_diametro_maximo_con_aletas(raw, rot):
    raw["geometria_base"]["aletas"]["rotacion_deg"] = rot
    caso = cargar(solo(raw, 350, 65)).casos[0]
    f = analizar_caso(caso, rellenos=[cargar(raw).rellenos[0]]).rellenos[0].fila
    r_tip = caso.geom.aletas.r_tip / MM
    assert f["D_parado_mm"] == pytest.approx(2 * r_tip)  # no depende de la rotación
    assert f["D_parado_mm"] >= max(f["h_env_mm"], f["w_env_mm"]) - 1e-9
    # guardado acostado con aletas diagonales: √2·r_tip, sea cual sea la rotación de vuelo
    assert f["D_acostado_mm"] == pytest.approx(math.sqrt(2) * r_tip)


def test_dibujo(raw, tmp_path):
    from sensor_lastre.figuras import fig_dibujo
    raw["geometria_base"]["aletas"].update(rotacion_deg=45)
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 40
    cfg = cargar(solo(raw, 350, 65))
    res = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]])
    ruta = tmp_path / "dibujo.png"
    fig_dibujo(res, res.rellenos[0], ruta)
    assert ruta.stat().st_size > 20_000


# --------------------------------------------------------------------------- amarre y trim


def _caso_recomendado(raw, **remolque):
    raw["geometria_base"]["aletas"].update(rotacion_deg=45)
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 40
    raw["remolque"].update(remolque)
    cfg = cargar(solo(raw, 350, 65))
    return cfg.casos[0], cfg.rellenos[0]


def test_amarre_en_CG_por_defecto(raw):
    caso, plomo = _caso_recomendado(raw)
    assert caso.remolque.modo == "en_CG"
    f = analizar_caso(caso, rellenos=[plomo]).rellenos[0].fila
    assert f["alpha_trim_deg"] == pytest.approx(0.0, abs=1e-9)
    assert "trim_no_lineal" not in f["banderas"]


def test_tolerancia_amarre(raw):
    from sensor_lastre import G0
    caso, plomo = _caso_recomendado(raw)
    f = analizar_caso(caso, rellenos=[plomo]).rellenos[0].fila
    m, xcg, xcp = f["m_total_g"] * 1e-3, f["x_CG_mm"] * MM, f["x_CP_mm"] * MM
    S, q = math.pi * caso.geom.D**2 / 4, caso.vuelo.q
    tol = math.radians(5) * q * S * f["CN_alpha_rad"] * (xcp - xcg) / (m * G0)
    assert f["tol_amarre_mm"] == pytest.approx(tol / MM, rel=1e-9)
    # consistencia: un error de x_T igual a la tolerancia da ≈ 5° de trim
    a = math.degrees(alpha_trim(m, xcg, xcg - tol, xcp, q, S, f["CN_alpha_rad"]))
    assert a == pytest.approx(5.0, rel=0.02)
    assert 0.5 < f["tol_amarre_mm"] < 3  # del orden del mm con ~7 kg a 30 m/s


def test_trim_no_lineal(raw):
    caso, plomo = _caso_recomendado(raw, x_T={"modo": "relativo_L", "valor": 0.15})
    f = analizar_caso(caso, rellenos=[plomo]).rellenos[0].fila
    assert abs(f["alpha_trim_deg"]) > 15 and "trim_no_lineal" in f["banderas"]


def test_validacion_alphas_remolque(raw):
    from sensor_lastre.config import ConfigError
    raw["remolque"].update(alpha_trim_max_deg=20, alpha_lineal_max_deg=15)
    with pytest.raises(ConfigError, match="alpha_trim_max_deg"):
        cargar(solo(raw, 350, 65))


def test_bahia_se_chequea_guardado_acostado(raw):
    """h_tip = 40 en vuelo + (alto 141.5 mm) cabe guardado en diagonal (100 mm) en 122 mm."""
    raw["geometria_base"]["aletas"].update(rotacion_deg=0)
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 40
    cfg = cargar(solo(raw, 350, 65))
    f = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0].fila
    assert f["h_env_mm"] == pytest.approx(141.48, abs=0.01)  # sección en vuelo (+)
    assert f["D_acostado_mm"] == pytest.approx(100.04, abs=0.01)
    assert f["cabe_alto"] is True
    raw["envolvente"]["rotacion_guardado_deg"] = 0  # si se guardara en +, no cabría
    cfg = cargar(solo(raw, 350, 65))
    f = analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0].fila
    assert f["cabe_alto"] is False


# --------------------------------------------------------------------------- tope de masa del sensor


def _fila_plomo(raw):
    cfg = cargar(solo(raw, 350, 65))
    return analizar_caso(cfg.casos[0], rellenos=[cfg.rellenos[0]]).rellenos[0].fila


def test_masa_max_sensor_limita_el_optimo(raw):
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 40
    raw["lastre"]["masa_max_sensor_g"] = None
    libre = _fila_plomo(raw)
    tope = libre["m_total_g"] - 1000
    raw["lastre"]["masa_max_sensor_g"] = tope
    f = _fila_plomo(raw)
    assert f["m_total_g"] == pytest.approx(tope, abs=0.5)
    assert f["m_relleno_g"] == pytest.approx(libre["m_relleno_g"] - 1000, abs=0.5)
    assert f["limitante_masa"] == "masa_max_sensor"
    assert f["SM_cal"] >= raw["estabilidad"]["SM_min_cal"] - 1e-6
    assert f["SM_cal"] > libre["SM_cal"]  # sale plomo del tapón trasero: el CG avanza


def test_masa_max_sensor_holgada_no_cambia_nada(raw):
    raw["lastre"]["masa_max_sensor_g"] = None
    libre = _fila_plomo(raw)
    raw["lastre"]["masa_max_sensor_g"] = 1e6
    f = _fila_plomo(raw)
    assert f["m_relleno_g"] == pytest.approx(libre["m_relleno_g"])
    assert f["limitante_masa"] == libre["limitante_masa"]


def test_masa_max_sensor_sin_SM(raw):
    raw["geometria_base"]["aletas"]["h_tip"]["valor"] = 20  # aleta chica: el SM depende del plomo
    raw["lastre"]["masa_max_sensor_g"] = 400  # casi sin plomo: el SM no llega a SM_min
    f = _fila_plomo(raw)
    assert "SM_inalcanzable_con_masa_max" in f["banderas"]
    assert not math.isfinite(f.get("m_relleno_g", math.nan))
