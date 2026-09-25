"""Pruebas analíticas de la sección 8.1 (tolerancia relativa 1e-4 salvo indicación)."""

import numpy as np
import pytest

from sensor_lastre.config import cargar
from sensor_lastre.geometria import construir_cavidad, erosionar, integrar, malla
from sensor_lastre.masas import masa_vacia
from sensor_lastre.perfiles import radio_cola, radio_exterior, radio_nariz

from conftest import solo

MM = 1e-3
DX = 0.1 * MM
REL = 1e-4


def vol_y_centroide(x, r):
    V = integrar(np.pi * r**2, x)
    return V, integrar(np.pi * r**2 * x, x) / V


def test_cilindro():
    x = malla(100 * MM, DX)
    r = np.full_like(x, 28 * MM)
    V, _ = vol_y_centroide(x, r)
    assert V / MM**3 == pytest.approx(246300.9, rel=REL)


@pytest.mark.parametrize("forma, k_vol, k_x", [("elipsoide", 2 / 3, 5 / 8), ("conica", 1 / 3, 3 / 4)])
def test_narices(forma, k_vol, k_x):
    R, Ln = 30 * MM, 60 * MM
    x = malla(Ln, DX)
    r = radio_nariz(x, forma, None, R, Ln)
    V, xc = vol_y_centroide(x, r)
    assert V == pytest.approx(k_vol * np.pi * R**2 * Ln, rel=REL)
    assert xc == pytest.approx(k_x * Ln, rel=REL)


def test_ogiva_tangente_y_haack_cierran_en_R():
    R, Ln = 30 * MM, 60 * MM
    for forma, p in [("ogiva_tangente", None), ("haack", 0.0), ("haack", 1 / 3), ("potencia", 0.5)]:
        r = radio_nariz(np.array([0.0, Ln]), forma, p, R, Ln)
        assert r[0] == pytest.approx(0.0, abs=1e-12)
        assert r[1] == pytest.approx(R, rel=1e-12)


def test_tronco_de_cono():
    R, Ra, Lt = 30 * MM, 18 * MM, 60 * MM
    u = malla(Lt, DX)
    r = radio_cola(u, "conica", R, Ra, Lt)
    V, _ = vol_y_centroide(u, r)
    assert V == pytest.approx(np.pi * Lt / 3 * (R**2 + R * Ra + Ra**2), rel=REL)


@pytest.mark.parametrize("forma", ["ogiva", "parabolica", "elipsoide", "haack"])
def test_colas_extensibles(forma):
    R, Ra, Lt = 30 * MM, 18 * MM, 60 * MM
    u = np.array([0.0, 1e-6, Lt])
    r = radio_cola(u, forma, R, Ra, Lt)
    assert r[0] == pytest.approx(R) and r[-1] == pytest.approx(Ra)
    assert (R - r[1]) / 1e-6 == pytest.approx(0.0, abs=1e-3)  # tangente al cuerpo


def test_erosion_cilindro_popa_cerrada():
    R, L, t = 28 * MM, 100 * MM, 2 * MM
    x = malla(L, DX)
    r = erosionar(np.full_like(x, R), t, DX, popa_cerrada=True)
    interior = (x >= t + DX) & (x <= L - t - DX)
    assert np.allclose(r[interior], R - t, atol=1e-12)
    ocup = x[r > 0]
    assert L - ocup[-1] == pytest.approx(t, abs=DX)  # pared de fondo
    assert ocup[0] == pytest.approx(t, abs=DX)  # pared de proa (el cilindro también está cerrado)


def test_erosion_popa_abierta_sin_fondo():
    """Popa abierta (v2 §3.5): el casco termina en un anillo de espesor t, sin pared de fondo."""
    R, L, t = 28 * MM, 100 * MM, 2 * MM
    x = malla(L, DX)
    r = erosionar(np.full_like(x, R), t, DX, popa_cerrada=False)
    assert r[-1] == pytest.approx(R - t)


def test_composicion_de_capas(raw):
    cfg = cargar(solo(raw, 350, 65))
    g = cfg.casos[0].geom
    x = malla(g.L, DX)
    r = radio_exterior(x, g)
    t1, t2 = 0.3 * MM, 1.5 * MM
    dos_pasos = erosionar(erosionar(r, t1, DX), t2, DX)
    un_paso = erosionar(r, t1 + t2, DX)
    assert np.max(np.abs(dos_pasos - un_paso)) <= DX


def test_erosion_exacta_vs_aproximada_en_el_cuerpo(raw):
    """Lejos de las esquinas la erosión exacta y la aproximada coinciden (r' = 0)."""
    from sensor_lastre.geometria import erosion_aproximada
    cfg = cargar(solo(raw, 350, 65))
    g = cfg.casos[0].geom
    x = malla(g.L, DX)
    r = radio_exterior(x, g)
    ex, ap = erosionar(r, 2 * MM, DX), erosion_aproximada(x, r, 2 * MM)
    cuerpo = (x > g.nariz.Ln + 5 * MM) & (x < g.x_cola - 5 * MM)
    assert np.allclose(ex[cuerpo], ap[cuerpo], atol=1e-9)


def test_conservacion_de_masa_una_capa(raw):
    raw["pared"]["por_defecto"] = [{"material": "PLA", "espesor_mm": 2.0, "fraccion_solida": 1.0}]
    cfg = cargar(solo(raw, 270, 60))
    caso = cfg.casos[0]
    cav = construir_cavidad(caso)
    mv = masa_vacia(caso, cav)
    rho = raw["materiales"]["PLA"]
    assert mv.m_casco == pytest.approx(rho * (cav.V_ext - cav.V_int_bruto), rel=1e-12)


def test_conservacion_de_masa_multicapa(raw):
    """Con varias capas y φ = 1, Σ m_k/ρ_k = V_ext − V_int."""
    cfg = cargar(solo(raw, 270, 60))
    caso = cfg.casos[0]
    cav = construir_cavidad(caso)
    mv = masa_vacia(caso, cav)
    rho = raw["materiales"]
    vol_capas = sum(c.m / rho[c.material] for c in mv.capas)
    assert vol_capas == pytest.approx(cav.V_ext - cav.V_int_bruto, rel=1e-12)


def test_paredes_distintas_por_estacion(raw):
    raw["pared"]["nariz"] = [{"material": "PLA", "espesor_mm": 4.0}]
    cfg = cargar(solo(raw, 350, 65))
    caso = cfg.casos[0]
    cav = construir_cavidad(caso)
    g = caso.geom
    i_c = np.searchsorted(cav.x, (g.nariz.Ln + g.x_cola) / 2)
    assert cav.r_i[i_c] == pytest.approx(g.R - 1.8 * MM, abs=1e-9)
    # la nariz, más gruesa, queda por dentro de lo que daría la pared por defecto
    i_n = np.searchsorted(cav.x, g.nariz.Ln / 2)
    assert cav.r_i[i_n] < erosionar(cav.r_e, 1.8 * MM, cav.dx)[i_n]
    # la cavidad es continua y no crece al pasar la unión
    assert np.all(np.isfinite(cav.r_i)) and cav.r_i.min() >= 0
