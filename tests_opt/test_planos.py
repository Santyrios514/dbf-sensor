"""Planos acotados (scripts/07_planos.py): cotas iguales al ranking y archivos reproducibles."""

import hashlib

import pytest

from sensor_opt import barrido, planos
from sensor_opt.config import cargar

from .conftest import malla_chica, raw_opt


@pytest.fixture(scope="module")
def cfg_y_fila():
    cfg = cargar(malla_chica(raw_opt()))
    rk = barrido.optimizar(cfg, refinar=False, procesos=1)
    fac = rk[rk["factible"]]
    assert len(fac) > 0
    return cfg, fac.iloc[0]


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_cotas_coinciden_con_el_ranking(cfg_y_fila, tmp_path):
    cfg, f = cfg_y_fila
    pl = planos.plano(cfg, f, tmp_path / "p", puesto=1, dpi=60)
    assert all(r.exists() for r in pl.rutas) and not pl.infactible
    for cota, col in planos.COL_RANKING.items():
        assert pl.cotas[cota][0] == pytest.approx(float(f[col]), abs=0.05), cota
    assert pl.cotas["L"][0] == pytest.approx(cfg.L * 1e3, abs=0.05)
    # el lastre dibujado es el del llenado con el mismo CP (sustituto sin validación)
    assert pl.cotas["plomo_delantero"][0] == pytest.approx(float(f["ell_mm"]), abs=0.05)
    assert pl.cotas["plomo_trasero"][0] == pytest.approx(float(f["ell_trasero_mm"]), abs=0.05)


def test_planos_reproducibles_byte_a_byte(cfg_y_fila, tmp_path):
    cfg, f = cfg_y_fila
    a = planos.plano(cfg, f, tmp_path / "a", dpi=60)
    b = planos.plano(cfg, f, tmp_path / "b", dpi=60)
    assert [_sha(p) for p in a.rutas] == [_sha(p) for p in b.rutas]


def test_altura_aparente_4_aletas_formula_cerrada():
    import math
    H, W, phi = planos.altura_aparente(35.0, 77.0, 2.0, 4)
    assert H == pytest.approx(math.sqrt(2) * (77.0 + 1.0), abs=1e-6)
    assert W == pytest.approx(H, abs=1e-6) and math.degrees(phi) == pytest.approx(45.0, abs=1e-4)


def test_altura_aparente_es_el_minimo_sobre_la_rotacion():
    import numpy as np
    for n, R, r_tip, t in ((3, 35.0, 77.0, 2.0), (4, 35.0, 40.0, 3.0), (6, 30.0, 60.0, 2.0)):
        H, _, _ = planos.altura_aparente(R, r_tip, t, n)
        bruto = min(planos._alto_ancho(R, r_tip, t, n, p)[0] for p in np.linspace(0, 2 * np.pi, 20001))
        assert H <= bruto + 1e-9 and H == pytest.approx(bruto, abs=1e-3)
        assert H >= 2 * R - 1e-12
