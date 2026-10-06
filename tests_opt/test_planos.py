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
