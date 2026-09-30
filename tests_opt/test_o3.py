"""Fase O3: barrido, ranking, Pareto, refinamiento, figuras y T9 (reproducibilidad)."""

import math

import numpy as np
import pandas as pd
import pytest

from sensor_opt import barrido, exportar
from sensor_opt.config import cargar

from .conftest import malla_chica


@pytest.fixture(scope="module")
def cfg_chica():
    from .conftest import raw_opt
    return cargar(malla_chica(raw_opt()))


@pytest.fixture(scope="module")
def ranking(cfg_chica):
    return barrido.optimizar(cfg_chica, procesos=2)


def test_T9_reproducible_y_sin_depender_de_procesos(cfg_chica, ranking):
    otra = barrido.optimizar(cfg_chica, procesos=1)
    pd.testing.assert_frame_equal(ranking, otra)


def test_ranking_ordenado_por_J_entre_factibles(ranking):
    fac = ranking[ranking["factible"]]
    assert len(fac) > 0 and list(fac.index) == list(range(len(fac)))  # factibles primero
    assert np.all(np.diff(fac["J"].to_numpy()) >= 0)
    assert list(fac["puesto"]) == list(range(1, len(fac) + 1))
    inf = ranking[~ranking["factible"]]
    assert inf["puesto"].isna().all() and (inf["motivos"] != "").all()


def test_factibles_cumplen_restricciones(cfg_chica, ranking):
    fac = ranking[ranking["factible"]]
    r = cfg_chica.restricciones
    assert (fac["SM_cal"] >= r.SM_min - 0.005).all()  # T5: ±0.005 cal (brentq con tol_raiz_mm) and (fac["SM_cal"] <= r.SM_max + 1e-9).all()
    assert (fac["m_total_g"] <= cfg_chica.m_max * 1e3 + 0.1).all()
    assert (fac["k"] >= r.k_min).all() and (fac["h_mm"] >= r.h_min * 1e3 - 1e-9).all()
    assert fac["restriccion_activa"].isin(["masa", "SM_min", "volumen", "tol_amarre"]).all()
    tol_min = cfg_chica.base_raw["remolque"].get("tol_amarre_min_mm")
    if tol_min is not None:
        assert (fac["tol_amarre_mm"] >= tol_min * (1 - 1e-3)).all()
    assert (fac["fineza_cola"] > 1.0).all()  # sin colas de base roma
    assert np.allclose(fac["D_ap_mm"], 2 * np.maximum(fac["D_mm"] / 2, fac["r_tip_mm"]))
    assert np.allclose(fac["h_45_mm"], np.maximum(fac["D_mm"], math.sqrt(2) * fac["r_tip_mm"]))


def test_pareto_entre_los_que_alcanzan_m_max(cfg_chica, ranking):
    p = ranking[ranking["en_pareto"]]
    assert len(p) > 0 and p["factible"].all()
    assert (p["m_total_g"] >= cfg_chica.m_max * 1e3 - 0.5).all()
    F = p[["f2", "f3", "f4"]].to_numpy()
    for i in range(len(F)):  # nadie del frente domina a otro del frente
        dom = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        assert not dom.any()


def test_refinamiento_a_medio_paso(cfg_chica, ranking):
    ref = ranking[ranking["refinamiento"]]
    assert len(ref) > 0 and not ranking["cand_id"].duplicated().any()
    assert set(ref["k"]) <= {0.3, 0.4, 0.5}  # k de la malla ± medio paso, dentro de [0.3, 0.5]
    assert ref["r_tip_rel_R"].between(1.6, 2.4).all() and ref["mu_cr"].between(0.75, 1.0).all()


def test_grupos_desde_ranking_reproduce_filas(cfg_chica, ranking):
    muestra = ranking.sample(20, random_state=1)
    gr = barrido.grupos_desde_ranking(muestra)
    otra = barrido.evaluar(cfg_chica, gr, procesos=1).set_index("cand_id")
    orig = muestra.set_index("cand_id")
    for col in ("J", "SM_cal", "m_total_g", "x_CP_sust_mm"):
        assert np.allclose(otra.loc[orig.index, col], orig[col], equal_nan=True)


def test_limite_de_evaluaciones(tmp_path):
    import subprocess
    import sys
    import yaml
    from .conftest import RAIZ, raw_opt
    raw = raw_opt()
    raw["base_config"] = str(RAIZ / "config" / "config.yaml")
    raw["ejecucion"]["max_evaluaciones"] = 10
    ruta = tmp_path / "opt.yaml"
    ruta.write_text(yaml.safe_dump(raw, allow_unicode=True))
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "04_optimizar_malla.py"), "--config", str(ruta)],
                       capture_output=True, text=True)
    assert r.returncode == 3 and "--forzar" in r.stderr


def test_salidas_y_figuras(cfg_chica, ranking, tmp_path):
    exportar.exportar_ranking(cfg_chica, ranking, tmp_path)
    leido = exportar.leer_ranking(tmp_path / "ranking.csv")
    assert len(leido) == len(ranking) and leido["factible"].dtype == bool
    assert (tmp_path / "pareto.csv").exists() and (tmp_path / "tolerancias_implicitas.json").exists()
    rutas = exportar.figuras_barrido(cfg_chica, ranking, tmp_path / "figs")
    nombres = {r.name for r in rutas}
    assert {"pareto.png", "factibilidad.png", "ganador_perfil.png", "ganador_dibujo.png",
            "sensibilidad.png"} <= nombres
