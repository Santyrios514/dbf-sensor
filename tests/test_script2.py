"""Script 2 de punta a punta y esquemas de los CSV (v1 §7 + v2 §7)."""

import subprocess
import sys

import pandas as pd
import pytest

from sensor_lastre import esquemas
from sensor_lastre.analisis import analizar_caso
from sensor_lastre.config import cargar
from sensor_lastre.masas import masas_por_componente
from sensor_lastre.perfiles import muestrear

from conftest import CONFIG, RAIZ

SCRIPT = RAIZ / "scripts" / "02_volumen_lastre.py"
MM = 1e-3


def _correr(*args, config=CONFIG):
    r = subprocess.run([sys.executable, str(SCRIPT), "--config", str(config), "--sin-figuras", *args],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


@pytest.fixture(scope="module")
def salida_sin_orlab(tmp_path_factory):
    d = tmp_path_factory.mktemp("datos")
    _correr("--sin-orlab", "--dir-datos", str(d))
    return d


def test_esquemas(salida_sin_orlab):
    d = salida_sin_orlab
    for nombre, cols in (("resultados.csv", esquemas.RESULTADOS), ("masas_capas.csv", esquemas.MASAS_CAPAS),
                         ("optimo.csv", esquemas.OPTIMO),
                         ("presupuestos.csv", esquemas.PRESUPUESTOS), ("ventanas.csv", esquemas.VENTANAS)):
        assert list(pd.read_csv(d / nombre).columns) == cols, nombre
    curvas = list((d / "curvas").glob("*.csv"))
    assert curvas
    for c in curvas:
        assert list(pd.read_csv(c).columns) == esquemas.CURVA


def test_esquema_sin_columnas_trapezoidales():
    for col in ("aletas_cr_mm", "aletas_ct_mm", "aletas_s_mm", "aletas_flecha_mm"):
        assert col not in esquemas.OR_RESUMEN


def test_cobertura_barrido(salida_sin_orlab):
    cfg = cargar(CONFIG)
    n = len(cfg.casos) * len(cfg.rellenos)
    res = pd.read_csv(salida_sin_orlab / "resultados.csv")
    assert len(res) == n == 25 * 5
    assert set(res["x_CP_fuente"]) == {"interno"}
    pres = pd.read_csv(salida_sin_orlab / "presupuestos.csv")
    assert len(pres) == n * 5
    banderas = {b for s in res["banderas"].fillna("") for b in s.split(";") if b}
    assert banderas <= set(esquemas.BANDERAS)
    capas = pd.read_csv(salida_sin_orlab / "masas_capas.csv")
    m_por_caso = capas.groupby("config_id")["m_g"].sum()
    m_res = res.groupby("config_id")["m_casco_g"].first()
    assert (m_por_caso - m_res).abs().max() < 1e-3


def _simular_script1(caso, cfg, d, dx_CP_mm=3.0):
    """Escribe los 4 CSV del script 1 a partir del modelo analítico (sin JVM)."""
    g = caso.geom
    muestrear(g, cfg.openrocket["n_muestras_perfil"], caso.id)[esquemas.OR_PERFILES].to_csv(
        d / "or_perfiles.csv", index=False)
    a = g.aletas
    pd.DataFrame({"config_id": caso.id, "orden": range(len(a.poligono)), "x_mm": a.poligono[:, 0] / MM,
                  "r_mm": a.poligono[:, 1] / MM, "origen": a.origen}).to_csv(d / "or_aletas.csv", index=False)
    res = analizar_caso(caso, rellenos=[])
    comp = masas_por_componente(caso, res.mv)
    pd.DataFrame([{"config_id": caso.id, "componente": k, "m_g": m / MM, "x_CG_mm": x / MM}
                  for k, (m, x) in comp.items()]).reindex(columns=esquemas.OR_COMPONENTES).to_csv(
        d / "or_componentes.csv", index=False)
    pd.DataFrame([{"config_id": caso.id, "estado": "ok", "x_CP_mm": res.cp_int.x_CP / MM + dx_CP_mm,
                   "CN_alpha_rad": res.cp_int.CNa, "aletas_x_r0_mm": a.x_LE / MM, "r_LE_mm": a.r_LE / MM,
                   "D_popa_mm": 2 * g.cola.Ra / MM}]).reindex(columns=esquemas.OR_RESUMEN).to_csv(
        d / "or_resumen.csv", index=False)


def test_modo_csv_openrocket(tmp_path):
    cfg = cargar(CONFIG)
    caso = cfg.caso("base_ork")
    _simular_script1(caso, cfg, tmp_path)
    _correr("--sin-orlab", "--solo", caso.id, "--dir-datos", str(tmp_path / "base"))
    ref = pd.read_csv(tmp_path / "base" / "resultados.csv").iloc[0]
    args = [f"--or-{k}={tmp_path / f'or_{k}.csv'}" for k in ("resumen", "perfiles", "aletas", "componentes")]
    _correr(*args, "--solo", caso.id, "--dir-datos", str(tmp_path / "or"))
    r = pd.read_csv(tmp_path / "or" / "resultados.csv").iloc[0]
    assert r["x_CP_fuente"] == "openrocket"
    assert r["dx_CP_or_vs_interno_mm"] == pytest.approx(3.0, abs=0.05)
    assert r["V_ext_cm3"] == pytest.approx(ref["V_ext_cm3"], rel=1e-3)
    assert r["V_util_geo_cm3"] == pytest.approx(ref["V_util_geo_cm3"], rel=2e-3)
    assert r["A_aleta_mm2"] == pytest.approx(ref["A_aleta_mm2"], rel=1e-9)
    assert abs(r["dif_masa_or_pct"]) < 0.5
    assert "dif_CP_alta" not in str(r["banderas"]) and "dif_masa_alta" not in str(r["banderas"])


def test_faltan_csv_de_openrocket(tmp_path):
    r = subprocess.run([sys.executable, str(SCRIPT), "--config", str(CONFIG), "--sin-figuras",
                        "--dir-datos", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 2 and "--sin-orlab" in r.stderr
