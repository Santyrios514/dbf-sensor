"""Script 2 de punta a punta y esquemas de los CSV (sección 7)."""

import subprocess
import sys

import pandas as pd
import pytest

from sensor_lastre import esquemas
from sensor_lastre.config import cargar
from sensor_lastre.perfiles import muestrear

from conftest import CONFIG, RAIZ

SCRIPT = RAIZ / "scripts" / "02_volumen_lastre.py"


def _correr(*args):
    r = subprocess.run([sys.executable, str(SCRIPT), "--config", str(CONFIG), "--sin-figuras", *args],
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
                         ("presupuestos.csv", esquemas.PRESUPUESTOS), ("ventanas.csv", esquemas.VENTANAS)):
        assert list(pd.read_csv(d / nombre).columns) == cols, nombre
    curvas = list((d / "curvas").glob("*.csv"))
    assert curvas
    for c in curvas:
        assert list(pd.read_csv(c).columns) == esquemas.CURVA


def test_cobertura_barrido(salida_sin_orlab):
    cfg = cargar(CONFIG)
    res = pd.read_csv(salida_sin_orlab / "resultados.csv")
    assert len(res) == len(cfg.casos) * len(cfg.rellenos) == 45
    assert set(res["x_CP_fuente"]) == {"interno"}
    pres = pd.read_csv(salida_sin_orlab / "presupuestos.csv")
    assert len(pres) == 45 * 5
    banderas = {b for s in res["banderas"].fillna("") for b in s.split(";") if b}
    assert banderas <= set(esquemas.BANDERAS)
    capas = pd.read_csv(salida_sin_orlab / "masas_capas.csv")
    m_por_caso = capas.groupby("config_id")["m_g"].sum()
    m_res = res.groupby("config_id")["m_casco_g"].first()
    assert (m_por_caso - m_res).abs().max() < 1e-3


def test_modo_csv_openrocket(tmp_path):
    """Simula la salida del script 1 con el perfil analítico y un x_CP desplazado 3 mm."""
    cfg = cargar(CONFIG)
    caso = cfg.caso("L350_D65")
    perf = muestrear(caso.geom, cfg.openrocket["n_muestras_perfil"], caso.id)
    perf[esquemas.OR_PERFILES].to_csv(tmp_path / "or_perfiles.csv", index=False)
    base = _correr("--sin-orlab", "--solo", caso.id, "--dir-datos", str(tmp_path / "base"))
    ref = pd.read_csv(tmp_path / "base" / "resultados.csv").iloc[0]
    resumen = pd.DataFrame([{"config_id": caso.id, "estado": "ok", "x_CP_mm": ref["x_CP_mm"] + 3.0,
                             "CN_alpha_rad": ref["CN_alpha_rad"]}]).reindex(columns=esquemas.OR_RESUMEN)
    resumen.to_csv(tmp_path / "or_resumen.csv", index=False)
    _correr("--or-resumen", str(tmp_path / "or_resumen.csv"), "--or-perfiles",
            str(tmp_path / "or_perfiles.csv"), "--solo", caso.id, "--dir-datos", str(tmp_path / "or"))
    r = pd.read_csv(tmp_path / "or" / "resultados.csv").iloc[0]
    assert r["x_CP_fuente"] == "openrocket"
    assert r["dx_CP_or_vs_interno_mm"] == pytest.approx(3.0, abs=0.05)
    # el perfil muestreado (800 pts/componente) reproduce el volumen analítico
    assert r["V_ext_cm3"] == pytest.approx(ref["V_ext_cm3"], rel=1e-3)
    assert r["V_util_geo_cm3"] == pytest.approx(ref["V_util_geo_cm3"], rel=2e-3)
    assert "dif_CP_alta" not in str(r["banderas"])
