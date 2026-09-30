"""Validación con OpenRocket de los ganadores del modelo propio (script 02).

La lógica de selección se prueba sin JVM con un exportador falso que devuelve el modelo analítico
con el CP corrido dx; la corrida real (scripts 01 → 02) se prueba con OpenRocket si está."""

import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from sensor_lastre.analisis import analizar_caso
from sensor_lastre.config import cargar
from sensor_lastre.masas import masas_por_componente
from sensor_lastre.optimizacion import ranking
from sensor_lastre.perfiles import muestrear
from sensor_lastre.validacion_or import ExportOR, _puede_superar, validar

from conftest import CONFIG, RAIZ, barrido

MM = 1e-3


def _exportador_falso(dx_mm: float, llamados: list):
    def exportar(caso):
        llamados.append(caso.id)
        g, a = caso.geom, caso.geom.aletas
        res = analizar_caso(caso, rellenos=[])
        comp = masas_por_componente(caso, res.mv)
        resumen = {"config_id": caso.id, "estado": "ok", "x_CP_mm": res.cp_int.x_CP / MM + dx_mm,
                   "CN_alpha_rad": res.cp_int.CNa, "aletas_x_r0_mm": a.x_LE / MM, "r_LE_mm": a.r_LE / MM,
                   "D_popa_mm": 2 * g.cola.Ra / MM, "or_warnings": ""}
        perf = muestrear(g, 800, caso.id).to_dict("records")
        ale = [{"config_id": caso.id, "orden": i, "x_mm": x / MM, "r_mm": r / MM, "origen": o}
               for i, ((x, r), o) in enumerate(zip(a.poligono, a.origen))]
        comps = [{"config_id": caso.id, "componente": k, "m_g": m / MM, "x_CG_mm": x / MM}
                 for k, (m, x) in comp.items()]
        return ExportOR(resumen=resumen, perfiles=perf, componentes=comps, aletas=ale)
    return exportar


@pytest.fixture(scope="module")
def cfg_y_resultados():
    import yaml
    raw = yaml.safe_load(open(CONFIG, encoding="utf-8"))
    barrido(raw, 350, [60, 65], aletas__h_tip__valor=[20, 30, 40], aletas__extension__valor=[0, 40])
    raw["lastre"]["masa_max_sensor_g"] = None  # sin tope: el SM decide la masa en las aletas chicas
    cfg = cargar(raw)
    filas = [rr.fila for c in cfg.casos for rr in analizar_caso(c).rellenos]
    return cfg, pd.DataFrame(filas)


def test_sin_diferencias_valida_solo_los_primeros(cfg_y_resultados):
    cfg, res = cfg_y_resultados
    llamados = []
    v = validar(cfg, res, _exportador_falso(0.0, llamados), n_verif=2, agrupar=[], max_iter=3, log=lambda *_: None)
    top = ranking(res, [])
    assert v.ganadores.iloc[0]["config_id"] == top.iloc[0]["config_id"]
    assert len(llamados) == 2  # el ganador validado no puede ser superado por ningún otro
    t = v.tabla
    assert t["ganador"].sum() == 1 and np.allclose(t["dx_CP_mm"], 0, atol=0.05)
    assert np.allclose(t["m_total_or_g"], t["m_total_int_g"], rtol=2e-3)


def test_cp_adelantado_en_openrocket_obliga_a_validar_mas(cfg_y_resultados):
    """Con el CP de OpenRocket 40 mm más adelante, los primeros internos pierden masa (el SM
    manda) y hay que seguir validando hasta que ningún no validado pueda superar al ganador."""
    cfg, res = cfg_y_resultados
    llamados = []
    v = validar(cfg, res, _exportador_falso(-40.0, llamados), n_verif=1, agrupar=[], max_iter=20,
                log=lambda *_: None)
    assert len(llamados) > 1
    g = v.ganadores.iloc[0]
    r_or = ranking(v.resultados_or, [])
    assert g["config_id"] == r_or.iloc[0]["config_id"]  # el mejor entre los validados, con OpenRocket
    no_val = ranking(res, []).query("config_id not in @llamados")
    assert not any(_puede_superar(f, g) for _, f in no_val.iterrows())


def test_agrupar_por_par_L_D(cfg_y_resultados):
    cfg, res = cfg_y_resultados
    v = validar(cfg, res, _exportador_falso(0.0, []), n_verif=1, agrupar=["L_mm", "D_mm"], max_iter=2,
                log=lambda *_: None)
    assert sorted(v.ganadores["D_mm"]) == [60, 65] and (v.ganadores["L_mm"] == 350).all()


def test_agrupar_por_clave_inexistente(cfg_y_resultados):
    cfg, res = cfg_y_resultados
    with pytest.raises(KeyError):
        validar(cfg, res, _exportador_falso(0.0, []), n_verif=1, agrupar=["no_existe"], max_iter=1)


def _hay_openrocket() -> bool:
    import os
    try:
        import jpype  # noqa: F401
        from orlab.jars import jar_cache_dir
    except Exception:  # noqa: BLE001
        return False
    return shutil.which("java") is not None and (bool(os.environ.get("ORLAB_JAR"))
                                                 or any(jar_cache_dir().glob("OpenRocket-*.jar")))


@pytest.mark.skipif(not _hay_openrocket(), reason="sin JVM/jar de OpenRocket")
def test_scripts_01_y_02_de_punta_a_punta(tmp_path):
    d, f = tmp_path / "datos", tmp_path / "figs"
    for s, extra in (("01_volumen_lastre.py", ["--sin-figuras"]), ("02_validar_openrocket.py", ["--n-verif", "2"])):
        r = subprocess.run([sys.executable, str(RAIZ / "scripts" / s), "--config", str(CONFIG),
                            "--dir-datos", str(d), "--dir-figuras", str(f), *extra],
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, s + "\n" + r.stdout[-3000:] + r.stderr[-3000:]
    t = pd.read_csv(d / "validacion_or.csv")
    gan = pd.read_csv(d / "ganadores.csv")
    assert len(gan) == 1 and t["ganador"].sum() == 1
    assert (t["estado_or"] == "ok").all() and (t["dx_CP_mm"].abs() < 0.05 * 350).all()
    assert gan.iloc[0]["x_CP_fuente"] == "openrocket" and gan.iloc[0]["SM_cal"] >= 1.0 - 0.005
    exp = pd.read_csv(d / "or_resumen.csv")
    assert set(exp["config_id"]) == set(t["config_id"])  # solo se exportan las validadas
    assert len(pd.read_csv(d / "resultados.csv")) > len(exp)
    assert any(p.name.startswith("ganador_or__") for p in f.iterdir())
