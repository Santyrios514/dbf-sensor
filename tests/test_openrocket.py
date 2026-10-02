"""Fase C2/D2 (v2 §8.2, §8.5 y v1 §8.3): pruebas contra OpenRocket. Requieren JVM y el jar de
OpenRocket 24.12 (`python -m orlab fetch 24.12` o ORLAB_JAR); si no están, se saltan."""

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from sensor_lastre import esquemas
from sensor_lastre.analisis import analizar_caso
from sensor_lastre.config import cargar
from sensor_lastre.masas import masas_por_componente
from sensor_lastre.perfiles import FORMAS, forma_or, radio_cola, radio_exterior

from conftest import CONFIG, ORK, RAIZ

MM = 1e-3


def _hay_openrocket() -> bool:
    import os
    import shutil
    try:
        import jpype  # noqa: F401
        from orlab.jars import jar_cache_dir
    except Exception:  # noqa: BLE001
        return False
    if shutil.which("java") is None:
        return False
    return bool(os.environ.get("ORLAB_JAR")) or any(jar_cache_dir().glob("OpenRocket-*.jar"))


pytestmark = pytest.mark.skipif(not _hay_openrocket(), reason="sin JVM/jar de OpenRocket")


@pytest.fixture(scope="module")
def puente():
    import orlab
    from sensor_lastre.or_bridge import PuenteOR
    cfg = cargar(CONFIG)
    with orlab.OpenRocketInstance() as inst:
        yield PuenteOR(inst, cfg.openrocket["nombres_componentes"]).cargar(ORK), cfg


def test_regresion_seccion_1_4():
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "02_validar_openrocket.py"), "--regresion"],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    assert r.stdout.count("OK") == 4


def test_base_leida_por_openrocket_coincide(puente):
    from sensor_lastre.verificacion import comparar
    pr, cfg = puente
    pr.cargar(ORK)
    assert comparar(pr.valores(), cfg.caso("base_ork"), 1e-6) == []


@pytest.mark.parametrize("recortada", [True, False])
@pytest.mark.parametrize("forma, p", [("conica", None), ("ogiva", 1.0), ("ogiva", 0.5), ("elipsoide", None),
                                      ("potencia", 0.5), ("parabolica", 0.5), ("parabolica", 1.0),
                                      ("haack", 0.0), ("haack", 1 / 3)])
def test_perfil_transicion_vs_openrocket(puente, forma, p, recortada):
    """v2 §3.1/§8.2: el perfil analítico coincide con el muestreado de OpenRocket (< 0.05 mm),
    con y sin recorte (clipped)."""
    pr, _ = puente
    pr.cargar(ORK)
    t = pr.comp["cola"]
    t.setShapeType(pr.core.rocketcomponent.Transition.Shape.valueOf(FORMAS[forma][0]))
    if p is not None:
        t.setShapeParameter(float(p))
    t.setClipped(recortada)
    Lt, Rf, Ra = float(t.getLength()), float(t.getForeRadius()), float(t.getAftRadius())
    u = np.linspace(0, Lt, 201)
    r_or = np.array([float(t.getRadius(float(x))) for x in u])
    r_an = radio_cola(u, forma, Rf, Ra, Lt, p, recortada)
    assert np.max(np.abs(r_or - r_an)) < 0.05 * MM


@pytest.mark.parametrize("forma, p", [("conica", None), ("ogiva", 1.0), ("elipsoide", None),
                                      ("potencia", 0.5), ("parabolica", 0.75), ("haack", 0.0)])
def test_perfil_nariz_vs_openrocket(puente, forma, p):
    pr, _ = puente
    pr.cargar(ORK)
    n = pr.comp["nariz"]
    n.setShapeType(pr.core.rocketcomponent.Transition.Shape.valueOf(FORMAS[forma][0]))
    if p is not None:
        n.setShapeParameter(float(p))
    L, R = float(n.getLength()), float(n.getBaseRadius())
    x = np.linspace(0, L, 201)
    r_or = np.array([float(n.getRadius(float(v))) for v in x])
    assert np.max(np.abs(r_or - forma_or(x, forma, R, L, p))) < 0.05 * MM


@pytest.fixture(scope="module")
def export(tmp_path_factory):
    d = tmp_path_factory.mktemp("or")
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "02_validar_openrocket.py"), "--todos",
                        "--config", str(CONFIG), "--solo", "base_ork,L350_D65_h_tip30_extension20_rotacion_deg0,L350_D65_h_tip40_extension0_rotacion_deg0",
                        "--dir-datos", str(d)], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    return d


def test_esquemas_script1(export):
    for nombre, cols in (("or_resumen.csv", esquemas.OR_RESUMEN), ("or_perfiles.csv", esquemas.OR_PERFILES),
                         ("or_componentes.csv", esquemas.OR_COMPONENTES), ("or_aletas.csv", esquemas.OR_ALETAS)):
        assert list(pd.read_csv(export / nombre).columns) == cols, nombre
    res = pd.read_csv(export / "or_resumen.csv")
    assert set(res["estado"]) == {"ok"} and res["or_version"].astype(str).eq("24.12").all()
    assert res["mach_usado"].iloc[0] == pytest.approx(0.0897, abs=1e-3)  # no 0.3
    # la advertencia de perfil irregular (escalón P4–P5) se conserva; sin extensión no hay escalón
    jag = res.set_index("config_id")["or_warnings"].fillna("").str.contains("Jagged")
    assert jag["base_ork"] and jag["L350_D65_h_tip30_extension20_rotacion_deg0"]
    assert not jag["L350_D65_h_tip40_extension0_rotacion_deg0"]


def test_integracion_vs_modelo_interno(export):
    """v1 §8.3 / v2 §6: CP < 5 % L, masas ±3 % (aletas ±2 %), perfil < 0.1 mm."""
    cfg = cargar(CONFIG)
    res = pd.read_csv(export / "or_resumen.csv").set_index("config_id")
    comp = pd.read_csv(export / "or_componentes.csv")
    perf = pd.read_csv(export / "or_perfiles.csv")
    for cid in res.index:
        caso = cfg.caso(cid)
        r = analizar_caso(caso, rellenos=[])
        assert abs(res.loc[cid, "x_CP_mm"] * MM - r.cp_int.x_CP) < 0.05 * caso.geom.L
        propias = masas_por_componente(caso, r.mv)
        c = comp[comp.config_id == cid].set_index("componente")
        for k, (m, _) in propias.items():
            lim = 0.02 if k == "aletas" else 0.03
            assert m / (c.loc[k, "m_g"] * 1e-3) - 1 == pytest.approx(0, abs=lim), k
        p = perf[perf.config_id == cid]
        r_an = radio_exterior(p["x_mm"].to_numpy() * MM, caso.geom)
        assert np.max(np.abs(r_an - p["r_ext_mm"].to_numpy() * MM)) < 0.1 * MM
        # r_LE leído de OpenRocket = r_LE analítico (P5 sobre la superficie)
        assert res.loc[cid, "r_LE_mm"] * MM == pytest.approx(caso.geom.aletas.r_LE, abs=0.01 * MM)
        assert res.loc[cid, "h_env_mm"] * MM == pytest.approx(
            2 * (caso.geom.aletas.r_tip), abs=0.01 * MM)


def test_poligono_or_igual_al_analitico(export):
    cfg = cargar(CONFIG)
    ale = pd.read_csv(export / "or_aletas.csv")
    for cid, d in ale.groupby("config_id"):
        g = cfg.caso(cid).geom.aletas
        pts = d[d.origen == "punto"][["x_mm", "r_mm"]].to_numpy() * MM
        assert np.max(np.abs(pts - (g.puntos + [g.x_LE, g.r_LE]))) < 1e-6
        from sensor_lastre.aletas import area_centroide
        A, _ = area_centroide(d.sort_values("orden")[["x_mm", "r_mm"]].to_numpy() * MM)
        assert A == pytest.approx(g.area, rel=2e-3)  # la raíz de OpenRocket es de baja resolución


def test_puntos_invalidos_se_detectan(puente):
    """OpenRocket no lanza excepción con puntos inválidos: los modifica. aplicar() lo detecta."""
    from sensor_lastre.or_bridge import ErrorAleta
    pr, cfg = puente
    pr.cargar(ORK)
    caso = cfg.caso("base_ork")
    from dataclasses import replace
    # P4 dentro del cuerpo (ε negativo): OpenRocket lo pega a la superficie sin avisar
    malo = replace(caso, params_aleta=replace(caso.params_aleta, eps=-0.02))
    with pytest.raises(ErrorAleta):
        pr.aplicar(malo)
