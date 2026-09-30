"""Fase O4 (requiere JVM y jar de OpenRocket; si no están, se salta). Las corridas con Java van en
subprocesos: JPype no puede reiniciar la JVM dentro del mismo proceso de pytest."""

import json
import os
import shutil
import subprocess
import sys
import textwrap

import pandas as pd
import pytest
import yaml

from .conftest import CONFIG_OPT, RAIZ, malla_chica, raw_opt


def _hay_openrocket() -> bool:
    try:
        import jpype  # noqa: F401
        from orlab.jars import jar_cache_dir
    except Exception:  # noqa: BLE001
        return False
    if shutil.which("java") is None:
        return False
    return bool(os.environ.get("ORLAB_JAR")) or any(jar_cache_dir().glob("OpenRocket-*.jar"))


pytestmark = pytest.mark.skipif(not _hay_openrocket(), reason="sin JVM/jar de OpenRocket")


def test_T7_perfil_de_la_cola_coincide_con_openrocket():
    codigo = textwrap.dedent(f"""
        import sys; sys.path.insert(0, {str(RAIZ / 'src')!r})
        import orlab
        from sensor_opt.config import cargar
        from sensor_opt.verificacion import PuenteOpt, comparar_perfiles
        cfg = cargar({str(CONFIG_OPT)!r})
        orc = cfg.base_raw["openrocket"]
        with orlab.OpenRocketInstance(log_level="ERROR") as inst:
            pr = PuenteOpt(inst, orc["nombres_componentes"], cfg.ruta_base(orc["ork_base"]))
            print(comparar_perfiles(cfg, pr).to_json(orient="records"))
    """)
    r = subprocess.run([sys.executable, "-c", codigo], capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr[-3000:]
    df = pd.DataFrame(json.loads(r.stdout.strip().splitlines()[-1]))
    assert len(df) == 4 * 2 * 2  # 4 formas × k extremos × L_t extremos
    assert (df["dif_max_mm"] < 0.05).all(), df.sort_values("dif_max_mm").tail()


@pytest.fixture(scope="module")
def corrida(tmp_path_factory):
    d = tmp_path_factory.mktemp("opt")
    raw = malla_chica(raw_opt())
    raw["base_config"] = str(RAIZ / "config" / "config.yaml")
    raw["salida"] = {"dir": str(d / "data"), "dir_figuras": str(d / "figs"), "exportar_ork_top": 2}
    raw["ejecucion"]["verificacion_or"] = {"N_verif": 5, "N_cal": 8, "tol_cp_mm": 2.0, "max_iter_calibracion": 2}
    ruta = d / "opt.yaml"
    ruta.write_text(yaml.safe_dump(raw, allow_unicode=True))
    for s in ("04_optimizar_malla.py", "05_verificar_openrocket.py"):
        r = subprocess.run([sys.executable, str(RAIZ / "scripts" / s), "--config", str(ruta), "--procesos", "2"],
                           capture_output=True, text=True, timeout=1800)
        assert r.returncode == 0, s + "\n" + r.stdout[-3000:] + r.stderr[-3000:]
    return d


def test_T10_ganador_verificado_y_ork_reproducible(corrida):
    rk = pd.read_csv(corrida / "data" / "ranking.csv")
    g = rk[rk["ganador"].astype(str).str.lower() == "true"]
    assert len(g) == 1
    g = g.iloc[0]
    assert str(g["verificado_or"]).lower() == "true"
    assert 1.0 <= g["SM_or_cal"] <= 2.0
    cal = json.loads((corrida / "data" / "calibracion.json").read_text())
    assert cal["ganador"] == g["cand_id"] and cal["reapertura_ganador_dif_CP_mm"] <= 0.5
    orks = sorted((corrida / "data" / "ork").glob("rank*.ork"))
    assert orks and orks[0].name == f"rank01_{g['cand_id']}.ork"


def test_verificacion_or_y_calibracion(corrida):
    ver = pd.read_csv(corrida / "data" / "verificacion_or.csv")
    assert {"verificacion", "calibracion"} <= set(r.split("+")[0] for r in ver["rol"])
    ok = ver[ver["factible_or"].astype(str).str.lower() == "true"]
    assert (ok["dif_masa_or_pct"].abs() < 1.0).all()  # masas y CG del modelo = OpenRocket con Mass components
    cal = json.loads((corrida / "data" / "calibracion.json").read_text())
    assert 0.8 < cal["beta"] < 1.2 and cal["R2"] > 0.9
    for f in ("calibracion.png", "ganador_perfil.png", "ganador_dibujo.png", "pareto.png"):
        assert (corrida / "figs" / f).exists()
