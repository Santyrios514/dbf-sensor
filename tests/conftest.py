from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

RAIZ = Path(__file__).resolve().parents[1]
CONFIG_REPO = RAIZ / "config" / "config.yaml"
ORK = RAIZ / "modelos" / "analisis_vol_int.ork"


ELECTRONICA_DETRAS = {"longitud_mm": 100, "masa_g": 150, "modo": "detras_del_lastre", "x_inicio_mm": None,
                      "holgura_mm": 2.0, "x_cg_relativo_mm": None}


def _config_tests() -> Path:
    """config.yaml con los valores con que se escribieron los tests del flujo principal.

    Los tests comprueban el código, no el diseño vigente: se fijan la electrónica detrás del lastre
    (100 mm), el presupuesto de 6.2 kg, la fracción máxima de 0.6 L y la densidad de pared del
    `.ork` base (la verificación contra el `.ork` la compara). El archivo vive en tests/ para que
    las rutas relativas del config se resuelvan igual que desde config/.
    """
    raw = yaml.safe_load(CONFIG_REPO.read_text(encoding="utf-8"))
    raw["electronica"] = dict(ELECTRONICA_DETRAS)
    raw["lastre"].update({"fraccion_max_L": 0.6, "masa_max_sensor_g": 6200})
    raw["materiales"]["MAT_PARED_EQ"] = 1342  # valor del .ork base
    ruta = RAIZ / "tests" / "_config_tests.yaml"
    ruta.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return ruta


CONFIG = _config_tests()


@pytest.fixture
def raw() -> dict:
    """Configuración base del repositorio como dict (se puede modificar libremente)."""
    with open(CONFIG, encoding="utf-8") as f:
        return copy.deepcopy(yaml.safe_load(f))


PERDIGON = {"nombre": "perdigon_pb_epoxy", "material": "plomo",
            "granular": {"empaquetamiento": 0.62, "matriz": "epoxy", "costo_grano_usd_kg": 6}}


def con_rellenos(raw: dict, *extra: dict) -> dict:
    """Agrega rellenos al del equipo (solo plomo) para probar granulares, costos, etc."""
    raw["rellenos"] = list(raw["rellenos"]) + [dict(r) for r in extra]
    return raw


def barrido(raw: dict, L_mm, D_mm, **rutas) -> dict:
    """Barrido cartesiano de L × D (listas o escalares) más rutas extra; sin overrides."""
    def lista(v):
        return list(v) if isinstance(v, (list, tuple)) else [v]
    params = {"L_mm": lista(L_mm), "D_mm": lista(D_mm)}
    params.update({k.replace("__", "."): lista(v) for k, v in rutas.items()})
    raw["barrido"] = {"modo": "cartesiano", "parametros": params}
    raw["configuraciones"] = {}
    return raw


def solo(raw: dict, L_mm: float, D_mm: float) -> dict:
    """Restringe el barrido a una sola configuración (id L{L}_D{D})."""
    return barrido(raw, L_mm, D_mm)
