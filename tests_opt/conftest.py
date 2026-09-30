from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

RAIZ = Path(__file__).resolve().parents[1]
CONFIG_OPT = RAIZ / "config" / "optimizacion.yaml"


def raw_opt() -> dict:
    with open(CONFIG_OPT, encoding="utf-8") as f:
        return copy.deepcopy(yaml.safe_load(f))


def malla_chica(raw: dict, **cambios) -> dict:
    """Malla reducida (1 D, 1 forma, 2 L_t, 2 k y pocas aletas) para tests rápidos."""
    raw["malla"] = {
        "D_mm": [65], "cola_forma": [{"forma": "parabolica", "parametro": 1.0}],
        "L_t_rel_D": [1.0, 1.5], "k": [0.3, 0.5], "lambda_LE": [0.0, 0.4], "mu_cr": [0.75, 1.0],
        "r_tip_rel_R": [1.6, 2.0, 2.4], "gamma_ct": [0.5], "sigma_flecha": [0.0, 1.0],
    }
    raw["malla"].update(cambios)
    return raw


@pytest.fixture
def raw() -> dict:
    return raw_opt()


@pytest.fixture(scope="session")
def cfg_default():
    from sensor_opt.config import cargar
    return cargar(CONFIG_OPT)
