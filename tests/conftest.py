from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

RAIZ = Path(__file__).resolve().parents[1]
CONFIG = RAIZ / "config" / "config.yaml"


@pytest.fixture
def raw() -> dict:
    """Configuración base del repositorio como dict (se puede modificar libremente)."""
    with open(CONFIG, encoding="utf-8") as f:
        return copy.deepcopy(yaml.safe_load(f))


def solo(raw: dict, L_mm: float, D_mm: float) -> dict:
    """Restringe el barrido a una sola configuración."""
    raw["barrido"] = {"L_mm": [L_mm], "D_mm": [D_mm]}
    raw["configuraciones"] = {}
    return raw
