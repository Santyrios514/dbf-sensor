from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

RAIZ = Path(__file__).resolve().parents[1]
CONFIG = RAIZ / "config" / "config.yaml"
ORK = RAIZ / "modelos" / "analisis_vol_int.ork"


@pytest.fixture
def raw() -> dict:
    """Configuración base del repositorio como dict (se puede modificar libremente)."""
    with open(CONFIG, encoding="utf-8") as f:
        return copy.deepcopy(yaml.safe_load(f))


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
