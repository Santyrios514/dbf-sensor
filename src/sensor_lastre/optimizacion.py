"""Criterio de optimización del equipo: **máxima masa de lastre con SM ≥ SM_min**.

Para cada configuración × relleno, `analisis` ya da la masa máxima admisible (`m_relleno_g`):

1. el tapón delantero crece hasta ℓ_SM,max (lo que el SM o la geometría permitan);
2. si llega a ℓ_geo y `lastre.lastre_trasero` está activo, se agrega el tapón trasero detrás
   de la electrónica hasta que SM = SM_min o se acaba la cavidad.

Aquí solo se filtran las configuraciones que caben en la bahía y se ordenan por esa masa. Con la
cavidad llena, varias configuraciones empatan; entonces se prefiere la de más SM (margen ante la
incertidumbre del CP) y luego la que ocupa menos guardada (D acostado, luego D parado).
"""

from __future__ import annotations

import pandas as pd

from .esquemas import OPTIMO


def ranking(resultados: pd.DataFrame, exigir: list[str] | tuple[str, ...] = ()) -> pd.DataFrame:
    """Configuraciones factibles ordenadas por masa de lastre, por relleno (puesto 1 = óptimo).

    Factible: m_relleno_g finito (cumple SM_min) y cada bandera de `exigir` verdadera o sin dato
    (una bandera vacía significa que falta la dimensión de la bahía, no que no quepa).
    """
    df = resultados.copy()
    ok = df["m_relleno_g"].notna()
    for col in exigir:
        v = df[col]
        ok &= v.isna() | v.astype(str).str.lower().isin(["true", "1"])
    df = df[ok].assign(_m=lambda d: d["m_relleno_g"].round(0))  # empate = misma masa al gramo
    claves = [(c, a) for c, a in (("relleno", True), ("_m", False), ("SM_cal", False),
                                  ("D_acostado_mm", True), ("D_parado_mm", True)) if c in df]
    df = df.sort_values([c for c, _ in claves], ascending=[a for _, a in claves], kind="stable")
    df["puesto"] = df.groupby("relleno").cumcount() + 1
    return df.reindex(columns=OPTIMO)
