"""Llenado de lastre hasta m_max (v3 §3.5): se reutiliza tal cual el de `sensor_lastre`.

El código actual ya implementa el orden que exige la spec (solución del programa lineal: llenar
en x creciente, contiguo desde la nariz, y pasar a la zona trasera solo con la delantera llena):

1. tapón delantero desde x_b0 hasta ℓ_SM,max (lo que permitan el SM o la geometría);
2. si llegó a ℓ_geo, tapón trasero detrás de la electrónica hasta SM = SM_min o cavidad llena;
3. si la masa total supera m_max (`lastre.masa_max_sensor_g`), el lastre se recorta a
   m_max − m_vacío y se ubica igual, delantero primero.

Aquí solo se arma el `ResultadoCaso` del candidato (cavidad y masas del casco del caché del
cuerpo, aleta del candidato, CP del sustituto o calibrado) y se llama a `_analizar_relleno`.
El SM_max no entra al llenado: si con el lastre máximo admisible el SM supera SM_max, el
candidato es infactible por sobreestabilidad (lo resuelven aletas más pequeñas).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from sensor_lastre.aletas import GeomAleta
from sensor_lastre.analisis import ResultadoCaso, ResultadoRelleno, _analizar_relleno
from sensor_lastre.barrowman import ResultadoCP

from .geometria import Cuerpo

G = 1e-3
MM = 1e-3

# limitante_masa de sensor_lastre → restriccion_activa de la spec v3
RESTRICCION = {"masa_max_sensor": "masa", "SM": "SM_min", "geometria": "volumen"}
# banderas del código actual que invalidan el candidato (restricciones heredadas, v3 §3.7.6)
BANDERAS_INFACTIBLES = ("inviable_geo", "SM_inalcanzable", "SM_inalcanzable_por_geometria",
                        "SM_inalcanzable_con_masa_max", "inestable_respecto_remolque", "trim_no_lineal")


@dataclass
class Llenado:
    res: ResultadoCaso
    rr: ResultadoRelleno
    restriccion_activa: str  # masa | SM_min | volumen | ninguna
    ell: float  # longitud del tapón delantero [m]; NaN si no hay llenado válido
    m_delantero: float  # kg
    motivos: list[str]

    @property
    def fila(self) -> dict:
        return self.rr.fila

    @property
    def ok(self) -> bool:
        return not self.motivos


def caso_candidato(cu: Cuerpo, g: GeomAleta):
    return replace(cu.caso.con_aletas(g), params_aleta=g.params)


def llenar(cu: Cuerpo, g: GeomAleta, cp: ResultadoCP, SM_max: float, fuente: str = "sustituto") -> Llenado:
    caso = caso_candidato(cu, g)
    mv = replace(cu.mv, m_aletas=g.masa, x_aletas=g.x_cg)
    res = ResultadoCaso(caso=caso, cav=cu.cav, mv=mv, cp_int=cp, x_CP=cp.x_CP, CNa=cp.CNa,
                        x_CP_fuente=fuente, x_b0=cu.x_b0, lim=cu.lim, banderas=[])
    rr = _analizar_relleno(res, caso.rellenos[0])
    f = rr.fila
    banderas = [b for b in (f.get("banderas") or "").split(";") if b]
    motivos = [b for b in banderas if b in BANDERAS_INFACTIBLES]
    m_rel = f.get("m_relleno_g", math.nan)
    if m_rel is None or not math.isfinite(m_rel):
        if not motivos:
            motivos.append("sin_llenado")
        return Llenado(res, rr, "ninguna", math.nan, math.nan, motivos)
    m_del = (m_rel - f.get("m_relleno_trasero_g", 0.0)) * G
    ell = rr.ell_delantero  # largo real del tapón delantero (recortado si manda el tope de masa)
    if f["SM_cal"] > SM_max + 1e-9:
        motivos.append("SM_sobre_max")
    return Llenado(res, rr, RESTRICCION.get(f.get("limitante_masa"), "ninguna"), ell, m_del, motivos)
