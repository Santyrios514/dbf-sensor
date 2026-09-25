"""Verificación de que la configuración base coincide con el `.ork` (v2 §5.1 y test 8.1).

Los valores del `.ork` se llevan a un diccionario neutro (`ValoresOrk`), ya sea desde el XML
(sin JVM) o desde OpenRocket (or_bridge), y se comparan con un `Caso` resuelto.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .config import Caso
from .materiales import densidad_equivalente, espesor_total
from .ork_xml import OrkXML

FORMA_OR_A_YAML = {"conical": "conica", "ogive": "ogiva", "ellipsoid": "elipsoide", "power": "potencia",
                   "parabolic": "parabolica", "haack": "haack"}
EQUIV = {"ogiva_tangente": "ogiva"}
MM = 1e-3


@dataclass
class ValoresOrk:
    Ln: float
    R_nariz: float
    t_nariz: float
    forma_nariz: str
    param_nariz: float | None
    Lc: float
    R_cuerpo: float
    t_cuerpo: float
    Lt: float
    Rf: float
    Ra: float
    t_cola: float
    forma_cola: str
    param_cola: float | None
    popa_abierta: bool
    fin_n: int
    fin_t: float
    fin_rot: float
    fin_offset_bottom: float
    fin_puntos: list[tuple[float, float]]
    rho: dict[str, float] = field(default_factory=dict)  # por componente
    recortada: bool | None = None  # solo si la forma de la cola es recortable


def valores_desde_xml(ork: OrkXML, nombres: dict[str, str]) -> ValoresOrk:
    n, c, t, f = (ork.por_nombre(nombres[k]) for k in ("nariz", "cuerpo", "cola", "aletas"))
    if f.campos.get("axialoffset_method") != "bottom":
        raise ValueError(f"aletas '{f.nombre}': se esperaba posición 'bottom', hay "
                         f"'{f.campos.get('axialoffset_method')}'")
    # popa cerrada = hombro trasero con tapa (así cierra la base OpenRocket)
    cerrada = t.campos.get("aftshouldercapped", "false") == "true" and \
        float(t.campos.get("aftshoulderlength", "0")) > 0

    def param(comp):
        return comp.num("shapeparameter") if "shapeparameter" in comp.campos else None

    return ValoresOrk(
        Ln=n.num("length"), R_nariz=n.num("aftradius"), t_nariz=n.num("thickness"),
        forma_nariz=FORMA_OR_A_YAML[n.campos["shape"]], param_nariz=param(n),
        Lc=c.num("length"), R_cuerpo=c.num("radius"), t_cuerpo=c.num("thickness"),
        Lt=t.num("length"), Rf=t.num("foreradius"), Ra=t.num("aftradius"), t_cola=t.num("thickness"),
        forma_cola=FORMA_OR_A_YAML[t.campos["shape"]], param_cola=param(t), popa_abierta=not cerrada,
        fin_n=int(f.num("fincount")), fin_t=f.num("thickness"), fin_rot=f.num("rotation"),
        fin_offset_bottom=f.num("axialoffset"), fin_puntos=list(f.puntos),
        recortada=(t.campos["shapeclipped"] == "true") if "shapeclipped" in t.campos else None,
        rho={k: comp.densidad for k, comp in zip(("nariz", "cuerpo", "cola", "aletas"), (n, c, t, f))},
    )


def comparar(v: ValoresOrk, caso: Caso, tol: float, tol_rel_rho: float = 1e-3,
             tol_puntos: float = 0.1 * MM) -> list[str]:
    """Lista de discrepancias legibles entre el `.ork` y el caso (vacía si coinciden)."""
    g, pa = caso.geom, caso.params_aleta
    dif: list[str] = []

    def chk(nombre, a, b, t=tol):
        if a is None and b is None:
            return
        if a is None or b is None or not math.isclose(a, b, abs_tol=t, rel_tol=0):
            dif.append(f"{nombre}: .ork = {a}, config = {b}")

    chk("nariz.longitud", v.Ln, g.nariz.Ln)
    chk("nariz.radio", v.R_nariz, g.R)
    chk("cuerpo.longitud", v.Lc, g.x_cola - g.nariz.Ln)
    chk("cuerpo.radio", v.R_cuerpo, g.R)
    chk("cola.longitud", v.Lt, g.cola.Lt)
    chk("cola.radio_delantero", v.Rf, g.R)
    chk("cola.radio_popa", v.Ra, g.cola.Ra)
    if EQUIV.get(g.nariz.forma, g.nariz.forma) != v.forma_nariz:
        dif.append(f"nariz.forma: .ork = {v.forma_nariz}, config = {g.nariz.forma}")
    if g.cola.forma != v.forma_cola:
        dif.append(f"cola.forma: .ork = {v.forma_cola}, config = {g.cola.forma}")
    if v.param_cola is not None:
        chk("cola.parametro", v.param_cola, g.cola.parametro, 1e-9)
    if v.recortada is not None and v.recortada != g.cola.recortada:
        dif.append(f"cola.recortada: .ork = {v.recortada}, config = {g.cola.recortada}")
    if v.popa_abierta != (g.cola.popa == "abierta"):
        dif.append(f"cola.popa: .ork abierta = {v.popa_abierta}, config = {g.cola.popa}")
    for est, t_or in (("nariz", v.t_nariz), ("cuerpo", v.t_cuerpo), ("cola", v.t_cola)):
        chk(f"pared.{est}.espesor_total", t_or, espesor_total(caso.pared[est]))
        rho_eq = densidad_equivalente(caso.pared[est])
        if not math.isclose(v.rho[est], rho_eq, rel_tol=tol_rel_rho):
            dif.append(f"pared.{est}: ρ .ork = {v.rho[est]}, ρ_eq config = {rho_eq:.2f}")
    if v.fin_n != pa.n:
        dif.append(f"aletas.n: .ork = {v.fin_n}, config = {pa.n}")
    chk("aletas.espesor", v.fin_t, pa.t)
    chk("aletas.rotacion", v.fin_rot, pa.rotacion, 1e-9)
    chk("aletas.offset_bottom", v.fin_offset_bottom, pa.delta_b)
    if not math.isclose(v.rho["aletas"], pa.rho * pa.phi, rel_tol=tol_rel_rho):
        dif.append(f"aletas: ρ .ork = {v.rho['aletas']}, config = {pa.rho * pa.phi}")
    pts = g.aletas.puntos
    if len(pts) != len(v.fin_puntos):
        dif.append(f"aletas: {len(v.fin_puntos)} puntos en el .ork, {len(pts)} en la config")
    else:
        for i, ((xo, yo), (xc, yc)) in enumerate(zip(v.fin_puntos, pts)):
            if abs(xo - xc) > tol_puntos or abs(yo - yc) > tol_puntos:
                dif.append(f"aletas.P{i}: .ork = ({xo / MM:.3f}, {yo / MM:.3f}) mm, "
                           f"config = ({xc / MM:.3f}, {yc / MM:.3f}) mm")
    return dif
