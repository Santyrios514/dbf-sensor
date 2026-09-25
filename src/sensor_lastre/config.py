"""Carga, validación y expansión del barrido de `config/config.yaml` (spec v1 + v2).

Todo lo que sale de aquí está en SI (m, kg, rad). Los mm/g/grados del YAML se
convierten en este módulo y en ningún otro.
"""

from __future__ import annotations

import copy
import itertools
import math
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import yaml

from . import aletas as mod_aletas
from .atmosfera import isa
from .materiales import Capa, Relleno, costo_granular, densidad_granular, espesor_total
from .perfiles import FORMAS, radio_superficie

MM = 1e-3
G = 1e-3

ESTACIONES = ("nariz", "cuerpo", "cola")
FORMAS_NARIZ = ("conica", "elipsoide", "ogiva_tangente", "ogiva", "potencia", "parabolica", "haack")
FORMAS_COLA = ("conica", "ogiva", "elipsoide", "potencia", "parabolica", "haack")
CON_PARAMETRO = ("potencia",)  # formas cuyo parámetro no tiene valor por defecto


class ConfigError(ValueError):
    """Error de configuración con mensajes legibles."""

    def __init__(self, errores: list[str]):
        self.errores = errores
        super().__init__("Configuración inválida:\n  - " + "\n  - ".join(errores))


# --------------------------------------------------------------------------- dataclasses


@dataclass(frozen=True)
class Nariz:
    forma: str
    parametro: float | None
    Ln: float


@dataclass(frozen=True)
class Cola:
    forma: str
    parametro: float | None
    Lt: float
    Ra: float
    popa: str  # abierta | cerrada
    recortada: bool = True  # clipped de OpenRocket; solo afecta a elipsoide, potencia y haack

    @property
    def popa_cerrada(self) -> bool:
        return self.popa == "cerrada"


@dataclass(frozen=True)
class Geometria:
    L: float  # longitud del cuerpo (nariz + cuerpo + cola)
    D: float
    nariz: Nariz
    cola: Cola
    aletas: "mod_aletas.GeomAleta | None" = None

    @property
    def R(self) -> float:
        return self.D / 2.0

    @property
    def x_cuerpo(self) -> float:
        return self.nariz.Ln

    @property
    def x_cola(self) -> float:
        return self.L - self.cola.Lt

    @property
    def L_total(self) -> float:
        """Longitud total con las aletas que sobresalen de la popa."""
        return max(self.L, self.aletas.x_TE) if self.aletas is not None else self.L


@dataclass(frozen=True)
class Mamparo:
    x: float
    e: float
    material: str
    rho: float
    nombre: str = "mamparo"


@dataclass(frozen=True)
class Electronica:
    Le: float
    me: float
    modo: str  # detras_del_lastre | fija
    x_inicio: float | None
    holgura: float
    x_cg_rel: float | None

    def x_cg(self, x_e):
        return x_e + (self.x_cg_rel if self.x_cg_rel is not None else self.Le / 2.0)


@dataclass(frozen=True)
class MasaPuntual:
    nombre: str
    m: float
    x: float


@dataclass(frozen=True)
class Lastre:
    x_inicio: float | None  # None = auto
    r_min_util: float
    fll: float
    fmax: float | None
    permitir_en_cola: bool
    margen_cola: float
    zonas_prohibidas: tuple[tuple[float, float], ...]
    trasero: bool = False  # tapón trasero detrás de la electrónica (hasta el final de la cavidad)
    margen_popa: float = 0.0  # distancia libre entre el tapón trasero y el final de la cavidad


@dataclass(frozen=True)
class Estabilidad:
    SM_min: float
    SM_max: float | None
    umbral_margen_bajo: float


@dataclass(frozen=True)
class Vuelo:
    V: float
    rho: float
    a: float
    altitud: float
    aoa: float
    mach: float  # el que se usa (auto = V / a(h))
    mach_regresion_or: float | None

    @property
    def q(self) -> float:
        return 0.5 * self.rho * self.V**2


@dataclass(frozen=True)
class Remolque:
    modo: str  # absoluto | en_CG
    x_T: float | None
    alpha_max: float = math.radians(5.0)  # trim admisible (para la tolerancia del amarre)
    alpha_lineal: float = math.radians(15.0)  # validez de la fórmula lineal de trim


@dataclass(frozen=True)
class Envolvente:
    largo_max: float | None
    alto_max: float | None  # de la bahía, con el sensor guardado (acostado)
    ancho_max: float | None
    rot_guardado: float = math.radians(45.0)  # aletas en diagonal al guardarlo acostado


@dataclass(frozen=True)
class Numerico:
    dx: float
    n_ell: int
    tol: float
    n_superficie_aleta: int
    n_franjas_aleta: int
    tol_superficie_aleta: float


@dataclass(frozen=True)
class Caso:
    """Una configuración completamente resuelta (SI)."""

    id: str
    geom: Geometria
    params_aleta: "mod_aletas.ParamsAleta"
    pared: dict[str, tuple[Capa, ...]]
    mamparos: tuple[Mamparo, ...]
    electronica: Electronica
    puntuales: tuple[MasaPuntual, ...]
    lastre: Lastre
    estabilidad: Estabilidad
    vuelo: Vuelo
    remolque: Remolque
    envolvente: Envolvente
    presupuesto: tuple[float, ...]
    numerico: Numerico
    rellenos: tuple[Relleno, ...]
    barrido: dict[str, Any] = field(default_factory=dict)  # ruta -> valor de este caso

    def estacion_de(self, x: float) -> str:
        if x < self.geom.x_cuerpo:
            return "nariz"
        if x < self.geom.x_cola:
            return "cuerpo"
        return "cola"

    def con_aletas(self, g_aleta: "mod_aletas.GeomAleta") -> "Caso":
        """Copia del caso con otra geometría de aleta (p. ej. la exportada por OpenRocket)."""
        return replace(self, geom=replace(self.geom, aletas=g_aleta))


@dataclass
class Config:
    raw: dict[str, Any]
    ruta: Path | None
    casos: list[Caso]
    rellenos: tuple[Relleno, ...]
    openrocket: dict[str, Any] = field(default_factory=dict)
    verificacion: dict[str, Any] = field(default_factory=dict)
    salida: dict[str, Any] = field(default_factory=dict)
    parametros_barrido: tuple[str, ...] = ()

    def caso(self, config_id: str) -> Caso:
        for c in self.casos:
            if c.id == config_id:
                return c
        raise KeyError(config_id)


# --------------------------------------------------------------------------- utilidades


def deep_merge(base: dict, override: dict) -> dict:
    """Fusión en profundidad: los dict se combinan, lo demás se reemplaza."""
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _longitud(spec: Any, bases: dict[str, float], donde: str, errores: list[str]) -> float:
    """Resuelve {modo, valor} a metros. `bases` mapea el sufijo del modo relativo a su base
    (p. ej. {'D': D, 'L': L, 'Lt': Lt})."""
    if isinstance(spec, (int, float)):
        return float(spec) * MM
    if not isinstance(spec, dict) or "modo" not in spec:
        errores.append(f"{donde}: se esperaba {{modo, valor}}, se obtuvo {spec!r}")
        return math.nan
    modo, valor = spec["modo"], spec.get("valor")
    if valor is None:
        errores.append(f"{donde}: falta 'valor'")
        return math.nan
    if modo == "absoluto_mm":
        return float(valor) * MM
    if modo.startswith("relativo_") and modo[len("relativo_"):] in bases:
        return float(valor) * bases[modo[len("relativo_"):]]
    validos = ["absoluto_mm"] + [f"relativo_{k}" for k in bases]
    errores.append(f"{donde}: modo '{modo}' no válido aquí ({' | '.join(validos)})")
    return math.nan


def _material(nombre: str, materiales: dict[str, float], donde: str, errores: list[str]) -> float:
    if nombre not in materiales:
        errores.append(f"{donde}: material '{nombre}' no existe en 'materiales'")
        return math.nan
    return float(materiales[nombre])


def _capas(lista: Any, materiales: dict, donde: str, errores: list[str]) -> tuple[Capa, ...]:
    if not lista:
        errores.append(f"{donde}: la pared necesita al menos una capa")
        return ()
    capas = []
    for i, c in enumerate(lista):
        d = f"{donde}[{i}]"
        rho = _material(c.get("material", ""), materiales, d, errores)
        t = float(c.get("espesor_mm", 0.0)) * MM
        phi = float(c.get("fraccion_solida", 1.0))
        if not t > 0:
            errores.append(f"{d}: espesor_mm debe ser > 0 (vale {t / MM})")
        if not 0 < phi <= 1:
            errores.append(f"{d}: fraccion_solida debe estar en (0, 1] (vale {phi})")
        capas.append(Capa(material=c.get("material", "?"), rho=rho, t=t, phi=phi))
    return tuple(capas)


def _parametro(forma: str, param: Any, donde: str, errores: list[str]) -> float | None:
    if param is None:
        if forma in CON_PARAMETRO:
            errores.append(f"{donde}: la forma '{forma}' necesita 'parametro'")
        return FORMAS.get(forma, (None, None))[1]
    param = float(param)
    if forma == "potencia" and not param > 0:
        errores.append(f"{donde}: el exponente de 'potencia' debe ser > 0")
    if forma == "parabolica" and not 0 <= param <= 1:
        errores.append(f"{donde}: K de 'parabolica' debe estar en [0, 1]")
    if forma in ("ogiva", "ogiva_tangente") and not 0 <= param <= 1:
        errores.append(f"{donde}: el parámetro de 'ogiva' debe estar en [0, 1]")
    return param


# --------------------------------------------------------------------------- resolución


def _rellenos(raw: dict, errores: list[str]) -> tuple[Relleno, ...]:
    mats = raw.get("materiales", {})
    costos = raw.get("costos_usd_kg") or {}
    out = []
    for i, r in enumerate(raw.get("rellenos", [])):
        d = f"rellenos[{i}] ({r.get('nombre', '?')})"
        mat = r.get("material", "")
        rho = _material(mat, mats, d, errores)
        c = float(r.get("costo_usd_kg", costos.get(mat, math.nan)))
        gran = r.get("granular")
        if gran:
            phi = float(gran.get("empaquetamiento", math.nan))
            if not 0 < phi < 1:
                errores.append(f"{d}: empaquetamiento debe estar en (0, 1) (vale {phi})")
            matriz = gran.get("matriz", "")
            rho_m = _material(matriz, mats, d + ".matriz", errores)
            c_grano = float(gran.get("costo_grano_usd_kg", costos.get(mat, math.nan)))
            c = costo_granular(phi, rho, rho_m, c_grano, float(costos.get(matriz, math.nan)))
            rho = densidad_granular(rho, rho_m, phi)
        if math.isfinite(c) and c < 0:
            errores.append(f"{d}: el costo por kg debe ser ≥ 0")
        out.append(Relleno(nombre=r["nombre"], rho_b=rho, costo_usd_kg=c))
    if not out:
        errores.append("rellenos: la lista está vacía")
    return tuple(out)


def _params_aleta(a: dict, L: float, D: float, Lt: float, mats: dict, errores: list[str]):
    if a.get("tipo", "freeform_cola") != "freeform_cola":
        errores.append(f"aletas.tipo '{a.get('tipo')}' no soportado (solo freeform_cola)")
    b = {"D": D, "L": L, "Lt": Lt}
    c_r = _longitud(a["c_r"], b, "aletas.c_r", errores)
    x_s = _longitud(a["x_tip_le"], b, "aletas.x_tip_le", errores)
    ht = a["h_tip"]
    h = r_tip = None
    if isinstance(ht, dict) and ht.get("modo") == "r_tip_absoluto_mm":
        r_tip = float(ht["valor"]) * MM
    else:
        h = _longitud(ht, {"D": D}, "aletas.h_tip", errores)
    e = _longitud(a["extension"], {"D": D}, "aletas.extension", errores)
    t = float(a["espesor_mm"]) * MM
    if not t > 0:
        errores.append("aletas.espesor_mm debe ser > 0")
    n = int(a["n"])
    if n < 1:
        errores.append("aletas.n debe ser ≥ 1")
    phi = float(a.get("fraccion_solida", 1.0))
    if not 0 < phi <= 1:
        errores.append("aletas.fraccion_solida debe estar en (0, 1]")
    return mod_aletas.ParamsAleta(
        n=n, rotacion=math.radians(float(a.get("rotacion_deg", 0.0))), c_r=c_r, x_s=x_s, h=h,
        r_tip_obj=r_tip, e=e, r_in=float(a.get("r_interior_mm", 0.0)) * MM,
        delta_b=float(a.get("offset_bottom_mm", 0.0)) * MM, eps=float(a.get("epsilon_mm", 0.05)) * MM,
        t=t, material=a["material"], rho=_material(a["material"], mats, "aletas.material", errores),
        phi=phi, descontar_solape_eje=bool(a.get("descontar_solape_eje", False)),
    )


def _resolver_caso(cid: str, L_mm: float, D_mm: float, sec: dict, rellenos, errores_glob,
                   valores_barrido: dict) -> Caso | None:
    errores: list[str] = []
    mats = sec["materiales"]
    D = float(D_mm) * MM
    R = D / 2
    gb = sec["geometria_base"]
    ref = gb.get("longitud_referencia", "cuerpo")
    if ref not in ("cuerpo", "total"):
        errores.append(f"geometria_base.longitud_referencia '{ref}' no válida (cuerpo | total)")
    L_in = float(L_mm) * MM

    # extensión de aletas primero: con referencia 'total', L_cuerpo = L − e
    e_aleta = _longitud(gb["aletas"]["extension"], {"D": D}, "aletas.extension", [])
    L = L_in - e_aleta if ref == "total" else L_in
    if gb.get("cuerpo", {}).get("longitud", "auto") != "auto":
        errores.append("geometria_base.cuerpo.longitud: solo se admite 'auto' (L − L_n − L_t)")

    # --- nariz
    n = gb["nariz"]
    forma_n = n.get("forma")
    if forma_n not in FORMAS_NARIZ:
        errores.append(f"nariz.forma '{forma_n}' no válida {FORMAS_NARIZ}")
    param_n = _parametro(forma_n, n.get("parametro"), "nariz", errores)
    Ln = _longitud(n["longitud"], {"D": D, "L": L}, "nariz.longitud", errores)
    nariz = Nariz(forma=forma_n, parametro=param_n, Ln=Ln)

    # --- cola
    c = gb["cola"]
    forma_c = c.get("forma")
    if forma_c not in FORMAS_COLA:
        errores.append(f"cola.forma '{forma_c}' no válida {FORMAS_COLA}")
    param_c = _parametro(forma_c, c.get("parametro"), "cola", errores)
    Lt = _longitud(c["longitud"], {"D": D, "L": L}, "cola.longitud", errores)
    Ra = _longitud(c["diametro_popa"], {"D": D}, "cola.diametro_popa", errores) / 2
    popa = c.get("popa", "abierta")
    if popa not in ("abierta", "cerrada"):
        errores.append(f"cola.popa '{popa}' no válida (abierta | cerrada)")
    cola = Cola(forma=forma_c, parametro=param_c, Lt=Lt, Ra=Ra, popa=popa,
                recortada=bool(c.get("recortada", True)))

    if not Ln > 0:
        errores.append(f"L_n debe ser > 0 (vale {Ln / MM:.2f} mm)")
    if not Lt > 0:
        errores.append(f"L_t debe ser > 0 (vale {Lt / MM:.2f} mm)")
    if not Ln + Lt < L:
        errores.append(f"L_n + L_t = {(Ln + Lt) / MM:.2f} mm debe ser < L_cuerpo = {L / MM:.2f} mm")
    if not 0 < Ra <= R:
        errores.append(f"R_a = {Ra / MM:.2f} mm debe cumplir 0 < R_a ≤ R = {R / MM:.2f} mm")

    geom = Geometria(L=L, D=D, nariz=nariz, cola=cola)

    # --- pared
    p = sec["pared"]
    por_defecto = _capas(p.get("por_defecto"), mats, "pared.por_defecto", errores)
    pared = {est: (_capas(p[est], mats, f"pared.{est}", errores) if p.get(est) else por_defecto)
             for est in ESTACIONES}
    if pared["cola"] and not espesor_total(pared["cola"]) < Ra:
        errores.append("Σ t_k de la cola debe ser < R_a (popa abierta: anillo de pared en la base)")
    for est in ESTACIONES:
        if pared[est] and not espesor_total(pared[est]) < R:
            errores.append(f"pared.{est}: Σ t_k debe ser < R")

    # --- aletas (se resuelven sobre el perfil analítico; en modo OpenRocket se reemplazan)
    pa = _params_aleta(gb["aletas"], L, D, Lt, mats, errores)
    g_aleta = None
    if not errores:
        g_aleta = mod_aletas.construir(
            pa, geom.x_cola, Lt, lambda x: radio_superficie(x, geom),
            int(sec["numerico"].get("n_superficie_aleta", 200)),
            float(sec["numerico"].get("tol_superficie_aleta_mm", 0.01)) * MM, errores)
        geom = replace(geom, aletas=g_aleta)

    # --- mamparos (+ tapa de popa si es cerrada)
    mamparos = []
    for i, m in enumerate(sec.get("mamparos") or []):
        d = f"mamparos[{i}]"
        e = float(m.get("espesor_mm", 0.0)) * MM
        if not e > 0:
            errores.append(f"{d}: espesor_mm debe ser > 0")
        mamparos.append(Mamparo(x=_longitud(m["x"], {"D": D, "L": L}, d + ".x", errores), e=e,
                                material=m["material"], rho=_material(m["material"], mats, d, errores),
                                nombre=m.get("nombre", f"mamparo_{i}")))
    if cola.popa_cerrada:
        mp = c.get("mamparo_popa") or {}
        e_b = float(mp.get("espesor_mm", 0.0)) * MM
        if not e_b > 0:
            errores.append("cola.mamparo_popa.espesor_mm debe ser > 0 con popa cerrada")
        mat = mp.get("material", "")
        mamparos.append(Mamparo(x=L - e_b, e=e_b, material=mat,
                                rho=_material(mat, mats, "cola.mamparo_popa", errores), nombre="tapa_popa"))

    # --- electrónica
    el = sec["electronica"]
    modo_e = el.get("modo")
    if modo_e not in ("detras_del_lastre", "fija"):
        errores.append(f"electronica.modo '{modo_e}' no válido (detras_del_lastre | fija)")
    x_ini_e = el.get("x_inicio_mm")
    if modo_e == "fija" and x_ini_e is None:
        errores.append("electronica.x_inicio_mm es obligatorio con modo 'fija'")
    xcg_rel = el.get("x_cg_relativo_mm")
    electronica = Electronica(
        Le=float(el["longitud_mm"]) * MM, me=float(el["masa_g"]) * G, modo=modo_e,
        x_inicio=None if x_ini_e is None else float(x_ini_e) * MM,
        holgura=float(el.get("holgura_mm", 0.0)) * MM,
        x_cg_rel=None if xcg_rel is None else float(xcg_rel) * MM,
    )

    puntuales = tuple(
        MasaPuntual(nombre=mp.get("nombre", f"p{i}"), m=float(mp["masa_g"]) * G,
                    x=_longitud(mp["x"], {"D": D, "L": L}, f"masas_puntuales[{i}].x", errores))
        for i, mp in enumerate(sec.get("masas_puntuales") or [])
    )

    # --- lastre
    la = sec["lastre"]
    xi = la.get("x_inicio", "auto")
    x_inicio = None if xi in (None, "auto") else _longitud(xi, {"D": D, "L": L}, "lastre.x_inicio", errores)
    fll = float(la.get("factor_llenado", 1.0))
    if not 0 < fll <= 1:
        errores.append(f"lastre.factor_llenado debe estar en (0, 1] (vale {fll})")
    fmax = la.get("fraccion_max_L")
    zonas = []
    for i, z in enumerate(la.get("zonas_prohibidas") or []):
        z0, z1 = float(z["desde_mm"]) * MM, float(z["hasta_mm"]) * MM
        if not z1 > z0:
            errores.append(f"lastre.zonas_prohibidas[{i}]: hasta_mm debe ser > desde_mm")
        zonas.append((z0, z1))
    lastre = Lastre(
        x_inicio=x_inicio, r_min_util=float(la.get("r_min_util_mm", 0.0)) * MM, fll=fll,
        fmax=None if fmax is None else float(fmax),
        permitir_en_cola=bool(la.get("permitir_en_cola", False)),
        margen_cola=float(la.get("margen_cola_mm", 0.0)) * MM, zonas_prohibidas=tuple(zonas),
        trasero=bool(la.get("lastre_trasero", False)), margen_popa=float(la.get("margen_popa_mm", 0.0)) * MM,
    )
    if lastre.trasero and electronica.modo != "detras_del_lastre":
        errores.append("lastre.lastre_trasero requiere electronica.modo = detras_del_lastre")

    # --- estabilidad
    es = sec["estabilidad"]
    SM_min = float(es["SM_min_cal"])
    SM_max = es.get("SM_max_cal")
    if SM_max is not None and not float(SM_max) > SM_min:
        errores.append(f"SM_max_cal ({SM_max}) debe ser > SM_min_cal ({SM_min})")
    estab = Estabilidad(SM_min=SM_min, SM_max=None if SM_max is None else float(SM_max),
                        umbral_margen_bajo=float(es.get("umbral_margen_bajo_cal", 0.0)))

    # --- vuelo
    cv = sec["condiciones_vuelo"]
    h = float(cv.get("altitud_m", 0.0))
    est_isa = isa(h)
    rho_aire = est_isa.rho if cv.get("rho_aire", "isa") == "isa" else float(cv["rho_aire"])
    V = float(cv["V_m_s"])
    mach = cv.get("mach", "auto")
    mach = V / est_isa.a if mach in (None, "auto") else float(mach)
    mr = cv.get("mach_regresion_or")
    vuelo = Vuelo(V=V, rho=rho_aire, a=est_isa.a, altitud=h,
                  aoa=math.radians(float(cv.get("aoa_deg", 0.0))), mach=mach,
                  mach_regresion_or=None if mr is None else float(mr))

    # --- remolque
    rq = sec["remolque"]
    xT = rq["x_T"]
    a_max = math.radians(float(rq.get("alpha_trim_max_deg", 5.0)))
    a_lin = math.radians(float(rq.get("alpha_lineal_max_deg", 15.0)))
    if not 0 < a_max <= a_lin:
        errores.append("remolque: se requiere 0 < alpha_trim_max_deg ≤ alpha_lineal_max_deg")
    if xT == "en_CG" or isinstance(xT, dict) and xT.get("modo") == "en_CG":
        remolque = Remolque(modo="en_CG", x_T=None, alpha_max=a_max, alpha_lineal=a_lin)
    else:
        remolque = Remolque(modo="absoluto", x_T=_longitud(xT, {"D": D, "L": L}, "remolque.x_T", errores),
                            alpha_max=a_max, alpha_lineal=a_lin)

    env = sec["envolvente"]

    def _mm_o_none(k):
        v = env.get(k)
        return None if v is None else float(v) * MM

    envolvente = Envolvente(largo_max=_mm_o_none("largo_max_mm"), alto_max=_mm_o_none("alto_max_mm"),
                            ancho_max=_mm_o_none("ancho_max_mm"),
                            rot_guardado=math.radians(float(env.get("rotacion_guardado_deg", 45.0))))
    nu = sec["numerico"]
    numerico = Numerico(
        dx=float(nu["dx_mm"]) * MM, n_ell=int(nu["n_barrido_ell"]), tol=float(nu["tol_raiz_mm"]) * MM,
        n_superficie_aleta=int(nu.get("n_superficie_aleta", 200)),
        n_franjas_aleta=int(nu.get("n_franjas_aleta", 48)),
        tol_superficie_aleta=float(nu.get("tol_superficie_aleta_mm", 0.01)) * MM,
    )
    if not numerico.dx > 0:
        errores.append("numerico.dx_mm debe ser > 0")

    presupuesto = tuple(float(m) * G for m in (sec.get("presupuesto_masa") or {}).get("lastre_g", []))

    if errores:
        errores_glob.extend(f"[{cid}] {e}" for e in errores)
        return None
    return Caso(
        id=cid, geom=geom, params_aleta=pa, pared=pared, mamparos=tuple(mamparos),
        electronica=electronica, puntuales=puntuales, lastre=lastre, estabilidad=estab, vuelo=vuelo,
        remolque=remolque, envolvente=envolvente, presupuesto=presupuesto, numerico=numerico,
        rellenos=rellenos, barrido=dict(valores_barrido),
    )


# Secciones que un override de `configuraciones` o una ruta del barrido pueden modificar.
SECCIONES_CASO = (
    "pared", "mamparos", "electronica", "masas_puntuales", "lastre", "estabilidad",
    "condiciones_vuelo", "remolque", "envolvente", "presupuesto_masa", "numerico",
)
ESPECIALES = ("L_mm", "D_mm")


def _fmt_num(v) -> str:
    return f"{v:g}" if isinstance(v, (int, float)) else str(v)


def _base_secciones(raw: dict) -> dict:
    base = {k: copy.deepcopy(raw.get(k)) for k in SECCIONES_CASO}
    base["geometria_base"] = copy.deepcopy(raw["geometria_base"])
    base["materiales"] = raw["materiales"]
    return base


def _ubicar(sec: dict, ruta: str) -> tuple[dict, str] | None:
    """(dict contenedor, clave final) de una ruta con puntos. Se busca primero en las secciones
    del caso y luego dentro de geometria_base. None si la ruta no existe."""
    partes = ruta.split(".")
    for raiz in (sec, sec["geometria_base"]):
        d = raiz
        ok = True
        for p in partes[:-1]:
            if isinstance(d, dict) and isinstance(d.get(p), dict):
                d = d[p]
            else:
                ok = False
                break
        if ok and isinstance(d, dict) and partes[-1] in d:
            return d, partes[-1]
    return None


def _etiqueta(ruta: str, rutas: list[str]) -> str:
    partes = [p for p in ruta.split(".") if p != "valor"]
    corta = partes[-1]
    if sum(1 for r in rutas if [p for p in r.split(".") if p != "valor"][-1] == corta) > 1:
        corta = "_".join(partes)
    return corta


def expandir(raw: dict) -> tuple[list[tuple[str, float, float, dict, dict]], list[str]]:
    """Barrido (cartesiano o zip) + overrides explícitos → [(id, L_mm, D_mm, secciones, valores)]."""
    errores: list[str] = []
    b = raw.get("barrido") or {}
    modo = b.get("modo", "cartesiano")
    params = dict(b.get("parametros") or {})
    if not params:
        errores.append("barrido.parametros está vacío")
        return [], errores
    for k in ESPECIALES:
        if k not in params:
            errores.append(f"barrido.parametros debe incluir '{k}'")
    rutas = [k for k in params if k not in ESPECIALES]
    base = _base_secciones(raw)
    for r in rutas:
        if _ubicar(base, r) is None:
            errores.append(f"barrido: la ruta '{r}' no existe en el YAML")
    if errores:
        return [], errores
    claves = list(params)
    listas = [list(v) if isinstance(v, (list, tuple)) else [v] for v in params.values()]
    if modo == "cartesiano":
        combos = list(itertools.product(*listas))
    elif modo == "zip":
        if len({len(x) for x in listas}) != 1:
            errores.append("barrido.modo 'zip' requiere listas de la misma longitud")
            return [], errores
        combos = list(zip(*listas))
    else:
        errores.append(f"barrido.modo '{modo}' no válido (cartesiano | zip)")
        return [], errores

    out = []
    for combo in combos:
        valores = dict(zip(claves, combo))
        sec = copy.deepcopy(base)
        for r in rutas:
            d, k = _ubicar(sec, r)
            d[k] = valores[r]
        cid = f"L{_fmt_num(valores['L_mm'])}_D{_fmt_num(valores['D_mm'])}"
        cid += "".join(f"_{_etiqueta(r, rutas)}{_fmt_num(valores[r])}" for r in rutas)
        out.append((cid, valores["L_mm"], valores["D_mm"], sec, {r: valores[r] for r in rutas}))

    for cid, ov in (raw.get("configuraciones") or {}).items():
        ov = dict(ov or {})
        L_mm, D_mm = ov.pop("L_mm"), ov.pop("D_mm")
        sec = copy.deepcopy(base)
        geo_ov = {k: v for k, v in ov.items() if k in raw["geometria_base"]}
        sec["geometria_base"] = deep_merge(sec["geometria_base"], geo_ov)
        for k, v in ov.items():
            if k in SECCIONES_CASO:
                sec[k] = deep_merge(sec[k], v) if isinstance(v, dict) and isinstance(sec[k], dict) else v
            elif k not in geo_ov:
                errores.append(f"configuraciones.{cid}: clave '{k}' desconocida")
        out.append((str(cid), L_mm, D_mm, sec, {}))
    return out, errores


REQUERIDAS = ("materiales", "rellenos", "geometria_base", "pared", "electronica", "lastre",
              "estabilidad", "condiciones_vuelo", "remolque", "envolvente", "numerico", "barrido")


def cargar(ruta: str | Path | dict) -> Config:
    """Carga y valida el YAML (o un dict ya cargado). Lanza ConfigError si algo falla."""
    if isinstance(ruta, dict):
        raw, path = copy.deepcopy(ruta), None
    else:
        path = Path(ruta)
        with open(path, encoding="utf-8") as f:
            raw = yaml.safe_load(f)
    faltan = [k for k in REQUERIDAS if k not in raw]
    if faltan:
        raise ConfigError([f"falta la sección '{k}'" for k in faltan])
    errores: list[str] = []
    rellenos = _rellenos(raw, errores)
    expandidos, err_exp = expandir(raw)
    errores += err_exp
    casos = []
    for cid, L, D, sec, vals in expandidos:
        c = _resolver_caso(cid, L, D, sec, rellenos, errores, vals)
        if c is not None:
            casos.append(c)
    ids = [c.id for c in casos]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        errores.append(f"ids de configuración repetidos: {sorted(dup)}")
    if not expandidos and not err_exp:
        errores.append("no hay configuraciones (barrido y configuraciones vacíos)")
    if errores:
        raise ConfigError(errores)
    rutas = tuple(k for k in (raw["barrido"].get("parametros") or {}) if k not in ESPECIALES)
    return Config(
        raw=raw, ruta=path, casos=casos, rellenos=rellenos,
        openrocket=raw.get("openrocket") or {},
        verificacion=raw.get("verificacion") or {},
        salida=raw.get("salida") or {},
        parametros_barrido=rutas,
    )

