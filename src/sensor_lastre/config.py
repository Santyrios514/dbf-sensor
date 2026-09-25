"""Carga, validación y expansión del barrido de `config/config.yaml`.

Todo lo que sale de aquí está en SI (m, kg, rad). Los mm/g/grados del YAML se
convierten en este módulo y en ningún otro.
"""

from __future__ import annotations

import copy
import itertools
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .atmosfera import isa
from .materiales import Capa, Relleno, densidad_granular, espesor_total

MM = 1e-3
G = 1e-3

ESTACIONES = ("nariz", "cuerpo", "cola")
FORMAS_NARIZ = ("conica", "elipsoide", "ogiva_tangente", "potencia", "haack")
FORMAS_COLA = ("conica", "ogiva", "parabolica")


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
    Lt: float
    Ra: float
    popa_cerrada: bool


@dataclass(frozen=True)
class Aletas:
    n: int
    cr: float
    ct: float
    s: float
    xs: float
    tf: float
    material: str
    rho: float
    phi: float
    x_r0: float


@dataclass(frozen=True)
class Geometria:
    L: float
    D: float
    nariz: Nariz
    cola: Cola
    aletas: Aletas

    @property
    def R(self) -> float:
        return self.D / 2.0

    @property
    def x_cuerpo(self) -> float:
        return self.nariz.Ln

    @property
    def x_cola(self) -> float:
        return self.L - self.cola.Lt


@dataclass(frozen=True)
class Mamparo:
    x: float
    e: float
    material: str
    rho: float


@dataclass(frozen=True)
class Electronica:
    Le: float
    me: float
    modo: str  # detras_del_lastre | fija
    x_inicio: float | None
    holgura: float
    x_cg_rel: float | None

    def x_cg(self, x_e: float) -> float:
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


@dataclass(frozen=True)
class Estabilidad:
    SM_min: float
    SM_max: float | None


@dataclass(frozen=True)
class Vuelo:
    V: float
    rho: float
    a: float
    altitud: float
    aoa: float

    @property
    def q(self) -> float:
        return 0.5 * self.rho * self.V**2

    @property
    def mach(self) -> float:
        return self.V / self.a


@dataclass(frozen=True)
class Remolque:
    modo: str  # absoluto | en_CG
    x_T: float | None


@dataclass(frozen=True)
class Envolvente:
    largo_max: float
    diametro_max: float


@dataclass(frozen=True)
class Numerico:
    dx: float
    n_ell: int
    tol: float


@dataclass(frozen=True)
class Caso:
    """Una configuración geométrica completamente resuelta (SI)."""

    id: str
    geom: Geometria
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

    def estacion_de(self, x: float) -> str:
        if x < self.geom.x_cuerpo:
            return "nariz"
        if x < self.geom.x_cola:
            return "cuerpo"
        return "cola"


@dataclass
class Config:
    raw: dict[str, Any]
    ruta: Path | None
    casos: list[Caso]
    rellenos: tuple[Relleno, ...]
    openrocket: dict[str, Any] = field(default_factory=dict)
    verificacion: dict[str, float] = field(default_factory=dict)
    salida: dict[str, Any] = field(default_factory=dict)

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


def _longitud(spec: Any, L: float, D: float, donde: str, errores: list[str]) -> float:
    """Resuelve {modo: absoluto_mm | relativo_D | relativo_L, valor} a metros."""
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
    if modo == "relativo_D":
        return float(valor) * D
    if modo == "relativo_L":
        return float(valor) * L
    errores.append(f"{donde}: modo '{modo}' desconocido (absoluto_mm | relativo_D | relativo_L)")
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


# --------------------------------------------------------------------------- resolución


def _rellenos(raw: dict, errores: list[str]) -> tuple[Relleno, ...]:
    mats = raw.get("materiales", {})
    out = []
    for i, r in enumerate(raw.get("rellenos", [])):
        d = f"rellenos[{i}] ({r.get('nombre', '?')})"
        rho = _material(r.get("material", ""), mats, d, errores)
        gran = r.get("granular")
        if gran:
            phi = float(gran.get("empaquetamiento", math.nan))
            if not 0 < phi < 1:
                errores.append(f"{d}: empaquetamiento debe estar en (0, 1) (vale {phi})")
            rho_m = _material(gran.get("matriz", ""), mats, d + ".matriz", errores)
            rho = densidad_granular(rho, rho_m, phi)
        out.append(Relleno(nombre=r["nombre"], rho_b=rho))
    if not out:
        errores.append("rellenos: la lista está vacía")
    return tuple(out)


def _resolver_caso(cid: str, L_mm: float, D_mm: float, sec: dict, rellenos, errores_glob) -> Caso:
    errores: list[str] = []
    mats = sec["materiales"]
    L, D = float(L_mm) * MM, float(D_mm) * MM
    R = D / 2
    gb = sec["geometria_base"]

    # --- nariz
    n = gb["nariz"]
    forma_n = n.get("forma")
    if forma_n not in FORMAS_NARIZ:
        errores.append(f"nariz.forma '{forma_n}' no válida {FORMAS_NARIZ}")
    param = n.get("parametro")
    if forma_n == "potencia" and (param is None or not param > 0):
        errores.append("nariz.parametro (n) debe ser > 0 para forma 'potencia'")
    if forma_n == "haack" and param is None:
        param = 0.0
    Ln = _longitud(n["longitud"], L, D, "nariz.longitud", errores)
    nariz = Nariz(forma=forma_n, parametro=None if param is None else float(param), Ln=Ln)

    # --- cola
    c = gb["cola"]
    if c.get("forma") not in FORMAS_COLA:
        errores.append(f"cola.forma '{c.get('forma')}' no válida {FORMAS_COLA}")
    Lt = _longitud(c["longitud"], L, D, "cola.longitud", errores)
    Ra = _longitud(c["diametro_popa"], L, D, "cola.diametro_popa", errores) / 2
    cola = Cola(forma=c.get("forma"), Lt=Lt, Ra=Ra, popa_cerrada=bool(c.get("popa_cerrada", True)))

    if not Ln > 0:
        errores.append(f"L_n debe ser > 0 (vale {Ln / MM:.2f} mm)")
    if not Lt >= 0:
        errores.append(f"L_t debe ser ≥ 0 (vale {Lt / MM:.2f} mm)")
    if not Ln + Lt < L:
        errores.append(f"L_n + L_t = {(Ln + Lt) / MM:.2f} mm debe ser < L = {L_mm} mm")
    if not Ra <= R:
        errores.append(f"R_a = {Ra / MM:.2f} mm debe ser ≤ R = {R / MM:.2f} mm")

    # --- aletas
    a = gb["aletas"]
    cr = _longitud(a["cuerda_raiz"], L, D, "aletas.cuerda_raiz", errores)
    ct = _longitud(a["cuerda_punta"], L, D, "aletas.cuerda_punta", errores)
    s = _longitud(a["semienvergadura"], L, D, "aletas.semienvergadura", errores)
    fl = a.get("flecha", {"modo": "borde_salida_recto"})
    if fl.get("modo") == "borde_salida_recto":
        xs = cr - ct
    else:
        xs = _longitud(fl, L, D, "aletas.flecha", errores)
    tf = float(a["espesor_mm"]) * MM
    if not tf > 0:
        errores.append("aletas.espesor_mm debe ser > 0")
    rho_f = _material(a["material"], mats, "aletas.material", errores)
    phi_f = float(a.get("fraccion_solida", 1.0))
    x_r0 = L - float(a.get("borde_salida_desde_popa_mm", 0.0)) * MM - cr
    aletas = Aletas(
        n=int(a["n"]), cr=cr, ct=ct, s=s, xs=xs, tf=tf, material=a["material"],
        rho=rho_f, phi=phi_f, x_r0=x_r0,
    )
    geom = Geometria(L=L, D=D, nariz=nariz, cola=cola, aletas=aletas)

    # --- pared
    p = sec["pared"]
    por_defecto = _capas(p.get("por_defecto"), mats, "pared.por_defecto", errores)
    pared = {}
    for est in ESTACIONES:
        pared[est] = _capas(p[est], mats, f"pared.{est}", errores) if p.get(est) else por_defecto
    if not cola.popa_cerrada and pared["cola"] and not espesor_total(pared["cola"]) < Ra:
        errores.append("con popa abierta, Σ t_k de la cola debe ser < R_a")
    for est in ESTACIONES:
        if pared[est] and not espesor_total(pared[est]) < R:
            errores.append(f"pared.{est}: Σ t_k debe ser < R")

    # --- mamparos
    mamparos = []
    for i, m in enumerate(sec.get("mamparos") or []):
        d = f"mamparos[{i}]"
        e = float(m.get("espesor_mm", 0.0)) * MM
        if not e > 0:
            errores.append(f"{d}: espesor_mm debe ser > 0")
        mamparos.append(Mamparo(
            x=_longitud(m["x"], L, D, d + ".x", errores), e=e, material=m["material"],
            rho=_material(m["material"], mats, d, errores),
        ))

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
                    x=_longitud(mp["x"], L, D, f"masas_puntuales[{i}].x", errores))
        for i, mp in enumerate(sec.get("masas_puntuales") or [])
    )

    # --- lastre
    la = sec["lastre"]
    xi = la.get("x_inicio", "auto")
    x_inicio = None if xi in (None, "auto") else _longitud(xi, L, D, "lastre.x_inicio", errores)
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
    )

    # --- estabilidad
    es = sec["estabilidad"]
    SM_min = float(es["SM_min_cal"])
    SM_max = es.get("SM_max_cal")
    if SM_max is not None and not float(SM_max) > SM_min:
        errores.append(f"SM_max_cal ({SM_max}) debe ser > SM_min_cal ({SM_min})")
    estab = Estabilidad(SM_min=SM_min, SM_max=None if SM_max is None else float(SM_max))

    # --- vuelo
    cv = sec["condiciones_vuelo"]
    h = float(cv.get("altitud_m", 0.0))
    est_isa = isa(h)
    rho_aire = est_isa.rho if cv.get("rho_aire", "isa") == "isa" else float(cv["rho_aire"])
    vuelo = Vuelo(V=float(cv["V_m_s"]), rho=rho_aire, a=est_isa.a, altitud=h,
                  aoa=math.radians(float(cv.get("aoa_deg", 0.0))))

    # --- remolque
    xT = sec["remolque"]["x_T"]
    if isinstance(xT, str) and xT == "en_CG" or isinstance(xT, dict) and xT.get("modo") == "en_CG":
        remolque = Remolque(modo="en_CG", x_T=None)
    else:
        remolque = Remolque(modo="absoluto", x_T=_longitud(xT, L, D, "remolque.x_T", errores))

    env = sec["envolvente"]
    envolvente = Envolvente(largo_max=float(env["largo_max_mm"]) * MM,
                            diametro_max=float(env["diametro_max_mm"]) * MM)
    nu = sec["numerico"]
    numerico = Numerico(dx=float(nu["dx_mm"]) * MM, n_ell=int(nu["n_barrido_ell"]),
                        tol=float(nu["tol_raiz_mm"]) * MM)
    if not numerico.dx > 0:
        errores.append("numerico.dx_mm debe ser > 0")

    presupuesto = tuple(float(m) * G for m in (sec.get("presupuesto_masa") or {}).get("lastre_g", []))

    if errores:
        errores_glob.extend(f"[{cid}] {e}" for e in errores)
    return Caso(
        id=cid, geom=geom, pared=pared, mamparos=tuple(mamparos), electronica=electronica,
        puntuales=puntuales, lastre=lastre, estabilidad=estab, vuelo=vuelo, remolque=remolque,
        envolvente=envolvente, presupuesto=presupuesto, numerico=numerico, rellenos=rellenos,
    )


# Secciones que un override de `configuraciones` puede modificar además de la geometría.
SECCIONES_CASO = (
    "pared", "mamparos", "electronica", "masas_puntuales", "lastre", "estabilidad",
    "condiciones_vuelo", "remolque", "envolvente", "presupuesto_masa", "numerico",
)


def _fmt_num(v: float) -> str:
    return f"{v:g}"


def expandir(raw: dict) -> list[tuple[str, float, float, dict]]:
    """Barrido cartesiano L×D + overrides explícitos → [(id, L_mm, D_mm, secciones)]."""
    base = {k: raw.get(k) for k in SECCIONES_CASO}
    base["geometria_base"] = raw["geometria_base"]
    base["materiales"] = raw["materiales"]
    out = []
    b = raw.get("barrido") or {}
    for L_mm, D_mm in itertools.product(b.get("L_mm", []), b.get("D_mm", [])):
        out.append((f"L{_fmt_num(L_mm)}_D{_fmt_num(D_mm)}", L_mm, D_mm, base))
    for cid, ov in (raw.get("configuraciones") or {}).items():
        ov = dict(ov or {})
        L_mm, D_mm = ov.pop("L_mm"), ov.pop("D_mm")
        sec = dict(base)
        geo_ov = {k: v for k, v in ov.items() if k in raw["geometria_base"]}
        sec["geometria_base"] = deep_merge(raw["geometria_base"], geo_ov)
        for k, v in ov.items():
            if k in SECCIONES_CASO:
                sec[k] = deep_merge(base[k], v) if isinstance(v, dict) and isinstance(base[k], dict) else v
            elif k not in geo_ov:
                raise ConfigError([f"configuraciones.{cid}: clave '{k}' desconocida"])
        out.append((str(cid), L_mm, D_mm, sec))
    return out


REQUERIDAS = ("materiales", "rellenos", "geometria_base", "pared", "electronica", "lastre",
              "estabilidad", "condiciones_vuelo", "remolque", "envolvente", "numerico")


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
    casos = [_resolver_caso(cid, L, D, sec, rellenos, errores) for cid, L, D, sec in expandir(raw)]
    ids = [c.id for c in casos]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        errores.append(f"ids de configuración repetidos: {sorted(dup)}")
    if not casos:
        errores.append("no hay configuraciones (barrido y configuraciones vacíos)")
    if errores:
        raise ConfigError(errores)
    return Config(
        raw=raw, ruta=path, casos=casos, rellenos=rellenos,
        openrocket=raw.get("openrocket") or {},
        verificacion=raw.get("verificacion") or {},
        salida=raw.get("salida") or {},
    )
