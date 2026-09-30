"""Carga de `config/optimizacion.yaml`, herencia de `base_config` y expansión de la malla (v3 §2, §5).

La malla se separa en dos niveles, porque todo lo que no depende de las aletas (perfil, erosión,
volúmenes, masas del casco, CP del cuerpo) se calcula una vez por cuerpo (v3 §4.2):

- cuerpo: (D, forma y parámetro de la cola, L_t/D, k);
- aleta:  (λ_LE, μ, r_tip/R, γ, σ).

Cada cuerpo se traduce a un YAML de `sensor_lastre` (`raw_cuerpo`) y se resuelve con su
`cargar`, de modo que materiales, pared, electrónica, lastre, remolque y vuelo salen del código
existente sin reinterpretarlos.
"""

from __future__ import annotations

import copy
import itertools
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from sensor_lastre.config import FORMAS_COLA, ConfigError, deep_merge

MM = 1e-3
G = 1e-3
RELLENO = "plomo_macizo"  # v3 §1: único relleno

VARS_CUERPO = ("D_mm", "L_t_rel_D", "k")
VARS_ALETA = ("lambda_LE", "mu_cr", "r_tip_rel_R", "gamma_ct", "sigma_flecha")


def _num(v: float) -> str:
    return f"{v:g}"


@dataclass(frozen=True)
class CuerpoSpec:
    D_mm: float
    forma: str
    parametro: float | None
    L_t_rel_D: float
    k: float

    @property
    def id(self) -> str:
        p = "" if self.parametro is None else _num(self.parametro)
        return f"D{_num(self.D_mm)}_{self.forma}{p}_Lt{_num(self.L_t_rel_D)}_k{_num(self.k)}"

    @property
    def D(self) -> float:
        return self.D_mm * MM

    @property
    def L_t(self) -> float:
        return self.L_t_rel_D * self.D


@dataclass(frozen=True)
class AletaSpec:
    lambda_LE: float
    mu_cr: float
    r_tip_rel_R: float
    gamma_ct: float
    sigma_flecha: float

    @property
    def id(self) -> str:
        return (f"l{_num(self.lambda_LE)}_m{_num(self.mu_cr)}_r{_num(self.r_tip_rel_R)}"
                f"_g{_num(self.gamma_ct)}_s{_num(self.sigma_flecha)}")


def cand_id(c: CuerpoSpec, a: AletaSpec) -> str:
    return f"{c.id}_{a.id}"


@dataclass(frozen=True)
class ParamsAletaOpt:
    n: int
    rotacion: float  # rad
    material: str
    rho: float  # kg/m³ (sin fracción sólida)
    t: float
    phi: float
    E: float  # Pa
    nu: float
    factor_flutter: float


@dataclass(frozen=True)
class Restricciones:
    SM_min: float
    SM_max: float
    k_min: float
    h_min: float
    c_min: float
    eps_CN: float
    angulo_cola_max: float | None  # rad


@dataclass(frozen=True)
class Objetivo:
    pesos: tuple[float, float, float, float]
    f4_modo: str
    SM_centro: float


@dataclass
class ConfigOpt:
    raw: dict[str, Any]
    ruta: Path | None
    base_raw: dict[str, Any]  # config de sensor_lastre ya fusionada con los materiales propios
    m_max: float  # kg
    factor_llenado: float
    L: float
    nariz: dict[str, Any]
    aleta: ParamsAletaOpt
    malla: dict[str, list]
    restricciones: Restricciones
    objetivo: Objetivo
    ejecucion: dict[str, Any]
    numerico: dict[str, Any]
    salida: dict[str, Any]
    raiz: Path = field(default_factory=Path.cwd)
    base_ruta: Path | None = None  # ruta de base_config; sus rutas relativas (p. ej. el .ork) cuelgan de su raíz

    def ruta_base(self, v: str) -> Path:
        """Ruta relativa de base_config (p. ej. openrocket.ork_base) resuelta desde su proyecto."""
        p = Path(v)
        if p.is_absolute():
            return p
        raiz = self.base_ruta.resolve().parents[1] if self.base_ruta is not None else self.raiz
        return raiz / p

    # ------------------------------------------------------------------ malla
    def cuerpos(self, malla: dict | None = None) -> list[CuerpoSpec]:
        m = malla or self.malla
        return [CuerpoSpec(float(D), f["forma"], None if f.get("parametro") is None else float(f["parametro"]),
                           float(lt), float(k))
                for D, f, lt, k in itertools.product(m["D_mm"], m["cola_forma"], m["L_t_rel_D"], m["k"])]

    def aletas(self, malla: dict | None = None) -> list[AletaSpec]:
        m = malla or self.malla
        return [AletaSpec(*(float(v) for v in combo))
                for combo in itertools.product(*(m[k] for k in VARS_ALETA))]

    @property
    def n_evaluaciones(self) -> int:
        return len(self.cuerpos()) * len(self.aletas())

    @property
    def procesos(self) -> int:
        p = self.ejecucion.get("procesos", "auto")
        return max(1, os.cpu_count() or 1) if p in (None, "auto") else max(1, int(p))

    # ------------------------------------------------------------------ cotas a priori (v3 §3.8)
    @property
    def cotas(self) -> dict[str, float]:
        """Cotas de normalización tomadas de la malla configurada, no de los resultados."""
        D = [float(d) * MM for d in self.malla["D_mm"]]
        k = [float(v) for v in self.malla["k"]]
        r_tip_max = max(float(v) for v in self.malla["r_tip_rel_R"]) * max(D) / 2
        return {"D_ap_lo": min(D), "D_ap_hi": 2 * r_tip_max, "k_lo": min(k), "k_hi": max(k)}

    def dir_salida(self) -> Path:
        return self._ruta(self.salida.get("dir", "data_opt"))

    def dir_figuras(self) -> Path:
        return self._ruta(self.salida.get("dir_figuras", "figs_opt"))

    def _ruta(self, v: str) -> Path:
        p = Path(v)
        return p if p.is_absolute() else self.raiz / p


# --------------------------------------------------------------------------- carga


def _resolver(ruta: str, *bases: Path) -> Path:
    p = Path(ruta)
    if p.is_absolute():
        return p
    for b in bases:
        if (b / p).exists():
            return b / p
    return bases[0] / p


def cargar(ruta: str | Path | dict, raiz: Path | None = None) -> ConfigOpt:
    """Carga y valida optimizacion.yaml (o un dict). Lanza ConfigError con todos los problemas."""
    if isinstance(ruta, dict):
        raw, path = copy.deepcopy(ruta), None
        raiz = raiz or Path.cwd()
    else:
        path = Path(ruta).resolve()
        with open(path, encoding="utf-8") as f:
            raw = yaml.safe_load(f)
        raiz = raiz or path.parents[1]
    errores: list[str] = []
    for sec in ("base_config", "masa", "geometria_fija", "malla", "restricciones", "objetivo"):
        if sec not in raw:
            errores.append(f"falta la sección '{sec}'")
    if errores:
        raise ConfigError(errores)

    base_path = _resolver(raw["base_config"], raiz, Path.cwd())
    with open(base_path, encoding="utf-8") as f:
        base = yaml.safe_load(f)
    base["materiales"] = {**base.get("materiales", {}), **(raw.get("materiales") or {})}
    mats = base["materiales"]

    ms = raw["masa"]
    m_max = float(ms["m_max_g"]) * G
    if not m_max > 0:
        errores.append("masa.m_max_g debe ser > 0")
    fll = float(ms.get("factor_llenado", 1.0))
    if not 0 < fll <= 1:
        errores.append("masa.factor_llenado debe estar en (0, 1]")

    gf = raw["geometria_fija"]
    L = float(gf["L_total_mm"]) * MM
    a = gf["aletas"]
    if a.get("material") not in mats:
        errores.append(f"geometria_fija.aletas.material '{a.get('material')}' no está en materiales")
    aleta = ParamsAletaOpt(
        n=int(a.get("n", 4)), rotacion=math.radians(float(a.get("rotacion_deg", 0.0))),
        material=a.get("material", "?"), rho=float(mats.get(a.get("material"), math.nan)),
        t=float(a["espesor_mm"]) * MM, phi=float(a.get("fraccion_solida", 1.0)),
        E=float(a.get("E_GPa", math.nan)) * 1e9, nu=float(a.get("nu", math.nan)),
        factor_flutter=float(a.get("factor_flutter", 3.0)))
    if not aleta.t > 0:
        errores.append("geometria_fija.aletas.espesor_mm debe ser > 0")
    if not 0 < aleta.phi <= 1:
        errores.append("geometria_fija.aletas.fraccion_solida debe estar en (0, 1]")

    m = raw["malla"]
    for k in ("D_mm", "cola_forma") + VARS_CUERPO[1:] + VARS_ALETA:
        if not m.get(k):
            errores.append(f"malla.{k} está vacía o falta")
    r = raw["restricciones"]
    rest = Restricciones(
        SM_min=float(r["SM_min_cal"]), SM_max=float(r["SM_max_cal"]), k_min=float(r["k_min"]),
        h_min=float(r.get("h_min_mm", 5.0)) * MM, c_min=float(r.get("c_min_mm", 5.0)) * MM,
        eps_CN=float(r.get("eps_CN", 0.5)),
        angulo_cola_max=None if r.get("angulo_cola_max_deg") is None else math.radians(float(r["angulo_cola_max_deg"])))
    if not rest.SM_max > rest.SM_min:
        errores.append("restricciones: SM_max_cal debe ser > SM_min_cal")
    if not rest.k_min > 0:
        errores.append("restricciones.k_min debe ser > 0 (la cola no se cierra)")
    if all(m.get(k) for k in ("D_mm", "cola_forma") + VARS_CUERPO[1:] + VARS_ALETA):
        for f in m["cola_forma"]:
            if f.get("forma") not in FORMAS_COLA:
                errores.append(f"malla.cola_forma: forma '{f.get('forma')}' no válida {FORMAS_COLA}")
        rangos = {"lambda_LE": (0, 1, "[0, 1)"), "mu_cr": (0, 1, "(0, 1]"), "gamma_ct": (0, 1, "(0, 1]"),
                  "sigma_flecha": (0, 1, "[0, 1]")}
        for v in m["lambda_LE"]:
            if not 0 <= v < 1:
                errores.append(f"malla.lambda_LE = {v} fuera de [0, 1)")
        for key in ("mu_cr", "gamma_ct"):
            for v in m[key]:
                if not 0 < v <= 1:
                    errores.append(f"malla.{key} = {v} fuera de {rangos[key][2]}")
        for v in m["sigma_flecha"]:
            if not 0 <= v <= 1:
                errores.append(f"malla.sigma_flecha = {v} fuera de [0, 1]")
        for v in m["k"]:
            if not v > 0:
                errores.append(f"malla.k = {v} debe ser > 0")
        for key in ("D_mm", "L_t_rel_D", "r_tip_rel_R"):
            for v in m[key]:
                if not v > 0:
                    errores.append(f"malla.{key} = {v} debe ser > 0")

    ob = raw["objetivo"]
    pesos = tuple(float(w) for w in ob.get("pesos", (1e6, 1e4, 1e2, 1.0)))
    if len(pesos) != 4 or not all(w > 0 for w in pesos):
        errores.append("objetivo.pesos: se esperan 4 pesos positivos")
    modo = ob.get("f4_modo", "centro")
    if modo not in ("centro", "max"):
        errores.append(f"objetivo.f4_modo '{modo}' no válido (centro | max)")
    if errores:
        raise ConfigError(errores)

    return ConfigOpt(
        raw=raw, ruta=path, base_raw=base, m_max=m_max, factor_llenado=fll, L=L, nariz=gf["nariz"],
        aleta=aleta, malla=m, restricciones=rest,
        objetivo=Objetivo(pesos=pesos, f4_modo=modo, SM_centro=float(ob.get("SM_centro_cal", 1.5))),
        ejecucion=raw.get("ejecucion") or {}, numerico=raw.get("numerico") or {},
        salida=raw.get("salida") or {}, raiz=raiz, base_ruta=base_path)


# --------------------------------------------------------------------------- puente con sensor_lastre


def raw_cuerpo(cfg: ConfigOpt, c: CuerpoSpec) -> dict:
    """YAML de `sensor_lastre` con una sola configuración: el cuerpo `c`.

    La aleta que lleva es solo un marcador válido (se reemplaza por cada candidato con
    `Caso.con_aletas`); su material y espesor ya son los de la optimización. El SM_max no se
    pasa al llenado (v3 §3.5): la sobreestabilidad se evalúa después, como restricción.
    """
    raw = copy.deepcopy(cfg.base_raw)
    a = cfg.aleta
    Lt_mm = c.L_t / MM
    raw["rellenos"] = [{"nombre": RELLENO, "material": "plomo"}]
    gb = raw["geometria_base"]
    gb["longitud_referencia"] = "cuerpo"
    gb["nariz"] = deep_merge(gb.get("nariz", {}), cfg.nariz)
    gb["cuerpo"] = {"longitud": "auto"}
    gb["cola"] = deep_merge(gb.get("cola", {}), {
        "forma": c.forma, "parametro": c.parametro, "popa": "abierta",
        "longitud": {"modo": "absoluto_mm", "valor": Lt_mm},
        "diametro_popa": {"modo": "absoluto_mm", "valor": c.k * c.D_mm},
    })
    gb["aletas"] = {
        "tipo": "freeform_cola", "n": a.n, "rotacion_deg": math.degrees(a.rotacion),
        "c_r": {"modo": "absoluto_mm", "valor": 0.5 * Lt_mm},
        "x_tip_le": {"modo": "absoluto_mm", "valor": 0.0},
        "h_tip": {"modo": "absoluto_mm", "valor": 0.25 * c.D_mm},
        "extension": {"modo": "absoluto_mm", "valor": 0.0},
        "r_interior_mm": 0.0, "offset_bottom_mm": 0.0, "epsilon_mm": 0.0,
        "espesor_mm": a.t / MM, "material": a.material, "fraccion_solida": a.phi,
    }
    raw["lastre"] = {**raw["lastre"], "factor_llenado": cfg.factor_llenado, "masa_max_sensor_g": cfg.m_max / G}
    raw["estabilidad"] = {**raw["estabilidad"], "SM_min_cal": cfg.restricciones.SM_min, "SM_max_cal": None}
    raw["presupuesto_masa"] = {"lastre_g": []}
    raw["envolvente"] = {**raw.get("envolvente", {}), "largo_max_mm": None, "alto_max_mm": None,
                         "ancho_max_mm": None}
    raw["barrido"] = {"modo": "cartesiano", "parametros": {"L_mm": [cfg.L / MM], "D_mm": [c.D_mm]}}
    raw["configuraciones"] = {}
    return raw
