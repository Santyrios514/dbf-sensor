#!/usr/bin/env python
"""Script 1: aplica cada configuración al .ork base con OpenRocket y exporta CP, C_Nα, masas,
perfil y polígono de aleta a CSV (spec v1 §5 + v2 §5).

    python scripts/01_orlab_export.py --config config/config.yaml [--ork ruta.ork]
        [--solo id1,id2] [--guardar-ork] [--regresion] [--verbose]
"""

from __future__ import annotations

import argparse
import math
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import pandas as pd  # noqa: E402

from sensor_lastre import esquemas  # noqa: E402
from sensor_lastre.aletas import envolvente  # noqa: E402
from sensor_lastre.config import ConfigError, cargar  # noqa: E402
from sensor_lastre.verificacion import comparar  # noqa: E402

MM, G = 1e-3, 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "config.yaml"))
    p.add_argument("--ork", help="sobrescribe openrocket.ork_base")
    p.add_argument("--solo", help="ids de configuración separados por comas")
    p.add_argument("--guardar-ork", action="store_true", help="guarda cada variante en data/ork/{id}.ork")
    p.add_argument("--regresion", action="store_true",
                   help="sin modificar la geometría: M de regresión y comparación con v2 §1.4")
    p.add_argument("--dir-datos", help="sobrescribe salida.dir_datos")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def _ruta(base: Path, v: str) -> Path:
    p = Path(v)
    return p if p.is_absolute() else base / p


def _escribir(filas, columnas, ruta: Path):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(filas).reindex(columns=columnas).to_csv(ruta, index=False, float_format="%.8g")


def regresion(pr, cfg) -> int:
    """v2 §1.4 / §8.5: masa, CG, CP y estabilidad del .ork sin tocar, a M = mach_regresion_or."""
    vuelo = cfg.casos[0].vuelo
    ref = cfg.verificacion["regresion_or"]
    aero = pr.aero(vuelo.mach_regresion_or, 0.0)
    m = pr.masas()
    D = float(pr.configuracion.getReferenceLength())
    obtenido = {"masa_g": m.m_total / G, "x_CG_mm": m.x_CG_total / MM, "x_CP_mm": aero.x_CP / MM,
                "estabilidad_cal": (aero.x_CP - m.x_CG_total) / D}
    ok = True
    print(f"Regresión OpenRocket {pr.version_or} (M = {vuelo.mach_regresion_or}, AoA = 0):")
    for k, v in obtenido.items():
        r = ref[k]
        pasa = abs(v - r["valor"]) <= r["tol"]
        ok &= pasa
        print(f"  {k:16s} {v:10.4f}   ref {r['valor']} ± {r['tol']}   {'OK' if pasa else 'FALLA'}")
    print(f"  advertencias: {aero.warnings}")
    return 0 if ok else 1


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    base = Path(a.config).resolve().parents[1]
    orc = cfg.openrocket
    ork = Path(a.ork) if a.ork else _ruta(base, orc["ork_base"])
    dir_datos = Path(a.dir_datos) if a.dir_datos else _ruta(base, cfg.salida.get("dir_datos", "data"))
    casos = cfg.casos
    if a.solo:
        ids = [s.strip() for s in a.solo.split(",")]
        casos = [c for c in casos if c.id in ids]

    import orlab
    from sensor_lastre.or_bridge import ErrorAleta, PuenteOR, estaciones_T_rho, sha256

    jar = orc.get("jar")
    tol_geo = cfg.verificacion.get("tol_geometria_ork_mm", 0.001) * MM
    with orlab.OpenRocketInstance(jar_path=jar) as inst:
        pr = PuenteOR(inst, orc["nombres_componentes"]).cargar(ork)
        if a.regresion:
            return regresion(pr, cfg)

        # 1. verificar que la base leída por OpenRocket coincide con la config (v2 §5.1)
        caso_base = next((c for c in cfg.casos if c.id == "base_ork"), None)
        if caso_base is not None:
            dif = comparar(pr.valores(), caso_base, tol_geo)
            if dif:
                print("La geometría del .ork no coincide con la configuración base:\n  - "
                      + "\n  - ".join(dif), file=sys.stderr)
                return 3
            print(f"Base verificada contra {ork.name} (OpenRocket {pr.version_or}).")
        else:
            print("Aviso: no hay configuración 'base_ork'; se omite la verificación de la base.",
                  file=sys.stderr)

        meta = {"or_version": pr.version_or, "orlab_version": pr.version_orlab(),
                "ork_sha256": sha256(ork), "fecha_iso": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        dif_cfg = orc.get("cna_diferencias") or {}
        delta = math.radians(dif_cfg.get("delta_aoa_deg", 1.0)) if dif_cfg.get("activar") else None
        resumen, perfiles, componentes, aletas = [], [], [], []
        for caso in casos:
            g, pa, vu = caso.geom, caso.params_aleta, caso.vuelo
            T = estaciones_T_rho(caso)
            fila = {
                "config_id": caso.id, "L_mm": g.L / MM, "D_mm": g.D / MM, "L_cuerpo_mm": g.L / MM,
                "nariz_forma": g.nariz.forma, "nariz_parametro": g.nariz.parametro,
                "L_nariz_mm": g.nariz.Ln / MM, "L_cola_mm": g.cola.Lt / MM, "D_popa_mm": 2 * g.cola.Ra / MM,
                "cola_forma": g.cola.forma, "cola_parametro": g.cola.parametro, "popa": g.cola.popa,
                "aletas_n": pa.n, "aletas_t_mm": pa.t / MM, "aleta_c_r_mm": pa.c_r / MM,
                "aleta_x_tip_le_mm": pa.x_s / MM, "aleta_extension_mm": pa.e / MM,
                "aleta_r_interior_mm": pa.r_in / MM, "aleta_rotacion_deg": math.degrees(pa.rotacion),
                **{f"T_pared_{e}_mm": T[e][0] / MM for e in T}, **{f"rho_eq_{e}_kg_m3": T[e][1] for e in T},
                "mach": vu.mach, "mach_usado": vu.mach, "aoa_deg": math.degrees(vu.aoa),
                "S_ref_mm2": math.pi * g.D**2 / 4 / MM**2, "D_ref_mm": g.D / MM, **meta,
            }
            try:
                pr.cargar(ork)  # cada variante parte del .ork base
                geo = pr.aplicar(caso)
                aero = pr.aero(vu.mach, vu.aoa, delta)
                m = pr.masas()
                P, orig, x_LE, r_LE = pr.poligono_aleta()
                h_env, w_env = envolvente(r_LE + geo["h"], g.R, pa.n, pa.rotacion)
                fila.update({
                    "estado": "ok", "mensaje": "", "L_total_mm": float(pr.configuracion.getLength()) / MM,
                    "aletas_x_r0_mm": x_LE / MM, "aleta_h_tip_mm": geo["h"] / MM, "r_LE_mm": r_LE / MM,
                    "r_tip_mm": (r_LE + geo["h"]) / MM, "x_CP_mm": aero.x_CP / MM, "CN_alpha_rad": aero.CNa,
                    "CN_alpha_dif_rad": aero.CNa_dif, "CD": aero.CD, "m_vacio_or_g": m.m_total / G,
                    "x_CG_vacio_or_mm": m.x_CG_total / MM, "h_env_mm": h_env / MM, "w_env_mm": w_env / MM,
                    "or_warnings": " | ".join(aero.warnings),
                })
                perfiles += [{"config_id": caso.id, "componente": k, "x_mm": x / MM, "r_ext_mm": r / MM}
                             for k, x, r in pr.perfil(int(orc.get("n_muestras_perfil", 800)))]
                for k in ("nariz", "cuerpo", "cola", "aletas"):
                    cna, xcp = aero.por_componente.get(k, (math.nan, math.nan))
                    mk, xk = m.por_componente[k]
                    componentes.append({"config_id": caso.id, "componente": k, "CN_alpha_rad": cna,
                                        "x_CP_mm": xcp / MM, "m_g": mk / G, "x_CG_mm": xk / MM})
                componentes.append({"config_id": caso.id, "componente": "total", "CN_alpha_rad": aero.CNa,
                                    "x_CP_mm": aero.x_CP / MM, "m_g": m.m_total / G,
                                    "x_CG_mm": m.x_CG_total / MM})
                aletas += [{"config_id": caso.id, "orden": i, "x_mm": x / MM, "r_mm": r / MM, "origen": o}
                           for i, ((x, r), o) in enumerate(zip(P, orig))]
                if a.guardar_ork:
                    destino = dir_datos / "ork" / f"{caso.id}.ork"
                    destino.parent.mkdir(parents=True, exist_ok=True)
                    pr.orl.save_doc(str(destino), pr.doc)
                print(f"{caso.id:>30}  x_CP = {aero.x_CP / MM:7.2f} mm  C_Nα = {aero.CNa:6.3f}  "
                      f"m = {m.m_total / G:6.1f} g  x_CG = {m.x_CG_total / MM:6.1f} mm")
            except ErrorAleta as e:
                fila.update({"estado": "error_aleta", "mensaje": str(e)})
                print(f"{caso.id:>30}  error_aleta: {e}", file=sys.stderr)
            except Exception as e:  # noqa: BLE001 - un fallo no detiene el barrido
                fila.update({"estado": "error", "mensaje": f"{type(e).__name__}: {e}"})
                print(f"{caso.id:>30}  error: {e}", file=sys.stderr)
                if a.verbose:
                    traceback.print_exc()
            resumen.append(fila)

    _escribir(resumen, esquemas.OR_RESUMEN, dir_datos / "or_resumen.csv")
    _escribir(perfiles, esquemas.OR_PERFILES, dir_datos / "or_perfiles.csv")
    _escribir(componentes, esquemas.OR_COMPONENTES, dir_datos / "or_componentes.csv")
    _escribir(aletas, esquemas.OR_ALETAS, dir_datos / "or_aletas.csv")
    n_ok = sum(1 for f in resumen if f.get("estado") == "ok")
    print(f"\n{n_ok}/{len(resumen)} configuraciones exportadas → {dir_datos}")
    return 0 if n_ok == len(resumen) else 1


if __name__ == "__main__":
    sys.exit(main())
