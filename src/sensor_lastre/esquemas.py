"""Esquemas contractuales de los CSV (v1 §7 + v2 §7). El orden de las columnas es parte del contrato."""

OR_RESUMEN = [
    "config_id", "estado", "mensaje", "L_mm", "D_mm", "L_cuerpo_mm", "L_total_mm",
    "nariz_forma", "nariz_parametro", "L_nariz_mm", "L_cola_mm", "D_popa_mm", "cola_forma",
    "cola_parametro", "popa", "aletas_n", "aletas_t_mm", "aletas_x_r0_mm", "aleta_c_r_mm",
    "aleta_x_tip_le_mm", "aleta_h_tip_mm", "aleta_extension_mm", "aleta_r_interior_mm", "r_LE_mm",
    "r_tip_mm", "aleta_rotacion_deg", "T_pared_nariz_mm", "T_pared_cuerpo_mm", "T_pared_cola_mm",
    "rho_eq_nariz_kg_m3", "rho_eq_cuerpo_kg_m3", "rho_eq_cola_kg_m3", "mach", "mach_usado",
    "aoa_deg", "x_CP_mm", "CN_alpha_rad", "CN_alpha_dif_rad", "CD", "S_ref_mm2", "D_ref_mm",
    "m_vacio_or_g", "x_CG_vacio_or_mm", "h_env_mm", "w_env_mm", "or_warnings", "or_version",
    "orlab_version", "ork_sha256", "fecha_iso",
]

OR_PERFILES = ["config_id", "componente", "x_mm", "r_ext_mm"]

OR_COMPONENTES = ["config_id", "componente", "CN_alpha_rad", "x_CP_mm", "m_g", "x_CG_mm"]

OR_ALETAS = ["config_id", "orden", "x_mm", "r_mm", "origen"]

RESULTADOS = [
    "config_id", "relleno", "rho_b_kg_m3", "costo_usd_kg", "L_mm", "D_mm", "V_ext_cm3", "V_int_cm3", "x_b0_mm",
    "ell_geo_mm", "V_util_geo_cm3", "eta_int_geo", "eta_ext_geo", "m_relleno_geo_g",
    "m_total_geo_g", "x_CG_geo_mm", "SM_geo_cal", "ell_star_mm", "x_CG_min_mm",
    "SM_max_alcanzable_cal", "SM_inf_cal", "margen_SM_inf_cal", "ell_SM_min_mm", "ell_SM_max_mm",
    "V_util_cm3", "eta_int", "eta_ext", "m_relleno_g", "m_relleno_trasero_g",
    "costo_relleno_usd", "ell_trasero_mm", "limitante_masa", "m_total_g", "x_CG_mm", "SM_cal", "x_CP_mm",
    "x_CP_fuente", "dx_CP_or_vs_interno_mm", "CN_alpha_rad", "x_T_mm", "alpha_trim_geo_deg",
    "alpha_trim_deg", "x_T_trim_cero_mm", "m_casco_g", "A_aleta_mm2", "m_aletas_g",
    "x_CG_aletas_mm", "m_mamparos_g", "m_electronica_g", "m_puntuales_g", "dif_masa_or_pct",
    "L_total_mm", "D_max_aletas_mm", "h_env_mm", "w_env_mm", "cabe_largo", "cabe_alto", "cabe_ancho",
    "banderas",
]

MASAS_CAPAS = ["config_id", "estacion", "capa_idx", "material", "espesor_mm", "fraccion_solida",
               "m_g", "x_cg_mm"]

PRESUPUESTOS = ["config_id", "relleno", "m_lastre_obj_g", "ell_mm", "ell_trasero_mm", "no_cabe", "m_total_g",
                "x_CG_mm", "SM_cal", "alpha_trim_deg", "costo_usd"]

CURVA = ["ell_mm", "V_b_cm3", "m_b_g", "m_total_g", "x_CG_mm", "SM_cal", "alpha_trim_deg"]

# Ranking por el criterio de optimización: máxima masa de lastre con SM ≥ SM_min.
OPTIMO = ["relleno", "puesto", "config_id", "m_relleno_g", "m_relleno_trasero_g", "costo_usd_kg",
          "costo_relleno_usd", "m_total_g", "SM_cal", "limitante_masa", "ell_SM_max_mm", "ell_trasero_mm",
          "x_CG_mm", "x_CP_mm", "x_T_trim_cero_mm", "alpha_trim_deg", "D_max_aletas_mm", "h_env_mm",
          "w_env_mm", "L_total_mm", "cabe_alto", "cabe_ancho", "cabe_largo"]

# Una fila por configuración: dimensiones exteriores (con aletas) y dibujo.
CONFIGURACIONES = ["config_id", "L_mm", "L_total_mm", "D_mm", "D_max_aletas_mm", "h_env_mm", "w_env_mm",
                   "aletas_n", "aleta_h_tip_mm", "aleta_extension_mm", "aleta_rotacion_deg", "r_tip_mm",
                   "x_CP_mm", "x_CP_fuente", "dibujo"]

# Complemento (no contractual): todos los intervalos de ℓ que cumplen el SM, útil cuando
# x_CG(ℓ) no es unimodal o la ventana está partida.
VENTANAS = ["config_id", "relleno", "desde_mm", "hasta_mm"]

BANDERAS = ("inviable_geo", "SM_inalcanzable", "SM_inalcanzable_por_geometria", "margen_SM_bajo",
            "no_unimodal", "inestable_respecto_remolque", "dif_CP_alta", "dif_masa_alta")
