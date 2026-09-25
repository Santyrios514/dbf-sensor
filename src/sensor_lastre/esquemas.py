"""Esquemas contractuales de los CSV (sección 7). El orden de las columnas es parte del contrato."""

OR_RESUMEN = [
    "config_id", "estado", "mensaje", "L_mm", "D_mm", "nariz_forma", "nariz_parametro",
    "L_nariz_mm", "L_cola_mm", "D_popa_mm", "aletas_n", "aletas_cr_mm", "aletas_ct_mm",
    "aletas_s_mm", "aletas_flecha_mm", "aletas_t_mm", "aletas_x_r0_mm", "T_pared_nariz_mm",
    "T_pared_cuerpo_mm", "T_pared_cola_mm", "rho_eq_nariz_kg_m3", "rho_eq_cuerpo_kg_m3",
    "rho_eq_cola_kg_m3", "mach", "aoa_deg", "x_CP_mm", "CN_alpha_rad", "CN_alpha_dif_rad", "CD",
    "S_ref_mm2", "D_ref_mm", "m_vacio_or_g", "x_CG_vacio_or_mm", "or_warnings", "or_version",
    "orlab_version", "ork_sha256", "fecha_iso",
]

OR_PERFILES = ["config_id", "componente", "x_mm", "r_ext_mm"]

RESULTADOS = [
    "config_id", "relleno", "rho_b_kg_m3", "L_mm", "D_mm", "V_ext_cm3", "V_int_cm3", "x_b0_mm",
    "ell_geo_mm", "V_util_geo_cm3", "eta_int_geo", "eta_ext_geo", "m_relleno_geo_g",
    "m_total_geo_g", "x_CG_geo_mm", "SM_geo_cal", "ell_star_mm", "x_CG_min_mm",
    "SM_max_alcanzable_cal", "ell_SM_min_mm", "ell_SM_max_mm", "V_util_cm3", "eta_int", "eta_ext",
    "m_relleno_g", "m_total_g", "x_CG_mm", "SM_cal", "x_CP_mm", "x_CP_fuente",
    "dx_CP_or_vs_interno_mm", "CN_alpha_rad", "x_T_mm", "alpha_trim_geo_deg", "alpha_trim_deg",
    "x_T_trim_cero_mm", "m_casco_g", "m_aletas_g", "m_mamparos_g", "m_electronica_g",
    "m_puntuales_g", "cabe_largo", "cabe_diametro", "banderas",
]

MASAS_CAPAS = ["config_id", "estacion", "capa_idx", "material", "espesor_mm", "fraccion_solida",
               "m_g", "x_cg_mm"]

PRESUPUESTOS = ["config_id", "relleno", "m_lastre_obj_g", "ell_mm", "no_cabe", "m_total_g",
                "x_CG_mm", "SM_cal", "alpha_trim_deg"]

CURVA = ["ell_mm", "V_b_cm3", "m_b_g", "m_total_g", "x_CG_mm", "SM_cal", "alpha_trim_deg"]

# Complemento (no contractual): todos los intervalos de ℓ que cumplen el SM, útil cuando
# x_CG(ℓ) no es unimodal o la ventana está partida.
VENTANAS = ["config_id", "relleno", "desde_mm", "hasta_mm"]

BANDERAS = ("inviable_geo", "SM_inalcanzable", "no_unimodal", "inestable_respecto_remolque",
            "dif_CP_alta")
