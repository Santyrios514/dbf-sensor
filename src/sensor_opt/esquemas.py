"""Columnas de las salidas de sensor_opt (v3 §7). Las extra van después de las de la spec."""

RANKING = [
    # spec v3 §7
    "cand_id", "factible", "motivos", "J", "f1", "f2", "f3", "f4", "D_mm", "cola_forma", "cola_parametro",
    "L_t_mm", "k", "theta_eq_deg", "fineza_cola", "f_base_roma", "k_efectivo", "lambda_LE", "x_LE_mm", "c_r_mm", "r_LE_mm", "r_TE_raiz_mm", "r_tip_mm",
    "h_mm", "c_t_mm", "x_s_mm", "A_aleta_mm2", "V_flutter_m_s", "flutter_margen_bajo", "m_aletas_g", "D_ap_mm",
    "h_45_mm", "m_total_g", "m_lastre_g", "ell_mm", "restriccion_activa", "x_CG_mm", "x_CP_sust_mm",
    "x_CP_cal_mm", "SM_cal", "CN_alpha_total", "alpha_trim_deg", "en_pareto", "verificado_or", "x_CP_or_mm",
    "SM_or_cal", "dif_masa_or_pct",
    # extra: electrónica, trim y amarre (código actual), distribución del lastre, malla
    "puesto", "cuerpo_id", "refinamiento", "L_t_rel_D", "mu_cr", "r_tip_rel_R", "gamma_ct", "sigma_flecha",
    "L_c_mm", "x_b0_mm", "ell_geo_mm", "m_casco_g", "m_electronica_g", "m_lastre_delantero_g",
    "m_lastre_trasero_g", "ell_trasero_mm", "x_T_mm", "x_T_trim_cero_mm", "tol_amarre_mm", "banderas_lastre",
    "J_or", "m_total_or_g", "x_CG_or_mm",
]

VERIFICACION = [
    "cand_id", "rol", "iteracion", "x_CP_sust_mm", "x_CP_cal_mm", "x_CP_or_mm", "dx_or_sust_mm", "dx_or_cal_mm",
    "CNa_or", "CNa_sust", "CNa_or_nariz", "x_CP_or_nariz_mm", "CNa_or_cuerpo", "x_CP_or_cuerpo_mm",
    "CNa_or_cola", "x_CP_or_cola_mm", "CNa_or_aletas", "x_CP_or_aletas_mm",
    "m_or_g", "x_CG_or_mm", "m_modelo_g", "x_CG_modelo_mm", "dif_masa_or_pct",
    "m_or_nariz_g", "m_or_cuerpo_g", "m_or_cola_g", "m_or_aletas_g", "SM_or_cal", "factible_or", "advertencias",
]
