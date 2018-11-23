/** @file facts.h
 *  @brief This file lists the constants and routines associated with the Facts data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FACTS_HEADER__
#define __FACTS_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// Variables
/** \defgroup facts_vars Facts Variable Masks 
 *  @{
 */
#define FACTS_VAR_VMAG_S 0x01    /**< @brief Variable: series voltage magnitude */
#define FACTS_VAR_VANG_S 0x02    /**< @brief Variable: series voltage angle */
#define FACTS_VAR_P 0x04         /**< @brief Variable: active power */
#define FACTS_VAR_Q 0x08         /**< @brief Variable: reactive power */
/** @} */

// Infinity
#define FACTS_INF_VMAG_S 1e8 /**< @brief Infinite series voltage magnitude (p.u.) */
#define FACTS_INF_VANG_S 1e2 /**< @brief Infinite series voltage angle (p.u.) */
#define FACTS_INF_P 1e8      /**< @brief Infinite active power (p.u.) */
#define FACTS_INF_Q 1e8      /**< @brief Infinite active power (p.u.) */

// Series link modes
#define FACTS_SERIES_MODE_DISABLED 0
#define FACTS_SERIES_MODE_NORMAL 1
#define FACTS_SERIES_MODE_BYPASS 2
#define FACTS_SERIES_MODE_CZ 3
#define FACTS_SERIES_MODE_CV 4

// Properties
/** \defgroup facts_props Facts Property Masks 
 *  @{
 */
#define FACTS_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup facts_const Facts Constants
 *  @{
 */
#define FACTS_BUFFER_SIZE 100    /**< @brief Constant: buffer size for strings */
#define FACTS_NUM_JSON_FIELDS 30 /**< @brief Constant: max number of json fields */
/** @} */

// Facts
typedef struct Facts Facts;

// Others
typedef struct Bus Bus;

void FACTS_array_del(Facts* facts_array, int size);
void* FACTS_array_get(void* facts_array, int index);
Facts* FACTS_array_new(int size, int num_periods);
void FACTS_array_show(Facts* facts_array, int num, int t);
void FACTS_clear_sensitivities(Facts* facts); 
void FACTS_clear_flags(Facts* facts, char flag_type);
void FACTS_copy_from_facts(Facts* facts, Facts* other);

char FACTS_get_flags_vars(Facts* facts);
char FACTS_get_flags_fixed(Facts* facts);
char FACTS_get_flags_bounded(Facts* facts);
char FACTS_get_flags_sparse(Facts* facts);

char* FACTS_get_name(Facts* facts);
int FACTS_get_num_periods(Facts* facts);
char FACTS_get_obj_type(void* facts);
Bus* FACTS_get_reg_bus(Facts* facts);
Bus* FACTS_get_bus_k(Facts* facts);
Bus* FACTS_get_bus_m(Facts* facts);
int FACTS_get_index(Facts* facts);

int FACTS_get_index_v_mag_s(Facts* facts, int t);
int FACTS_get_index_v_ang_s(Facts* facts, int t);
int FACTS_get_index_P_k(Facts* facts, int t);
int FACTS_get_index_P_m(Facts* facts, int t);
int FACTS_get_index_Q_k(Facts* facts, int t);
int FACTS_get_index_Q_m(Facts* facts, int t);
int FACTS_get_index_P_dc(Facts* facts, int t);
int FACTS_get_index_Q_s(Facts* facts, int t);
int FACTS_get_index_Q_sh(Facts* facts, int t);

int* FACTS_get_index_v_mag_s_array(Facts* facts);
int* FACTS_get_index_v_ang_s_array(Facts* facts);
int* FACTS_get_index_P_k_array(Facts* facts);
int* FACTS_get_index_P_m_array(Facts* facts);
int* FACTS_get_index_Q_k_array(Facts* facts);
int* FACTS_get_index_Q_m_array(Facts* facts);
int* FACTS_get_index_P_dc_array(Facts* facts);
int* FACTS_get_index_Q_s_array(Facts* facts);
int* FACTS_get_index_Q_sh_array(Facts* facts);

Facts* FACTS_get_reg_next(Facts* facts);
Facts* FACTS_get_next_k(Facts* facts);
Facts* FACTS_get_next_m(Facts* facts);

REAL FACTS_get_v_max_s(Facts* facts);
REAL FACTS_get_g(Facts* facts);
REAL FACTS_get_b(Facts* facts);
char FACTS_get_mode_s(Facts* facts);
REAL FACTS_get_Q_par(Facts* facts);
REAL FACTS_get_Q_max_s(Facts* facts);
REAL FACTS_get_Q_max_sh(Facts* facts);
REAL FACTS_get_Q_min_s(Facts* facts);
REAL FACTS_get_Q_min_sh(Facts* facts);
REAL FACTS_get_i_max_s(Facts* facts);
REAL FACTS_get_i_max_sh(Facts* facts);
REAL FACTS_get_P_max_dc(Facts* facts);
REAL FACTS_get_v_min_m(Facts* facts);
REAL FACTS_get_v_max_m(Facts* facts);

REAL FACTS_get_v_mag_s(Facts* facts, int t);
REAL FACTS_get_v_ang_s(Facts* facts, int t);
REAL FACTS_get_P_k(Facts* facts, int t);
REAL FACTS_get_P_m(Facts* facts, int t);
REAL FACTS_get_Q_k(Facts* facts, int t);
REAL FACTS_get_Q_m(Facts* facts, int t);
REAL FACTS_get_Q_sh(Facts* facts, int t);
REAL FACTS_get_Q_s(Facts* facts, int t);
REAL FACTS_get_P_dc(Facts* facts, int t);
REAL FACTS_get_P_set(Facts* facts, int t);
REAL FACTS_get_Q_set(Facts* facts, int t);

REAL* FACTS_get_v_mag_s_array(Facts* facts);
REAL* FACTS_get_v_ang_s_array(Facts* facts);
REAL* FACTS_get_P_k_array(Facts* facts);
REAL* FACTS_get_P_m_array(Facts* facts);
REAL* FACTS_get_Q_k_array(Facts* facts);
REAL* FACTS_get_Q_m_array(Facts* facts);
REAL* FACTS_get_Q_sh_array(Facts* facts);
REAL* FACTS_get_Q_s_array(Facts* facts);
REAL* FACTS_get_P_dc_array(Facts* facts);
REAL* FACTS_get_P_set_array(Facts* facts);
REAL* FACTS_get_Q_set_array(Facts* facts);

void FACTS_get_var_values(Facts* facts, Vec* values, int code);
char* FACTS_get_var_info_string(Facts* facts, int index);
int FACTS_get_num_vars(void* facts, unsigned char var, int t_start, int t_end);
Vec* FACTS_get_var_indices(void* facts, unsigned char var, int t_start, int t_end);
char* FACTS_get_json_string(Facts* facts, char* output);
BOOL FACTS_has_flags(void* facts, char flag_type, unsigned char mask);
BOOL FACTS_has_properties(void* facts, char prop);
void FACTS_init(Facts* facts, int num_periods);

BOOL FACTS_is_equal(Facts* facts, Facts* other);
BOOL FACTS_is_regulator(Facts* facts);
BOOL FACTS_is_STATCOM(Facts* facts);
BOOL FACTS_is_SSSC(Facts* facts);
BOOL FACTS_is_UPFC(Facts* facts);
BOOL FACTS_is_series_link_disabled(Facts* facts);
BOOL FACTS_is_series_link_bypassed(Facts* facts);
BOOL FACTS_is_in_normal_series_mode(Facts* facts);
BOOL FACTS_is_in_constant_series_z_mode(Facts* facts);
BOOL FACTS_is_in_constant_series_v_mode(Facts* facts);

Facts* FACTS_list_k_add(Facts* facts_list, Facts* facts);
Facts* FACTS_list_k_del(Facts* facts_list, Facts* facts);
int FACTS_list_k_len(Facts* facts_list);
Facts* FACTS_list_m_add(Facts* facts_list, Facts* facts);
Facts* FACTS_list_m_del(Facts* facts_list, Facts* facts);
int FACTS_list_m_len(Facts* facts_list);
Facts* FACTS_list_reg_add(Facts* facts_list, Facts* facts);
Facts* FACTS_list_reg_del(Facts* facts_list, Facts* facts);
int FACTS_list_reg_len(Facts* facts_list);

Facts* FACTS_new(int num_periods);
void FACTS_propagate_data_in_time(Facts* facts, int start, int end);

void FACTS_set_name(Facts* facts, char* name);
void FACTS_set_reg_bus(Facts* facts, Bus* bus);
void FACTS_set_bus_k(Facts* facts, Bus* bus);
void FACTS_set_bus_m(Facts* facts, Bus* bus);
void FACTS_set_index(Facts* facts, int index);
void FACTS_set_v_mag_s(Facts* facts, REAL v_mag_s, int t);
void FACTS_set_v_ang_s(Facts* facts, REAL v_ang_s, int t);
void FACTS_set_v_max_s(Facts* facts, REAL v_max);
void FACTS_set_mode_s(Facts* facts, char mode);
void FACTS_set_P_k(Facts* facts, REAL P, int t);
void FACTS_set_P_m(Facts* facts, REAL P, int t);
void FACTS_set_Q_k(Facts* facts, REAL Q, int t);
void FACTS_set_Q_m(Facts* facts, REAL Q, int t);
void FACTS_set_Q_s(Facts* facts, REAL Q, int t);
void FACTS_set_Q_sh(Facts* facts, REAL Q, int t);
void FACTS_set_P_dc(Facts* facts, REAL P, int t);
void FACTS_set_Q_par(Facts* facts, REAL Q_par);
void FACTS_set_P_set(Facts* facts, REAL P, int t);
void FACTS_set_Q_set(Facts* facts, REAL Q, int t);
void FACTS_set_Q_max_s(Facts* facts, REAL Q_max);
void FACTS_set_Q_max_sh(Facts* facts, REAL Q_max);
void FACTS_set_Q_min_s(Facts* facts, REAL Q_min);
void FACTS_set_Q_min_sh(Facts* facts, REAL Q_min);
void FACTS_set_i_max_s(Facts* facts, REAL i_max);
void FACTS_set_i_max_sh(Facts* facts, REAL i_max);
void FACTS_set_P_max_dc(Facts* facts, REAL P_max);
void FACTS_set_v_min_m(Facts* facts, REAL v_min);
void FACTS_set_v_max_m(Facts* facts, REAL v_max);
void FACTS_set_g(Facts* facts, REAL g);
void FACTS_set_b(Facts* facts, REAL b);

int FACTS_set_flags(void* facts, char flag_type, unsigned char mask, int index);
void FACTS_set_var_values(Facts* facts, Vec* values);
void FACTS_show(Facts* facts, int t);

#endif
