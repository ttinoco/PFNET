/** @file branch.h
 *  @brief This file lists the constants and routines associated with the Branch data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BRANCH_HEADER__
#define __BRANCH_HEADER__

#include "types.h"
#include "list.h"
#include "vector.h"

// Branch types
#define BRANCH_TYPE_LINE 0       /**< @brief Type: transmission line */
#define BRANCH_TYPE_TRAN_FIXED 1 /**< @brief Type: fixed transformer */
#define BRANCH_TYPE_TRAN_TAP_V 2 /**< @brief Type: tap-changing transformer that regulates bus voltage magnitude */
#define BRANCH_TYPE_TRAN_TAP_Q 3 /**< @brief Type: tap-changing transformer that regulates reactive power flow */
#define BRANCH_TYPE_TRAN_PHASE 4 /**< @brief Type: phase-shifting transformer that regulates active power flow*/

// Branch flows
#define BRANCH_P_KM 0           /**< @brief Type: real power flow at the 'k' bus */
#define BRANCH_Q_KM 1           /**< @brief Type: reactive power flow at the 'k' bus */
#define BRANCH_P_KM_SERIES 2    /**< @brief Type: real power flow on the series element from 'k' to 'm' */
#define BRANCH_Q_KM_SERIES 3    /**< @brief Type: reactive power flow on the series element from 'k' to 'm' */
#define BRANCH_P_K_SHUNT 4      /**< @brief Type: real power flow on the shunt element from 'k' */
#define BRANCH_Q_K_SHUNT 5      /**< @brief Type: reactive power flow on the shunt element from 'k' */
#define BRANCH_P_MK 6           /**< @brief Type: real power flow at the 'm' bus */
#define BRANCH_Q_MK 7           /**< @brief Type: reactive power flow at the 'm' bus */
#define BRANCH_P_MK_SERIES 8    /**< @brief Type: real power flow on the series element from 'm' to 'k' */
#define BRANCH_Q_MK_SERIES 9    /**< @brief Type: reactive power flow on the series element from 'm' to 'k' */
#define BRANCH_P_M_SHUNT 10     /**< @brief Type: real power flow on the shunt element from 'm' */
#define BRANCH_Q_M_SHUNT 11     /**< @brief Type: reactive power flow on the shunt element from 'm' */
#define BRANCH_FLOW_SIZE 12     /**< @brief Type: the number of branch flows for a single branch including both directions */

// Variables
/** \defgroup branch_vars Branch Variable Masks
 *  @{
 */
#define BRANCH_VAR_RATIO 0x01     /**< @brief Variable: transformer taps ratio */
#define BRANCH_VAR_PHASE 0x02     /**< @brief Variable: transformer phase shift */
/** @} */

// Infinity
#define BRANCH_INF_RATIO 1e8 /**< @brief Infinite tap ratio */
#define BRANCH_INF_PHASE 1e8 /**< @brief Infinite phase shift */

// Properties
/** \defgroup branch_props Branch Property Masks
 *  @{
 */
#define BRANCH_PROP_ANY 0x00           /**< @brief Property: any */
#define BRANCH_PROP_TAP_CHANGER 0x01   /**< @brief Property: any tap-changing transformer */
#define BRANCH_PROP_TAP_CHANGER_V 0x02 /**< @brief Property: tap-changing transformer that regulates voltage */
#define BRANCH_PROP_TAP_CHANGER_Q 0x04 /**< @brief Property: tap-changing transformer that regulates reactive flow */
#define BRANCH_PROP_PHASE_SHIFTER 0x08 /**< @brief Property: phase-shiting transformer that regulates active flow */
#define BRANCH_PROP_NOT_OUT 0x10       /**< @brief Property: branch not on outage */
/** @} */

// Constants
/** \defgroup branch_const Branch Constants
 *  @{
 */
#define BRANCH_BUFFER_SIZE 100    /**< @brief Constant: buffer size for strings */
#define BRANCH_NUM_JSON_FIELDS 30 /**< @brief Constant: max number of json fields */
/** @} */

// Branch
typedef struct Branch Branch;

// Other
typedef struct Bus Bus;
typedef struct Vec Vec;

// Prototypes
void BRANCH_array_del(Branch* br_array, int size);
void* BRANCH_array_get(void* br, int index);
Branch* BRANCH_array_new(int size, int num_periods);
void BRANCH_array_show(Branch* br, int size, int t);
void BRANCH_clear_sensitivities(Branch* br);
void BRANCH_clear_flags(Branch* br, char flag_type);
void BRANCH_copy_from_branch(Branch* br, Branch* other);

char BRANCH_get_flags_vars(Branch* br);
char BRANCH_get_flags_fixed(Branch* br);
char BRANCH_get_flags_bounded(Branch* br);
char BRANCH_get_flags_sparse(Branch* br);

char* BRANCH_get_name(Branch* br);
int BRANCH_get_num_periods(Branch* br);
char BRANCH_get_type(Branch* br);
char BRANCH_get_obj_type(void* br);
REAL BRANCH_get_sens_P_u_bound(Branch* br, int t);
REAL BRANCH_get_sens_P_l_bound(Branch* br, int t);
REAL BRANCH_get_sens_ratio_u_bound(Branch* br, int t);
REAL BRANCH_get_sens_ratio_l_bound(Branch* br, int t);
REAL BRANCH_get_sens_phase_u_bound(Branch* br, int t);
REAL BRANCH_get_sens_phase_l_bound(Branch* br, int t);
int BRANCH_get_index(Branch* br);
int BRANCH_get_index_ratio(Branch* br, int t);
int BRANCH_get_index_phase(Branch* br, int t);
REAL BRANCH_get_ratio(Branch* br, int t);
REAL BRANCH_get_ratio_max(Branch* br);
REAL BRANCH_get_ratio_min(Branch* br);
REAL BRANCH_get_b(Branch* br);
REAL BRANCH_get_b_k(Branch* br);
REAL BRANCH_get_b_m(Branch* br);
REAL BRANCH_get_g(Branch* br);
REAL BRANCH_get_g_k(Branch* br);
REAL BRANCH_get_g_m(Branch* br);
Bus* BRANCH_get_bus_k(Branch* br);
Bus* BRANCH_get_bus_m(Branch* br);
Bus* BRANCH_get_reg_bus(Branch* br);
Branch* BRANCH_get_reg_next(Branch* br);
Branch* BRANCH_get_next_k(Branch* br);
Branch* BRANCH_get_next_m(Branch* br);
REAL BRANCH_get_phase(Branch* br, int t);
REAL BRANCH_get_phase_max(Branch* br);
REAL BRANCH_get_phase_min(Branch* br);
REAL BRANCH_get_P_max(Branch* br);
REAL BRANCH_get_P_min(Branch* br);
REAL BRANCH_get_Q_max(Branch* br);
REAL BRANCH_get_Q_min(Branch* br);
void BRANCH_compute_flows(Branch* br, Vec* var_values, int t, REAL* flows);
REAL BRANCH_get_i_km_mag(Branch* br, Vec* var_values, int t, REAL eps);
REAL BRANCH_get_i_mk_mag(Branch* br, Vec* var_values, int t, REAL eps);
REAL BRANCH_get_S_km_mag(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_S_mk_mag(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_km(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_km(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_mk(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_mk(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_km_series(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_km_series(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_mk_series(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_mk_series(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_k_shunt(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_k_shunt(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_P_m_shunt(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_Q_m_shunt(Branch* br, Vec* var_values, int t);
REAL BRANCH_get_ratingA(Branch* br);
REAL BRANCH_get_ratingB(Branch* br);
REAL BRANCH_get_ratingC(Branch* br);
REAL BRANCH_get_P_km_DC(Branch* br, int t);
REAL BRANCH_get_P_mk_DC(Branch* br, int t);
void BRANCH_get_var_values(Branch* br, Vec* values, int code);
char* BRANCH_get_var_info_string(Branch* br, int index);
int BRANCH_get_num_vars(void* br, unsigned char var, int t_start, int t_end);
Vec* BRANCH_get_var_indices(void* br, unsigned char var, int t_start, int t_end);
char* BRANCH_get_json_string(Branch* br, char* output);
BOOL BRANCH_has_flags(void* br, char flag_type, unsigned char mask);
BOOL BRANCH_has_pos_ratio_v_sens(Branch* br);
BOOL BRANCH_has_properties(void* br, char prop);
void BRANCH_init(Branch* br, int num_periods);
BOOL BRANCH_is_equal(Branch* br, Branch* other);
BOOL BRANCH_is_on_outage(Branch* br);
BOOL BRANCH_is_fixed_tran(Branch* br);
BOOL BRANCH_is_line(Branch* br);
BOOL BRANCH_is_phase_shifter(Branch* br);
BOOL BRANCH_is_tap_changer(Branch* br);
BOOL BRANCH_is_tap_changer_v(Branch* br);
BOOL BRANCH_is_tap_changer_Q(Branch* br);
Branch* BRANCH_list_reg_add(Branch* reg_br_list, Branch* br);
Branch* BRANCH_list_reg_del(Branch* reg_br_list, Branch* br);
int BRANCH_list_reg_len(Branch* reg_br_list);
Branch* BRANCH_list_k_add(Branch* k_br_list, Branch* br);
Branch* BRANCH_list_k_del(Branch* k_br_list, Branch* br);
int BRANCH_list_k_len(Branch* k_br_list);
Branch* BRANCH_list_m_add(Branch* m_br_list, Branch* br);
Branch* BRANCH_list_m_del(Branch* m_br_list, Branch* br);
int BRANCH_list_m_len(Branch* m_br_list);
Branch* BRANCH_new(int num_periods);
void BRANCH_propagate_data_in_time(Branch* br, int start, int end);
void BRANCH_set_name(Branch* br, char* name);
void BRANCH_set_outage(Branch* br, BOOL outage);
void BRANCH_set_sens_P_u_bound(Branch* br, REAL value, int t);
void BRANCH_set_sens_P_l_bound(Branch* br, REAL value, int t);
void BRANCH_set_sens_ratio_u_bound(Branch* br, REAL value, int t);
void BRANCH_set_sens_ratio_l_bound(Branch* br, REAL value, int t);
void BRANCH_set_sens_phase_u_bound(Branch* br, REAL value, int t);
void BRANCH_set_sens_phase_l_bound(Branch* br, REAL value, int t);
void BRANCH_set_index(Branch* br, int index);
void BRANCH_set_type(Branch* br, int type);
void BRANCH_set_bus_k(Branch* br, Bus* bus_k);
void BRANCH_set_bus_m(Branch* br, Bus* bus_m);
void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus);
void BRANCH_set_g(Branch* br, REAL g);
void BRANCH_set_g_k(Branch* br, REAL g_k);
void BRANCH_set_g_m(Branch* br, REAL g_m);
void BRANCH_set_b(Branch* br, REAL b);
void BRANCH_set_b_k(Branch* br, REAL b_k);
void BRANCH_set_b_m(Branch* br, REAL b_m);
void BRANCH_set_ratio(Branch* br, REAL ratio, int t);
void BRANCH_set_ratio_max(Branch* br, REAL ratio);
void BRANCH_set_ratio_min(Branch* br, REAL ratio);
void BRANCH_set_pos_ratio_v_sens(Branch* br, BOOL flag);
void BRANCH_set_phase(Branch* br, REAL phase, int t);
void BRANCH_set_phase_max(Branch* br, REAL phase);
void BRANCH_set_phase_min(Branch* br, REAL phase);
void BRANCH_set_P_max(Branch* br, REAL P_max);
void BRANCH_set_P_min(Branch* br, REAL P_min);
void BRANCH_set_Q_max(Branch* br, REAL Q_max);
void BRANCH_set_Q_min(Branch* br, REAL Q_min);
void BRANCH_set_ratingA(Branch* br, REAL r);
void BRANCH_set_ratingB(Branch* br, REAL r);
void BRANCH_set_ratingC(Branch* br, REAL r);
void BRANCH_set_var_values(Branch* br, Vec* values);
int BRANCH_set_flags(void* vbr, char flag_type, unsigned char mask, int index);
void BRANCH_show(Branch* br, int t);

#endif
