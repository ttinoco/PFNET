/** @file branch.h
 *  @brief This file lists the constants and routines associated with the Branch data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
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

// Variables
/** \defgroup branch_vars Branch Variable Masks
 *  @{
 */
#define BRANCH_VAR_RATIO 0x01     /**< @brief Variable: transformer taps ratio */
#define BRANCH_VAR_RATIO_DEV 0x02 /**< @brief Variable: branch taps ratio pos/neg deviations from current value */
#define BRANCH_VAR_PHASE 0x04     /**< @brief Variable: transformer phase shift */
#define BRANCH_VAR_P 0x08         /**< @brief Variable: branch active flow */
#define BRANCH_VAR_Q 0x10         /**< @brief Variable: branch ractive flow */
/** @} */

// Infinity
#define BRANCH_INF_RATIO 100. /**< @brief Infinite tap ratio */
#define BRANCH_INF_FLOW 1000. /**< @brief Infinite power flow (p.u.) */

// Properties
/** \defgroup branch_props Branch Property Masks
 *  @{
 */
#define BRANCH_PROP_ANY 0x00           /**< @brief Property: any */
#define BRANCH_PROP_TAP_CHANGER 0x01   /**< @brief Property: any tap-changing transformer */
#define BRANCH_PROP_TAP_CHANGER_V 0x02 /**< @brief Property: tap-changing transformer that regulates voltage */
#define BRANCH_PROP_TAP_CHANGER_Q 0x04 /**< @brief Property: tap-changing transformer that regulates reactive flow */
#define BRANCH_PROP_PHASE_SHIFTER 0x08 /**< @brief Property: phase-shiting transformer that regulates active flow */
/** @} */

// Branch
typedef struct Branch Branch;

// Other
typedef struct Bus Bus;
typedef struct Vec Vec;

// Prototypes
void* BRANCH_array_get(void* br, int index);
Branch* BRANCH_array_new(int num);
void BRANCH_array_show(Branch* br, int num);
void BRANCH_clear_sensitivities(Branch* br);
void BRANCH_clear_flags(Branch* br, char flag_type);
char BRANCH_get_type(Branch* br);
char BRANCH_get_obj_type(void* br);
REAL BRANCH_get_sens_P_u_bound(Branch* br);
REAL BRANCH_get_sens_P_l_bound(Branch* br);
int BRANCH_get_index(Branch* br);
int BRANCH_get_index_ratio(Branch* br);
int BRANCH_get_index_ratio_y(Branch* br);
int BRANCH_get_index_ratio_z(Branch* br);
int BRANCH_get_index_phase(Branch* br);
REAL BRANCH_get_ratio(Branch* br);
REAL BRANCH_get_ratio_max(Branch* br);
REAL BRANCH_get_ratio_min(Branch* br);
REAL BRANCH_get_b(Branch* br);
REAL BRANCH_get_b_from(Branch* br);
REAL BRANCH_get_b_to(Branch* br);
REAL BRANCH_get_g(Branch* br);
REAL BRANCH_get_g_from(Branch* br);
REAL BRANCH_get_g_to(Branch* br);
Bus* BRANCH_get_bus_from(Branch* br);
Bus* BRANCH_get_bus_to(Branch* br);
Bus* BRANCH_get_reg_bus(Branch* br);
Branch* BRANCH_get_reg_next(Branch* br);
Branch* BRANCH_get_from_next(Branch* br);
Branch* BRANCH_get_to_next(Branch* br);
REAL BRANCH_get_phase(Branch* br);
REAL BRANCH_get_phase_max(Branch* br);
REAL BRANCH_get_phase_min(Branch* br);
REAL BRANCH_get_ratingA(Branch* br);
REAL BRANCH_get_ratingB(Branch* br);
REAL BRANCH_get_ratingC(Branch* br);
REAL BRANCH_get_P_flow_DC(Branch* br);
void BRANCH_get_var_values(Branch* br, Vec* values, int code);
int BRANCH_get_var_index(void* br, char var);
BOOL BRANCH_has_flags(void* br, char flag_type, char mask);
BOOL BRANCH_has_pos_ratio_v_sens(Branch* br);
BOOL BRANCH_has_properties(void* br, char prop);
void BRANCH_init(Branch* br);
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
Branch* BRANCH_list_from_add(Branch* from_br_list, Branch* br);
Branch* BRANCH_list_from_del(Branch* from_br_list, Branch* br);
int BRANCH_list_from_len(Branch* from_br_list);
Branch* BRANCH_list_to_add(Branch* to_br_list, Branch* br);
Branch* BRANCH_list_to_del(Branch* to_br_list, Branch* br);
int BRANCH_list_to_len(Branch* to_br_list);
Branch* BRANCH_new(void);
void BRANCH_set_outage(Branch* branch, BOOL outage);
void BRANCH_set_sens_P_u_bound(Branch* br, REAL value);
void BRANCH_set_sens_P_l_bound(Branch* br, REAL value);
void BRANCH_set_index(Branch* br, int index);
void BRANCH_set_type(Branch* br, int type);
void BRANCH_set_bus_from(Branch* br, Bus* bus_from);
void BRANCH_set_bus_to(Branch* br, Bus* bus_to);
void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus);
void BRANCH_set_g(Branch* br, REAL g);
void BRANCH_set_g_from(Branch* br, REAL g_from);
void BRANCH_set_g_to(Branch* br, REAL g_to);
void BRANCH_set_b(Branch* br, REAL b);
void BRANCH_set_b_from(Branch* br, REAL b_from);
void BRANCH_set_b_to(Branch* br, REAL b_to);
void BRANCH_set_ratio(Branch* br, REAL ratio);
void BRANCH_set_ratio_max(Branch* br, REAL ratio);
void BRANCH_set_ratio_min(Branch* br, REAL ratio);
void BRANCH_set_pos_ratio_v_sens(Branch* br, BOOL flag);
void BRANCH_set_phase(Branch* br, REAL phase);
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
int BRANCH_set_flags(void* vbr, char flag_type, char mask, int index);
void BRANCH_show(Branch* br);

#endif
