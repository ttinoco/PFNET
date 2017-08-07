/** @file load.h
 *  @brief This file lists the constants and routines associated with the Load data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __LOAD_HEADER__
#define __LOAD_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// Variables
/** \defgroup load_vars Load Variable Masks 
 *  @{
 */
#define LOAD_VAR_P 0x01    /**< @brief Variable: load active power */
#define LOAD_VAR_Q 0x02    /**< @brief Variable: load reactive power */
/** @} */

// Infinity
#define LOAD_INF_P 1e8 /**< @brief Infinite active power */
#define LOAD_INF_Q 1e8 /**< @brief Infinite reactive power */

// Others
#define LOAD_MIN_TARGET_PF 1e-4 /**< @brief Minimum target power factor */

// Properties
/** \defgroup load_props Load Property Masks 
 *  @{
 */
#define LOAD_PROP_ANY 0x00     /**< @brief Property: any */
#define LOAD_PROP_P_ADJUST 0x1 /**< @brief Property: P adjustable (Pmin < Pmax) */
/** @} */

// Constants
/** \defgroup load_const Load Constants
 *  @{
 */
#define LOAD_BUFFER_SIZE 100    /**< @brief Constant: buffer size for strings */
#define LOAD_NUM_JSON_FIELDS 15 /**< @brief Constant: max number of json fields */
/** @} */

// Load
typedef struct Load Load;

// Others
typedef struct Bus Bus;

void LOAD_array_del(Load* load_array, int size);
void* LOAD_array_get(void* load_array, int index);
Load* LOAD_array_new(int size, int num_periods);
void LOAD_array_show(Load* load_array, int num, int t);
void LOAD_clear_sensitivities(Load* load); 
void LOAD_clear_flags(Load* load, char flag_type);
void LOAD_propagate_data_in_time(Load* load, int start, int end);
int LOAD_get_num_periods(Load* load);
REAL LOAD_get_sens_P_u_bound(Load* load, int t);
REAL LOAD_get_sens_P_l_bound(Load* load, int t);
char LOAD_get_obj_type(void* load);
Bus* LOAD_get_bus(Load* load);
REAL LOAD_get_P_util(Load* load, int t);
REAL LOAD_get_P_util_for(Load* load, REAL P);
REAL LOAD_get_util_coeff_Q0(Load* load);
REAL LOAD_get_util_coeff_Q1(Load* load);
REAL LOAD_get_util_coeff_Q2(Load* load);
REAL LOAD_get_power_factor(Load* load, int t);
REAL LOAD_get_target_power_factor(Load* load);
int LOAD_get_index(Load* load);
int LOAD_get_index_P(Load* load, int t);
int LOAD_get_index_Q(Load* load, int t);
Load* LOAD_get_next(Load* load);
REAL LOAD_get_P(Load* load, int t);
REAL LOAD_get_P_max(Load* load, int t);
REAL LOAD_get_P_min(Load* load, int t);
REAL LOAD_get_Q(Load* load, int t);
void LOAD_get_var_values(Load* load, Vec* values, int code);
int LOAD_get_num_vars(void* load, unsigned char var, int t_start, int t_end);
Vec* LOAD_get_var_indices(void* load, unsigned char var, int t_start, int t_end);
char* LOAD_get_json_string(Load* load, char* output);
BOOL LOAD_has_flags(void* load, char flag_type, unsigned char mask);
BOOL LOAD_has_properties(void* load, char prop);
void LOAD_init(Load* load, int num_periods);
BOOL LOAD_is_P_adjustable(Load* load);
Load* LOAD_list_add(Load *load_list, Load* load);
int LOAD_list_len(Load* load_list);
Load* LOAD_new(int num_periods);
void LOAD_set_target_power_factor(Load* load, REAL pf);
void LOAD_set_sens_P_u_bound(Load* load, REAL value, int t);
void LOAD_set_sens_P_l_bound(Load* load, REAL value, int t);
void LOAD_set_util_coeff_Q0(Load* load, REAL q);
void LOAD_set_util_coeff_Q1(Load* load, REAL q);
void LOAD_set_util_coeff_Q2(Load* load, REAL q);
void LOAD_set_bus(Load* load, Bus* bus);
void LOAD_set_index(Load* load, int index);
void LOAD_set_P(Load* load, REAL P, int t);
void LOAD_set_P_max(Load* load, REAL P, int t);
void LOAD_set_P_min(Load* load, REAL P, int t);
void LOAD_set_Q(Load* load, REAL Q, int t);
int LOAD_set_flags(void* load, char flag_type, unsigned char mask, int index);
void LOAD_set_var_values(Load* load, Vec* values);
void LOAD_show(Load* load, int t);

#endif
