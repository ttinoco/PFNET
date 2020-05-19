/** @file gen.h
 *  @brief This file lists the constants and routines associated with the Gen data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __GEN_HEADER__
#define __GEN_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// Variables
/** \defgroup gen_vars Generator Variable Masks 
 *  @{
 */
#define GEN_VAR_P 0x01    /**< @brief Variable: generator active power */
#define GEN_VAR_Q 0x02    /**< @brief Variable: generator reactive power */
/** @} */

// Infinity
#define GEN_INF_P 1e8 /**< @brief Infinite active power */
#define GEN_INF_Q 1e8 /**< @brief Infinite reactive power */

// Properties
/** \defgroup gen_props Generator Property Masks 
 *  @{
 */
#define GEN_PROP_ANY 0x00       /**< @brief Property: any */
#define GEN_PROP_SLACK 0x01     /**< @brief Property: slack generator */
#define GEN_PROP_REG 0x02       /**< @brief Property: regulating generator */
#define GEN_PROP_NOT_REG 0x04   /**< @brief Property: non-regulating generator */
#define GEN_PROP_NOT_SLACK 0x08 /**< @brief Property: non-slack generator */
#define GEN_PROP_P_ADJUST 0x20  /**< @brief Property: P adjustable (Pmin < Pmax) */
#define GEN_PROP_REDISP 0x40    /**< @brief Property: P re-dispatchable */
/** @} */

// Constants
/** \defgroup gen_const Gen Constants
 *  @{
 */
#define GEN_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define GEN_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define GEN_NUM_JSON_FIELDS 20   /**< @brief Constant: max number of json fields */
/** @} */

// Generator
typedef struct Gen Gen;

// Other
typedef struct Bus Bus;

// Prototypes
void GEN_array_del(Gen* gen_array, int size);
void* GEN_array_get(void* gen_array, int index);
Gen* GEN_array_new(int size, int num_periods);
void GEN_array_show(Gen* gen_array, int size, int t);
void GEN_clear_sensitivities(Gen* gen);
void GEN_clear_flags(Gen* gen, char flag_type);
void GEN_copy_from_gen(Gen* gen, Gen* other);

char GEN_get_flags_vars(Gen* gen);
char GEN_get_flags_fixed(Gen* gen);
char GEN_get_flags_bounded(Gen* gen);
char GEN_get_flags_sparse(Gen* gen);

char* GEN_get_name(Gen* gen);
int GEN_get_num_periods(Gen* gen);

REAL GEN_get_sens_P_u_bound(Gen* gen, int t);
REAL* GEN_get_sens_P_u_bound_array(Gen* gen);
REAL GEN_get_sens_P_l_bound(Gen* gen, int t);
REAL* GEN_get_sens_P_l_bound_array(Gen* gen);
REAL GEN_get_sens_Q_u_bound(Gen* gen, int t);
REAL* GEN_get_sens_Q_u_bound_array(Gen* gen);
REAL GEN_get_sens_Q_l_bound(Gen* gen, int t);
REAL* GEN_get_sens_Q_l_bound_array(Gen* gen);

Bus* GEN_get_bus(Gen* gen);
Bus* GEN_get_reg_bus(Gen* gen);
REAL GEN_get_P_cost(Gen* gen, int t);
REAL GEN_get_P_cost_for(Gen* gen, REAL P);
REAL GEN_get_cost_coeff_Q0(Gen* gen);
REAL GEN_get_cost_coeff_Q1(Gen* gen);
REAL GEN_get_cost_coeff_Q2(Gen* gen);
char GEN_get_obj_type(void* gen);
int GEN_get_index(Gen* gen);
int GEN_get_index_P(Gen* gen, int t);
int GEN_get_index_Q(Gen* gen, int t);
int* GEN_get_index_P_array(Gen* gen);
int* GEN_get_index_Q_array(Gen* gen);
Gen* GEN_get_next(Gen* gen);
Gen* GEN_get_reg_next(Gen* gen);
REAL GEN_get_P(Gen* gen, int t);
REAL GEN_get_Q(Gen* gen, int t);
REAL* GEN_get_P_array(Gen* gen);
REAL* GEN_get_Q_array(Gen* gen);
REAL GEN_get_P_prev(Gen* gen);
REAL GEN_get_dP_max(Gen* gen);
REAL GEN_get_P_max(Gen* gen);
REAL GEN_get_P_min(Gen* gen);
REAL GEN_get_Q_max(Gen* gen);
REAL GEN_get_Q_min(Gen* gen);
REAL GEN_get_Q_par(Gen* gen);
void GEN_get_var_values(Gen* gen, Vec* values, int code);
char* GEN_get_var_info_string(Gen* gen, int index);
int GEN_get_num_vars(void* gen, unsigned char var, int t_start, int t_end);
Vec* GEN_get_var_indices(void* gen, unsigned char var, int t_start, int t_end);
char* GEN_get_json_string(Gen* gen, char* output);
BOOL GEN_has_flags(void* gen, char flag_type, unsigned char mask);
BOOL GEN_has_properties(void* gen, char prop);
void GEN_init(Gen* gen, int num_periods);
BOOL GEN_is_in_service(void* gen);
BOOL GEN_is_equal(Gen* gen, Gen* other);
BOOL GEN_is_P_adjustable(Gen* gen);
BOOL GEN_is_regulator(Gen* gen);
BOOL GEN_is_slack(Gen* gen);
BOOL GEN_is_redispatchable(Gen* gen);
Gen* GEN_list_add(Gen* gen_list, Gen* gen);
Gen* GEN_list_del(Gen* gen_list, Gen* gen);
int GEN_list_len(Gen* gen_list);
Gen* GEN_list_reg_add(Gen* reg_gen_list, Gen* reg_gen);
Gen* GEN_list_reg_del(Gen* reg_gen_list, Gen* reg_gen);
int GEN_list_reg_len(Gen* reg_gen_list);
Gen* GEN_new(int num_periods);
void GEN_propagate_data_in_time(Gen* gen, int start, int end);
void GEN_set_in_service(Gen* gen, BOOL in_service);
void GEN_set_redispatchable(Gen* gen, BOOL redisp);
void GEN_set_name(Gen* gen, char* name);
void GEN_set_sens_P_u_bound(Gen* gen, REAL value, int t);
void GEN_set_sens_P_l_bound(Gen* gen, REAL value, int t);
void GEN_set_sens_Q_u_bound(Gen* gen, REAL value, int t);
void GEN_set_sens_Q_l_bound(Gen* gen, REAL value, int t);
void GEN_set_cost_coeff_Q0(Gen* gen, REAL q);
void GEN_set_cost_coeff_Q1(Gen* gen, REAL q);
void GEN_set_cost_coeff_Q2(Gen* gen, REAL q);
void GEN_set_network(Gen* gen, void* net);
void GEN_set_bus(Gen* gen, Bus* bus);
void GEN_set_reg_bus(Gen* gen, Bus* reg_bus);
void GEN_set_index(Gen* gen, int index);
void GEN_set_P(Gen* gen, REAL P, int t);
void GEN_set_dP_max(Gen* gen, REAL P);
void GEN_set_P_max(Gen* gen, REAL P);
void GEN_set_P_min(Gen* gen, REAL P);
void GEN_set_P_prev(Gen* gen, REAL P);
void GEN_set_Q(Gen* gen, REAL Q, int t);
void GEN_set_Q_max(Gen* gen, REAL Q);
void GEN_set_Q_min(Gen* gen, REAL Q);
void GEN_set_Q_par(Gen* gen, REAL Q);
int GEN_set_flags(void* gen, char flag_type, unsigned char mask, int index);
void GEN_set_var_values(Gen* gen, Vec* values);
void GEN_show(Gen* gen, int t);

#endif
