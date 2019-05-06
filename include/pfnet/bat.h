/** @file bat.h
 *  @brief This file lists the constants and routines associated with the Bat data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BAT_HEADER__
#define __BAT_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// Variables
/** \defgroup bat_vars Battery Variable Masks 
 *  @{
 */
#define BAT_VAR_P 0x01    /**< @brief Variable: battery charing/discharging power */
#define BAT_VAR_E 0x02    /**< @brief Variable: battery energy level */
/** @} */

// Infinity
#define BAT_INF_P 1e8 /**< @brief Infinite charging/discharging power */
#define BAT_INF_E 1e8 /**< @brief Infinite energy level */

// Properties
/** \defgroup bat_props Battery Property Masks 
 *  @{
 */
#define BAT_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup bat_const Bat Constants
 *  @{
 */
#define BAT_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define BAT_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define BAT_NUM_JSON_FIELDS 15   /**< @brief Constant: max number of json fields */
/** @} */

// Battery
typedef struct Bat Bat;

// Others
typedef struct Bus Bus;

void BAT_array_del(Bat* bat_array, int size);
void* BAT_array_get(void* bat_array, int index);
Bat* BAT_array_new(int size, int num_periods);
void BAT_array_show(Bat* bat_array, int size, int t);
void BAT_clear_sensitivities(Bat* bat);
void BAT_clear_flags(Bat* bat, char flag_type);
void BAT_copy_from_bat(Bat* bat, Bat* other);

char BAT_get_flags_vars(Bat* bat);
char BAT_get_flags_fixed(Bat* bat);
char BAT_get_flags_bounded(Bat* bat);
char BAT_get_flags_sparse(Bat* bat);

char* BAT_get_name(Bat* bat);
int BAT_get_num_periods(Bat* bat);
char BAT_get_obj_type(void* bat);
Bus* BAT_get_bus(Bat* bat);
int BAT_get_index(Bat* bat);
int BAT_get_index_Pc(Bat* bat, int t);
int BAT_get_index_Pd(Bat* bat, int t);
int BAT_get_index_E(Bat* bat, int t);
int* BAT_get_index_Pc_array(Bat* bat);
int* BAT_get_index_Pd_array(Bat* bat);
int* BAT_get_index_E_array(Bat* bat);
Bat* BAT_get_next(Bat* bat);
REAL BAT_get_P(Bat* bat, int t);
REAL* BAT_get_P_array(Bat* bat);
REAL BAT_get_P_max(Bat* bat);
REAL BAT_get_P_min(Bat* bat);
REAL BAT_get_E(Bat* bat, int t);
REAL* BAT_get_E_array(Bat* bat);
REAL BAT_get_E_init(Bat* bat);
REAL BAT_get_E_final(Bat* bat);
REAL BAT_get_E_max(Bat* bat);
REAL BAT_get_eta_c(Bat* bat);
REAL BAT_get_eta_d(Bat* bat);
void BAT_get_var_values(Bat* bat, Vec* values, int code);
char* BAT_get_var_info_string(Bat* bat, int index);
int BAT_get_num_vars(void* bat, unsigned char var, int t_start, int t_end);
Vec* BAT_get_var_indices(void* bat, unsigned char var, int t_start, int t_end);
char* BAT_get_json_string(Bat* bat, char* output);
BOOL BAT_has_flags(void* bat, char flag_type, unsigned char mask);
BOOL BAT_has_properties(void* bat, char prop);
void BAT_init(Bat* bat, int num_periods);
BOOL BAT_is_in_service(void* bat);
BOOL BAT_is_equal(Bat* bat, Bat* other);
Bat* BAT_list_add(Bat* bat_list, Bat* bat);
Bat* BAT_list_del(Bat* bat_list, Bat* bat);
int BAT_list_len(Bat* bat_list);
Bat* BAT_new(int num_periods);
void BAT_propagate_data_in_time(Bat* bat, int start, int end);
void BAT_set_network(Bat* bat, void* net);
void BAT_set_in_service(Bat* bat, BOOL in_service);
void BAT_set_name(Bat* bat, char* name);
void BAT_set_bus(Bat* bat, Bus* bus);
void BAT_set_index(Bat* bat, int index);
void BAT_set_P(Bat* bat, REAL P, int t);
void BAT_set_P_max(Bat* bat, REAL P);
void BAT_set_P_min(Bat* bat, REAL P);
void BAT_set_E(Bat* bat, REAL E, int t);
void BAT_set_E_init(Bat* bat, REAL E);
void BAT_set_E_final(Bat* bat, REAL E);
void BAT_set_E_max(Bat* bat, REAL E);
void BAT_set_eta_c(Bat* bat, REAL eta_c);
void BAT_set_eta_d(Bat* bat, REAL eta_d);
int BAT_set_flags(void* bat, char flag_type, unsigned char mask, int index);
void BAT_set_var_values(Bat* bat, Vec* values);
void BAT_show(Bat* bat, int t);

#endif
