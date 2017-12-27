/** @file vargen.h
 *  @brief This file lists the constants and routines associated with the Vargen data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __VARGEN_HEADER__
#define __VARGEN_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// Variables
/** \defgroup vargen_vars Variable Generator Variable Masks 
 *  @{
 */
#define VARGEN_VAR_P 0x01    /**< @brief Variable: variable generator active power */
#define VARGEN_VAR_Q 0x02    /**< @brief Variable: variable generator reactive power */
/** @} */

// Infinity
#define VARGEN_INF_P 1e8 /**< @brief Infinite active power */
#define VARGEN_INF_Q 1e8 /**< @brief Infinite reactive power */

// Variable generator types
#define VARGEN_TYPE_WIND 0       /**< @brief Type: wind farm */
#define VARGEN_TYPE_SOLAR 1       /**< @brief Type: solar plant */

// Properties
/** \defgroup vargen_props Variable Generator Property Masks 
 *  @{
 */
#define VARGEN_PROP_ANY 0x00       /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup vargen_const Variable Generator Constants
 *  @{
 */
#define VARGEN_BUFFER_SIZE 100    /**< @brief Constant: buffer size for name */
#define VARGEN_NUM_JSON_FIELDS 15 /**< @brief Constant: max number of json fields */
/** @} */

// Variable generator
typedef struct Vargen Vargen;

// Other
typedef struct Bus Bus;

// Prototypes
void VARGEN_array_del(Vargen* gen_array, int size);
void* VARGEN_array_get(void* gen_array, int index);
Vargen* VARGEN_array_new(int size, int num_periods);
void VARGEN_array_show(Vargen* gen_array, int size, int t);
void VARGEN_clear_sensitivities(Vargen* vargen);
void VARGEN_clear_flags(Vargen* gen, char flag_type);
void VARGEN_copy_from_vargen(Vargen* vargen, Vargen* other);

char VARGEN_get_flags_vars(Vargen* gen);
char VARGEN_get_flags_fixed(Vargen* gen);
char VARGEN_get_flags_bounded(Vargen* gen);
char VARGEN_get_flags_sparse(Vargen* gen);

int VARGEN_get_num_periods(Vargen* gen);
char* VARGEN_get_name(Vargen* gen);
Bus* VARGEN_get_bus(Vargen* gen);
char VARGEN_get_obj_type(void* gen);
int VARGEN_get_index(Vargen* gen);
int VARGEN_get_index_P(Vargen* gen, int t);
int VARGEN_get_index_Q(Vargen* gen, int t);
Vargen* VARGEN_get_next(Vargen* gen);
REAL VARGEN_get_P(Vargen* gen, int t);
REAL VARGEN_get_P_ava(Vargen* gen, int t);
REAL VARGEN_get_P_max(Vargen* gen);
REAL VARGEN_get_P_min(Vargen* gen);
REAL VARGEN_get_P_std(Vargen* gen, int t);
REAL VARGEN_get_Q(Vargen* gen, int t);
REAL VARGEN_get_Q_max(Vargen* gen);
REAL VARGEN_get_Q_min(Vargen* gen);
void VARGEN_get_var_values(Vargen* gen, Vec* values, int code);
char* VARGEN_get_var_info_string(Vargen* gen, int index);
int VARGEN_get_num_vars(void* gen, unsigned char var, int t_start, int t_end);
Vec* VARGEN_get_var_indices(void* gen, unsigned char var, int t_start, int t_end);
char* VARGEN_get_json_string(Vargen* gen, char* output);
BOOL VARGEN_has_flags(void* gen, char flag_type, unsigned char mask);
BOOL VARGEN_has_properties(void* gen, char prop);
void VARGEN_init(Vargen* gen, int num_periods);
BOOL VARGEN_is_wind(Vargen* gen);
BOOL VARGEN_is_solar(Vargen* gen);
Vargen* VARGEN_list_add(Vargen* gen_list, Vargen* gen);
Vargen* VARGEN_list_del(Vargen* gen_list, Vargen* gen);
int VARGEN_list_len(Vargen* gen_list);
Vargen* VARGEN_new(int num_periods);
void VARGEN_propagate_data_in_time(Vargen* gen, int start, int end);
void VARGEN_set_name(Vargen* gen, char* name);
void VARGEN_set_type(Vargen* gen, int type);
void VARGEN_set_bus(Vargen* gen, Bus* bus);
void VARGEN_set_index(Vargen* gen, int index);
void VARGEN_set_P(Vargen* gen, REAL P, int t);
void VARGEN_set_P_ava(Vargen* gen, REAL P, int t);
void VARGEN_set_P_max(Vargen* gen, REAL P);
void VARGEN_set_P_min(Vargen* gen, REAL P);
void VARGEN_set_P_std(Vargen* gen, REAL P, int t);
void VARGEN_set_Q(Vargen* gen, REAL Q, int t);
void VARGEN_set_Q_max(Vargen* gen, REAL Q);
void VARGEN_set_Q_min(Vargen* gen, REAL Q);
int VARGEN_set_flags(void* gen, char flag_type, unsigned char mask, int index);
void VARGEN_set_var_values(Vargen* gen, Vec* values);
void VARGEN_show(Vargen* gen, int t);

#endif
