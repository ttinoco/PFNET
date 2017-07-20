/** @file shunt.h
 *  @brief This file list the constants and routines associated with the Shunt data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __SHUNT_HEADER__
#define __SHUNT_HEADER__

#include "stdio.h"
#include "types.h"
#include "vector.h"
#include "list.h"

// Variables
/** \defgroup shunt_vars Shunt Variable Masks
 *  @{
 */
#define SHUNT_VAR_SUSC 0x01     /**< @brief Variable: susceptance **/
/** @} */

// Infinity
#define SHUNT_INF_SUSC 1e8 /**< @brief Infinite susceptance */

// Properties
/** \defgroup shunt_props Shunt Property Masks
 *  @{
 */
#define SHUNT_PROP_ANY 0x00        /**< @brief Property: any **/
#define SHUNT_PROP_SWITCHED_V 0x01 /**< @brief Property: switched that regulates bus voltage **/
/** @} */

// Constants
/** \defgroup shunt_const Shunt Constants
 *  @{
 */
#define SHUNT_BUFFER_SIZE 100  /**< @brief Constant: buffer size for strings */
/** @} */

// Struct
typedef struct Shunt Shunt;

// Other
typedef struct Bus Bus;

// Function prototypes
void* SHUNT_array_get(void* shunt_array, int index);
void SHUNT_array_del(Shunt* shunt_array, int size);
Shunt* SHUNT_array_new(int size, int num_periods);
void SHUNT_array_show(Shunt* shunt_array, int size, int t);
void SHUNT_clear_flags(Shunt* shunt, char flag_type);
void SHUNT_propagate_data_in_time(Shunt* shunt);
int SHUNT_get_num_periods(Shunt* shunt);
char SHUNT_get_obj_type(void* shunt);
int SHUNT_get_index(Shunt* shunt);
int SHUNT_get_index_b(Shunt* shunt, int t);
Bus* SHUNT_get_bus(Shunt* shunt);
Bus* SHUNT_get_reg_bus(Shunt* shunt);
REAL SHUNT_get_g(Shunt* shunt);
REAL SHUNT_get_b(Shunt* shunt, int t);
REAL SHUNT_get_b_max(Shunt* shunt);
REAL SHUNT_get_b_min(Shunt* shunt);
Shunt* SHUNT_get_next(Shunt* shunt);
Shunt* SHUNT_get_reg_next(Shunt* shunt);
void SHUNT_get_var_values(Shunt* shunt, Vec* values, int code);
int SHUNT_get_num_vars(void* shunt, unsigned char var, int t_start, int t_end);
Vec* SHUNT_get_var_indices(void* shunt, unsigned char var, int t_start, int t_end);
char* SHUNT_get_json_string(Shunt* shunt);
BOOL SHUNT_has_flags(void* shunt, char flag_type, unsigned char mask);
BOOL SHUNT_has_properties(void* shunt, char prop);
void SHUNT_init(Shunt* shunt, int num_periods);
BOOL SHUNT_is_fixed(Shunt* shunt);
BOOL SHUNT_is_switched(Shunt* shunt);
BOOL SHUNT_is_switched_v(Shunt* shunt);
Shunt* SHUNT_list_add(Shunt *shunt_list, Shunt* shunt);
int SHUNT_list_len(Shunt* shunt_list);
Shunt* SHUNT_list_reg_add(Shunt *reg_shunt_list, Shunt* reg_shunt);
int SHUNT_list_reg_len(Shunt* reg_shunt_list);
Shunt* SHUNT_new(int num_periods);
void SHUNT_set_bus(Shunt* shunt, Bus* bus);
void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus);
void SHUNT_set_index(Shunt* shunt, int index);
void SHUNT_set_g(Shunt* shunt, REAL g);
void SHUNT_set_b(Shunt* shunt, REAL b, int t);
void SHUNT_set_b_max(Shunt* shunt, REAL b_max);
void SHUNT_set_b_min(Shunt* shunt, REAL b_min);
void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num, REAL norm);
int SHUNT_set_flags(void* shunt, char flag_type, unsigned char mask, int index);
void SHUNT_set_var_values(Shunt* shunt, Vec* values);
void SHUNT_show(Shunt* shunt, int t);

#endif
