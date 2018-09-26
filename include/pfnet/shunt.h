/** @file shunt.h
 *  @brief This file list the constants and routines associated with the Shunt data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __SHUNT_HEADER__
#define __SHUNT_HEADER__

#include "stdio.h"
#include "types.h"
#include "vector.h"
#include "list.h"

// Shunt types
#define SHUNT_TYPE_FIXED 0       /**< @brief Type: Fixed shunt */
#define SHUNT_TYPE_SWITCHED 1    /**< @brief Type: Switched shunt that is locked */
#define SHUNT_TYPE_SWITCHED_V 2  /**< @brief Type: Switched shunt that provides voltage regulation*/

// Switched shunt modes
#define SHUNT_MODE_CONT 0 /**< @brief Mode: Continuous */
#define SHUNT_MODE_DIS 1  /**< @brief Mode: Discrete */

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
#define SHUNT_BUFFER_SIZE 100    /**< @brief Constant: buffer size for strings */
#define SHUNT_NUM_JSON_FIELDS 15 /**< @brief Constant: max number of json fields */
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
void SHUNT_clear_sensitivities(Shunt* shunt);
void SHUNT_clear_flags(Shunt* shunt, char flag_type);
void SHUNT_copy_from_shunt(Shunt* shunt, Shunt* other);

char SHUNT_get_flags_vars(Shunt* shunt);
char SHUNT_get_flags_fixed(Shunt* shunt);
char SHUNT_get_flags_bounded(Shunt* shunt);
char SHUNT_get_flags_sparse(Shunt* shunt);

REAL SHUNT_get_sens_b_u_bound(Shunt* shunt, int t);
REAL* SHUNT_get_sens_b_u_bound_array(Shunt* shunt);
REAL SHUNT_get_sens_b_l_bound(Shunt* shunt, int t);
REAL* SHUNT_get_sens_b_l_bound_array(Shunt* shunt);

char SHUNT_get_type(Shunt* shunt);
char SHUNT_get_mode(Shunt* shunt);
char* SHUNT_get_name(Shunt* shunt);
int SHUNT_get_num_periods(Shunt* shunt);
char SHUNT_get_obj_type(void* shunt);
int SHUNT_get_index(Shunt* shunt);
int SHUNT_get_index_b(Shunt* shunt, int t);
int* SHUNT_get_index_b_array(Shunt* shunt);
Bus* SHUNT_get_bus(Shunt* shunt);
Bus* SHUNT_get_reg_bus(Shunt* shunt);
REAL SHUNT_get_g(Shunt* shunt);
REAL SHUNT_get_b(Shunt* shunt, int t);
REAL* SHUNT_get_b_array(Shunt* shunt);
REAL SHUNT_get_b_max(Shunt* shunt);
REAL SHUNT_get_b_min(Shunt* shunt);
REAL* SHUNT_get_b_values(Shunt* shunt);
int SHUNT_get_num_b_values(Shunt* shunt);
Shunt* SHUNT_get_next(Shunt* shunt);
Shunt* SHUNT_get_reg_next(Shunt* shunt);
void SHUNT_get_var_values(Shunt* shunt, Vec* values, int code);
char* SHUNT_get_var_info_string(Shunt* shunt, int index);
int SHUNT_get_num_vars(void* shunt, unsigned char var, int t_start, int t_end);
Vec* SHUNT_get_var_indices(void* shunt, unsigned char var, int t_start, int t_end);
char* SHUNT_get_json_string(Shunt* shunt, char* output);
BOOL SHUNT_has_flags(void* shunt, char flag_type, unsigned char mask);
BOOL SHUNT_has_properties(void* shunt, char prop);
void SHUNT_init(Shunt* shunt, int num_periods);
BOOL SHUNT_is_equal(Shunt* shunt, Shunt* other);
BOOL SHUNT_is_fixed(Shunt* shunt);
BOOL SHUNT_is_switched(Shunt* shunt);
BOOL SHUNT_is_switched_v(Shunt* shunt);
BOOL SHUNT_is_switched_locked(Shunt* shunt);
BOOL SHUNT_is_continuous(Shunt* shunt);
BOOL SHUNT_is_discrete(Shunt* shunt);
Shunt* SHUNT_list_add(Shunt* shunt_list, Shunt* shunt);
Shunt* SHUNT_list_del(Shunt* shunt_list, Shunt* shunt);
int SHUNT_list_len(Shunt* shunt_list);
Shunt* SHUNT_list_reg_add(Shunt* reg_shunt_list, Shunt* reg_shunt);
Shunt* SHUNT_list_reg_del(Shunt* reg_shunt_list, Shunt* reg_shunt);
int SHUNT_list_reg_len(Shunt* reg_shunt_list);
Shunt* SHUNT_new(int num_periods);
void SHUNT_propagate_data_in_time(Shunt* shunt, int start, int end);
void SHUNT_round_b(Shunt* shunt, int t);
void SHUNT_set_sens_b_u_bound(Shunt* shunt, REAL value, int t);
void SHUNT_set_sens_b_l_bound(Shunt* shunt, REAL value, int t);
void SHUNT_set_type(Shunt* shunt, char type);
void SHUNT_set_mode(Shunt* shunt, char mode);
void SHUNT_set_name(Shunt* shunt, char* name);
void SHUNT_set_bus(Shunt* shunt, Bus* bus);
void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus);
void SHUNT_set_index(Shunt* shunt, int index);
void SHUNT_set_g(Shunt* shunt, REAL g);
void SHUNT_set_b(Shunt* shunt, REAL b, int t);
void SHUNT_set_b_max(Shunt* shunt, REAL b_max);
void SHUNT_set_b_min(Shunt* shunt, REAL b_min);
void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num);
int SHUNT_set_flags(void* shunt, char flag_type, unsigned char mask, int index);
void SHUNT_set_var_values(Shunt* shunt, Vec* values);
void SHUNT_show(Shunt* shunt, int t);

#endif
