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
#define SHUNT_VAR_SUSC_DEV 0x02 /**< @brief Variable: suscptance pos/neg deviation from current value */
/** @} */

// Infinity
#define SHUNT_INF_SUSC 1000. /**< @brief Infinite susceptance */

// Properties
/** \defgroup shunt_props Shunt Property Masks
 *  @{
 */
#define SHUNT_PROP_ANY 0x00        /**< @brief Property: any **/
#define SHUNT_PROP_SWITCHED_V 0x01 /**< @brief Property: switched that regulates bus voltage **/
/** @} */

// Struct
typedef struct Shunt Shunt;

// Other
typedef struct Bus Bus;

// Function prototypes
void* SHUNT_array_get(void* shunt, int index);
void SHUNT_array_del(Shunt* shunt, int num);
Shunt* SHUNT_array_new(int num);
void SHUNT_array_show(Shunt* shunt, int num);
void SHUNT_clear_flags(Shunt* shunt, char flag_type);
char SHUNT_get_obj_type(void* shunt);
int SHUNT_get_index(Shunt* shunt);
int SHUNT_get_index_b(Shunt* shunt);
int SHUNT_get_index_y(Shunt* shunt);
int SHUNT_get_index_z(Shunt* shunt);
Bus* SHUNT_get_bus(Shunt* shunt);
Bus* SHUNT_get_reg_bus(Shunt* shunt);
REAL SHUNT_get_g(Shunt* shunt);
REAL SHUNT_get_b(Shunt* shunt, int t);
REAL SHUNT_get_b_max(Shunt* shunt);
REAL SHUNT_get_b_min(Shunt* shunt);
Shunt* SHUNT_get_next(Shunt* shunt);
Shunt* SHUNT_get_reg_next(Shunt* shunt);
void SHUNT_get_var_values(Shunt* shunt, Vec* values, int code);
Vec* SHUNT_get_var_indices(void* shunt, char var);
BOOL SHUNT_has_flags(void* shunt, char flag_type, char mask);
BOOL SHUNT_has_properties(void* shunt, char prop);
void SHUNT_init(Shunt* shunt);
BOOL SHUNT_is_fixed(Shunt* shunt);
BOOL SHUNT_is_switched(Shunt* shunt);
BOOL SHUNT_is_switched_v(Shunt* shunt);
Shunt* SHUNT_list_add(Shunt *shunt_list, Shunt* shunt);
int SHUNT_list_len(Shunt* shunt_list);
Shunt* SHUNT_list_reg_add(Shunt *reg_shunt_list, Shunt* reg_shunt);
int SHUNT_list_reg_len(Shunt* reg_shunt_list);
Shunt* SHUNT_new(void);
void SHUNT_set_bus(Shunt* shunt, Bus* bus);
void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus);
void SHUNT_set_index(Shunt* shunt, int index);
void SHUNT_set_g(Shunt* shunt, REAL g);
void SHUNT_set_b(Shunt* shunt, REAL b);
void SHUNT_set_b_max(Shunt* shunt, REAL b_max);
void SHUNT_set_b_min(Shunt* shunt, REAL b_min);
void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num, REAL norm);
int SHUNT_set_flags(void* shunt, char flag_type, char mask, int index);
void SHUNT_set_var_values(Shunt* shunt, Vec* values);
void SHUNT_show(Shunt* shunt);

#endif
