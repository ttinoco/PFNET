/** @file bat.h
 *  @brief This file lists the constants and routines associated with the Bat data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
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
#define BAT_VAR_P 0x01    /**< @brief Variable: battery charing power */
#define BAT_VAR_E 0x02    /**< @brief Variable: battery energy level */
/** @} */

// Infinity
#define BAT_INF_P 1000. /**< @brief Infinite charging power */
#define BAT_INF_E 1000. /**< @brief Infinite energy level */

// Properties
/** \defgroup bat_props Battery Property Masks 
 *  @{
 */
#define BAT_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Battery
typedef struct Bat Bat;

// Others
typedef struct Bus Bus;

void* BAT_array_get(void* bat, int index);
Bat* BAT_array_new(int num);
void BAT_array_show(Bat* bat, int num);
void BAT_clear_flags(Bat* bat, char flag_type);
char BAT_get_obj_type(void* bat);
Bus* BAT_get_bus(Bat* bat);
int BAT_get_index(Bat* bat);
int BAT_get_index_P(Bat* bat);
int BAT_get_index_E(Bat* bat);
Bat* BAT_get_next(Bat* bat);
REAL BAT_get_P(Bat* bat);
REAL BAT_get_P_max(Bat* bat);
REAL BAT_get_P_min(Bat* bat);
REAL BAT_get_E(Bat* bat);
REAL BAT_get_E_max(Bat* bat);
void BAT_get_var_values(Bat* bat, Vec* values, int code);
int BAT_get_var_index(void* bat, char var);
BOOL BAT_has_flags(void* bat, char flag_type, char mask);
BOOL BAT_has_properties(void* bat, char prop);
void BAT_init(Bat* bat);
Bat* BAT_list_add(Bat *bat_list, Bat* bat);
int BAT_list_len(Bat* bat_list);
Bat* BAT_new(void);
void BAT_set_bus(Bat* bat, Bus* bus);
void BAT_set_index(Bat* bat, int index);
void BAT_set_P(Bat* bat, REAL P);
void BAT_set_P_max(Bat* bat, REAL P);
void BAT_set_P_min(Bat* bat, REAL P);
void BAT_set_E(Bat* bat, REAL E);
void BAT_set_E_max(Bat* bat, REAL E);
int BAT_set_flags(void* bat, char flag_type, char mask, int index);
void BAT_set_var_values(Bat* bat, Vec* values);
void BAT_show(Bat* bat);

#endif
