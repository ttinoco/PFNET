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
/** @} */

// Variable generator types
#define VARGEN_TYPE_WIND 0       /**< @brief Type: wind farm */
#define VARGEN_TYPE_SOLAR 1       /**< @brief Type: solar plant */

// Properties
/** \defgroup vargen_props Variable Generator Property Masks 
 *  @{
 */
#define VARGEN_PROP_ANY 0x00       /**< @brief Property: any */
/** @} */

// Variable generator
typedef struct Vargen Vargen;

// Prototypes
void* VARGEN_array_get(void* gen, int index);
Vargen* VARGEN_array_new(int num);
void VARGEN_array_show(Vargen* gen, int num);
void VARGEN_clear_flags(Vargen* gen, char flag_type);
void* VARGEN_get_bus(Vargen* gen);
int VARGEN_get_index(Vargen* gen);
int VARGEN_get_index_P(Vargen* gen);
Vargen* VARGEN_get_next(Vargen* gen);
REAL VARGEN_get_P(Vargen* gen);
REAL VARGEN_get_P_max(Vargen* gen);
REAL VARGEN_get_P_std(Vargen* gen);
void VARGEN_get_var_values(Vargen* gen, Vec* values, int code);
int VARGEN_get_var_index(void* gen, char var);
BOOL VARGEN_has_flags(void* gen, char flag_type, char mask);
BOOL VARGEN_has_properties(void* gen, char prop);
void VARGEN_init(Vargen* gen);
BOOL VARGEN_is_wind_farm(Vargen* gen);
BOOL VARGEN_is_solar_plant(Vargen* gen);
Vargen* VARGEN_list_add(Vargen* gen_list, Vargen* gen);
int VARGEN_list_len(Vargen* gen_list);
Vargen* VARGEN_new(void);
void VARGEN_set_type(Vargen* gen, int type);
void VARGEN_set_bus(Vargen* gen, void* bus);
void VARGEN_set_index(Vargen* gen, int index);
void VARGEN_set_P(Vargen* gen, REAL P);
void VARGEN_set_P_max(Vargen* gen, REAL P_max);
void VARGEN_set_P_std(Vargen* gen, REAL P_std);
int VARGEN_set_flags(void* gen, char flag_type, char mask, int index);
void VARGEN_set_var_values(Vargen* gen, Vec* values);
void VARGEN_show(Vargen* gen);

#endif
