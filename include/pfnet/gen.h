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

// Properties
/** \defgroup gen_props Generator Property Masks 
 *  @{
 */
#define GEN_PROP_ANY 0x00       /**< @brief Property: any */
#define GEN_PROP_SLACK 0x01     /**< @brief Property: slack generator */
#define GEN_PROP_REG 0x02       /**< @brief Property: regulating generator */
#define GEN_PROP_NOT_REG 0x04   /**< @brief Property: non-regulating generator */
#define GEN_PROP_NOT_SLACK 0x08 /**< @brief Property: non-slack generator */
/** @} */

// Generator
typedef struct Gen Gen;

// Prototypes
void* GEN_array_get(void* gen, int index);
Gen* GEN_array_new(int num);
void GEN_array_show(Gen* gen, int num);
void GEN_clear_flags(Gen* gen, char flag_type);
void* GEN_get_bus(Gen* gen);
void* GEN_get_reg_bus(Gen* gen);
REAL GEN_get_cost_coeff_Q0(Gen* gen);
REAL GEN_get_cost_coeff_Q1(Gen* gen);
REAL GEN_get_cost_coeff_Q2(Gen* gen);
int GEN_get_index(Gen* gen);
int GEN_get_index_P(Gen* gen);
int GEN_get_index_Q(Gen* gen);
Gen* GEN_get_next(Gen* gen);
Gen* GEN_get_reg_next(Gen* gen);
REAL GEN_get_P(Gen* gen);
REAL GEN_get_Q(Gen* gen);
REAL GEN_get_P_max(Gen* gen);
REAL GEN_get_P_min(Gen* gen);
REAL GEN_get_Q_max(Gen* gen);
REAL GEN_get_Q_min(Gen* gen);
void GEN_get_var_values(Gen* gen, Vec* values);
int GEN_get_var_index(void* gen, char var);
BOOL GEN_has_flags(void* gen, char flag_type, char mask);
BOOL GEN_has_properties(void* gen, char prop);
void GEN_init(Gen* gen);
BOOL GEN_is_regulator(Gen* gen);
BOOL GEN_is_slack(Gen* gen);
Gen* GEN_list_add(Gen* gen_list, Gen* gen);
int GEN_list_len(Gen* gen_list);
Gen* GEN_list_reg_add(Gen* reg_gen_list, Gen* reg_gen);
int GEN_list_reg_len(Gen* reg_gen_list);
Gen* GEN_new(void);
void GEN_set_bus(Gen* gen, void* bus);
void GEN_set_reg_bus(Gen* gen, void* reg_bus);
void GEN_set_regulator(Gen* gen, BOOL regulator);
void GEN_set_index(Gen* gen, int index);
void GEN_set_P(Gen* gen, REAL P);
void GEN_set_P_max(Gen* gen, REAL P_max);
void GEN_set_P_min(Gen* gen, REAL P_min);
void GEN_set_Q(Gen* gen, REAL Q);
void GEN_set_Q_max(Gen* gen, REAL Q_max);
void GEN_set_Q_min(Gen* gen, REAL Q_min);
int GEN_set_flags(void* gen, char flag_type, char mask, int index);
void GEN_set_var_values(Gen* gen, Vec* values);
void GEN_show(Gen* gen);

#endif
