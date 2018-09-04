/** @file func_REG_VAR.h
 *  @brief This file lists the constants and routines associated with the function of type REG_VAR.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_VAR_HEADER__
#define __FUNC_REG_VAR_HEADER__

#include <math.h>
#include "func.h"

// Data
typedef struct Func_REG_VAR_Data Func_REG_VAR_Data;

// Function prototypes
Func* FUNC_REG_VAR_new(REAL weight, Net* net);
void FUNC_REG_VAR_init(Func* f);
void FUNC_REG_VAR_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_REG_VAR_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_REG_VAR_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values);
void FUNC_REG_VAR_free(Func* f);
void FUNC_REG_VAR_set_parameter(Func* f, char* key, void* value);

#endif
