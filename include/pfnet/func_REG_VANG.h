/** @file func_REG_VANG.h
 *  @brief This file lists the constants and routines associated with the function of type REG_VANG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_VANG_HEADER__
#define __FUNC_REG_VANG_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_VANG_PARAM 3.1416

// Function prototypes
Func* FUNC_REG_VANG_new(REAL weight, Net* net);
void FUNC_REG_VANG_count_step(Func* f, Bus* bus, int t);
void FUNC_REG_VANG_analyze_step(Func* f, Bus* bus, int t);
void FUNC_REG_VANG_eval_step(Func* f, Bus* bus, int t, Vec* var_values);

#endif
