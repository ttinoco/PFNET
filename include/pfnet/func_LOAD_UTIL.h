/** @file func_LOAD_UTIL.h
 *  @brief This file lists the constants and routines associated with the function of type LOAD_UTIL.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_LOAD_UTIL_HEADER__
#define __FUNC_LOAD_UTIL_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_LOAD_UTIL_new(REAL weight, Net* net);
void FUNC_LOAD_UTIL_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_LOAD_UTIL_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_LOAD_UTIL_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values);

#endif
