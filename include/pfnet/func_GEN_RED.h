/** @file func_GEN_RED.h
 *  @brief This file lists the constants and routines associated with the function of type GEN_RED.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2019, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_GEN_RED_HEADER__
#define __FUNC_GEN_RED_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_GEN_RED_new(REAL weight, Net* net);
void FUNC_GEN_RED_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_GEN_RED_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_GEN_RED_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values);

#endif
