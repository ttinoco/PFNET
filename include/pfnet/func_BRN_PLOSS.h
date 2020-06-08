/** @file func_BRN_PLOSS.h
 *  @brief This file lists the constants and routines associated with the function of type BRN_PLOSS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_BRN_PLOSS_HEADER__
#define __FUNC_BRN_PLOSS_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_BRN_PLOSS_new(REAL weight, Net* net);
void FUNC_BRN_PLOSS_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_BRN_PLOSS_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_BRN_PLOSS_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values);

#endif
