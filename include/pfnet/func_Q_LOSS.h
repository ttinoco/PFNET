/** @file func_Q_LOSS.h
 *  @brief This file lists the constants and routines associated with the function of type Q_LOSS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_Q_LOSS_HEADER__
#define __FUNC_Q_LOSS_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_Q_LOSS_new(REAL weight, Net* net);
void FUNC_Q_LOSS_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_Q_LOSS_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_Q_LOSS_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values);

#endif
