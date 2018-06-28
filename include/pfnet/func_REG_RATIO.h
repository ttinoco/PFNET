/** @file func_REG_RATIO.h
 *  @brief This file lists the constants and routines associated with the function of type REG_RATIO.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_RATIO_HEADER__
#define __FUNC_REG_RATIO_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_RATIO_PARAM 1e-4

// Function prototypes
Func* FUNC_REG_RATIO_new(REAL weight, Net* net);
void FUNC_REG_RATIO_count_step(Func* f, Branch* br, int t);
void FUNC_REG_RATIO_analyze_step(Func* f, Branch* br, int t);
void FUNC_REG_RATIO_eval_step(Func* f, Branch* br, int t, Vec* var_values);

#endif
