/** @file func_SLIM_VMAG.h
 *  @brief This file lists the constants and routines associated with the function of type SLIM_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_SLIM_VMAG_HEADER__
#define __FUNC_SLIM_VMAG_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_SLIM_VMAG_PARAM 1e-4

// Function prototypes
Func* FUNC_SLIM_VMAG_new(REAL weight, Net* net);
void FUNC_SLIM_VMAG_count_step(Func* f, Branch* br, int t);
void FUNC_SLIM_VMAG_analyze_step(Func* f, Branch* br, int t);
void FUNC_SLIM_VMAG_eval_step(Func* f, Branch* br, int t, Vec* var_values);

#endif
