/** @file func_REG_PHASE.h
 *  @brief This file lists the constants and routines associated with the function of type REG_PHASE.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_PHASE_HEADER__
#define __FUNC_REG_PHASE_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_PHASE_PARAM 1e-4

// Function prototypes
Func* FUNC_REG_PHASE_new(REAL weight, Net* net);
void FUNC_REG_PHASE_init(Func* f);
void FUNC_REG_PHASE_count_step(Func* f, Branch* br, int t);
void FUNC_REG_PHASE_allocate(Func* f);
void FUNC_REG_PHASE_clear(Func* f);
void FUNC_REG_PHASE_analyze_step(Func* f, Branch* br, int t);
void FUNC_REG_PHASE_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_REG_PHASE_free(Func* f);

#endif
