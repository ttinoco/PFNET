/** @file func_REG_PHASE.h
 *  @brief This file lists the constants and routines associated with the function of type REG_PHASE.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_PHASE_HEADER__
#define __FUNC_REG_PHASE_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_PHASE_PARAM 1e-4

// Function prototypes
void FUNC_REG_PHASE_init(Func* f);
void FUNC_REG_PHASE_count_branch(Func* f, Branch *branch);
void FUNC_REG_PHASE_allocate(Func* f);
void FUNC_REG_PHASE_clear(Func* f);
void FUNC_REG_PHASE_analyze_branch(Func* f, Branch *branch);
void FUNC_REG_PHASE_eval_branch(Func* f, Branch* branch, Vec* var_values);
void FUNC_REG_PHASE_free(Func* f);

#endif
