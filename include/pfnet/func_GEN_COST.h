/** @file func_GEN_COST.h
 *  @brief This file lists the constants and routines associated with the function of type GEN_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_GEN_COST_HEADER__
#define __FUNC_GEN_COST_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
void FUNC_GEN_COST_init(Func* f);
void FUNC_GEN_COST_count_branch(Func* f, Branch *branch);
void FUNC_GEN_COST_allocate(Func* f);
void FUNC_GEN_COST_clear(Func* f);
void FUNC_GEN_COST_analyze_branch(Func* f, Branch *branch);
void FUNC_GEN_COST_eval_branch(Func* f, Branch* branch, Vec* var_values);
void FUNC_GEN_COST_free(Func* f);

#endif
