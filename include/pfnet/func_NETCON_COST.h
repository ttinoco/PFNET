/** @file func_NETCON_COST.h
 *  @brief This file lists the constants and routines associated with the function of type NETCON_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_NETCON_COST_HEADER__
#define __FUNC_NETCON_COST_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
void FUNC_NETCON_COST_init(Func* f);
void FUNC_NETCON_COST_count_step(Func* f, Branch* br, int t);
void FUNC_NETCON_COST_allocate(Func* f);
void FUNC_NETCON_COST_clear(Func* f);
void FUNC_NETCON_COST_analyze_step(Func* f, Branch* br, int t);
void FUNC_NETCON_COST_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_NETCON_COST_free(Func* f);

#endif
