/** @file func_SP_CONTROLS.h
 *  @brief This file lists the constants and routines associated with the function of type SP_CONTROLS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_SP_CONTROLS_HEADER__
#define __FUNC_SP_CONTROLS_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_SP_CONTROLS_EPS 1e-6   /**< @brief Parameter for making function smooth. */
#define FUNC_SP_CONTROLS_CEPS 1e-4  /**< @brief Parameter for avoiding small control ranges. */

// Function prototypes
Func* FUNC_SP_CONTROLS_new(REAL weight, Net* net);
void FUNC_SP_CONTROLS_init(Func* f);
void FUNC_SP_CONTROLS_count_step(Func* f, Branch* br, int t);
void FUNC_SP_CONTROLS_allocate(Func* f);
void FUNC_SP_CONTROLS_clear(Func* f);
void FUNC_SP_CONTROLS_analyze_step(Func* f, Branch* br, int t);
void FUNC_SP_CONTROLS_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_SP_CONTROLS_free(Func* f);

#endif
