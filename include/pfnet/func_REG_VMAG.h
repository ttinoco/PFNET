/** @file func_REG_VMAG.h
 *  @brief This file lists the constants and routines associated with the function of type REG_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_VMAG_HEADER__
#define __FUNC_REG_VMAG_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_VMAG_PARAM 0.2

// Function prototypes
Func* FUNC_REG_VMAG_new(REAL weight, Net* net);
void FUNC_REG_VMAG_init(Func* f);
void FUNC_REG_VMAG_count_step(Func* f, Branch* br, int t);
void FUNC_REG_VMAG_allocate(Func* f);
void FUNC_REG_VMAG_clear(Func* f);
void FUNC_REG_VMAG_analyze_step(Func* f, Branch* br, int t);
void FUNC_REG_VMAG_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_REG_VMAG_free(Func* f);

#endif
