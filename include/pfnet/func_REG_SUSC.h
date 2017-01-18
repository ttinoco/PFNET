/** @file func_REG_SUSC.h
 *  @brief This file lists the constants and routines associated with the function of type REG_SUSC.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_REG_SUSC_HEADER__
#define __FUNC_REG_SUSC_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_REG_SUSC_PARAM 1e-4

// Function prototypes
void FUNC_REG_SUSC_init(Func* f);
void FUNC_REG_SUSC_count_step(Func* f, Branch* br, int t);
void FUNC_REG_SUSC_allocate(Func* f);
void FUNC_REG_SUSC_clear(Func* f);
void FUNC_REG_SUSC_analyze_step(Func* f, Branch* br, int t);
void FUNC_REG_SUSC_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_REG_SUSC_free(Func* f);

#endif
