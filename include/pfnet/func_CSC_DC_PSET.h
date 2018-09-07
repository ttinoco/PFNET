/** @file func_CSC_PSET.h
 *  @brief This file lists the constants and routines associated with the function of type CSC_PSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_CSC_DC_PSET_HEADER__
#define __FUNC_CSC_DC_PSET_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_CSC_DC_PSET_new(REAL weight, Net* net);
void FUNC_CSC_DC_PSET_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_CSC_DC_PSET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_CSC_DC_PSET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* v);

#endif
