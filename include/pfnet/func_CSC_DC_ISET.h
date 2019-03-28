/** @file func_CSC_ISET.h
 *  @brief This file lists the constants and routines associated with the function of type CSC_ISET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2019, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_CSC_DC_ISET_HEADER__
#define __FUNC_CSC_DC_ISET_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_CSC_DC_ISET_new(REAL weight, Net* net);
void FUNC_CSC_DC_ISET_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_CSC_DC_ISET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_CSC_DC_ISET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* v);

#endif
