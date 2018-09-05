/** @file func_FACTS_QSET.h
 *  @brief This file lists the constants and routines associated with the function of type FACTS_QSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_FACTS_QSET_HEADER__
#define __FUNC_FACTS_QSET_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
Func* FUNC_FACTS_QSET_new(REAL weight, Net* net);
void FUNC_FACTS_QSET_count_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_FACTS_QSET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t);
void FUNC_FACTS_QSET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* v);

#endif
