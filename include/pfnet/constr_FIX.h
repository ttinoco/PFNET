/** @file constr_FIX.h
 *  @brief This file lists the constants and routines associated with the constraint of type FIX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_FIX_HEADER__
#define __CONSTR_FIX_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_FIX_new(Net* net);
void CONSTR_FIX_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_FIX_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_FIX_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_FIX_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
