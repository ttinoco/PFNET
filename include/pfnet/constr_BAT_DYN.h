/** @file constr_BAT_DYN.h
 *  @brief This file lists the constants and routines associated with the constraint of type BAT_DYN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_BAT_DYN_HEADER__
#define __CONSTR_BAT_DYN_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_BAT_DYN_new(Net* net);
void CONSTR_BAT_DYN_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_BAT_DYN_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_BAT_DYN_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_BAT_DYN_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
