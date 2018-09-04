/** @file constr_GEN_RAMP.h
 *  @brief This file lists the constants and routines associated with the constraint of type GEN_RAMP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_GEN_RAMP_HEADER__
#define __CONSTR_GEN_RAMP_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_GEN_RAMP_new(Net* net);
void CONSTR_GEN_RAMP_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_GEN_RAMP_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_GEN_RAMP_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_GEN_RAMP_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
#endif
