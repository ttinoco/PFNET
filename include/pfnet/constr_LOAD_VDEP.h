/** @file constr_LOAD_VDEP.h
 *  @brief This file lists the constants and routines associated with the constraint of type LOAD_VDEP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_LOAD_VDEP_HEADER__
#define __CONSTR_LOAD_VDEP_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_LOAD_VDEP_new(Net* net);
void CONSTR_LOAD_VDEP_count_step(Constr* c, Bus* bus, int t);
void CONSTR_LOAD_VDEP_analyze_step(Constr* c, Bus* bus, int t);
void CONSTR_LOAD_VDEP_eval_step(Constr* c, Bus* bus, int t, Vec* v, Vec* ve);
void CONSTR_LOAD_VDEP_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
