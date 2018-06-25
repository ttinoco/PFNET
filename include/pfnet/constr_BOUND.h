/** @file constr_BOUND.h
 *  @brief This file lists the constants and routines associated with the constraint of type BOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_BOUND_HEADER__
#define __CONSTR_BOUND_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_BOUND_new(Net* net);
void CONSTR_BOUND_count_step(Constr* c, Branch* br, int t);
void CONSTR_BOUND_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_BOUND_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_BOUND_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
