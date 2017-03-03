/** @file constr_BOUND.h
 *  @brief This file lists the constants and routines associated with the constraint of type BOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_BOUND_HEADER__
#define __CONSTR_BOUND_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_BOUND_PARAM 1e-4

// Function prototypes
void CONSTR_BOUND_init(Constr* c);
void CONSTR_BOUND_count_step(Constr* c, Branch* br, int t);
void CONSTR_BOUND_allocate(Constr* c);
void CONSTR_BOUND_clear(Constr* c);
void CONSTR_BOUND_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_BOUND_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ev);
void CONSTR_BOUND_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_BOUND_free(Constr* c);

#endif
