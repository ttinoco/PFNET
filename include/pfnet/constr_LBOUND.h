/** @file constr_LBOUND.h
 *  @brief This file lists the constants and routines associated with the constraint of type LBOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_LBOUND_HEADER__
#define __CONSTR_LBOUND_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_LBOUND_new(Net* net);
void CONSTR_LBOUND_init(Constr* c);
void CONSTR_LBOUND_count_step(Constr* c, Branch* br, int t);
void CONSTR_LBOUND_allocate(Constr* c);
void CONSTR_LBOUND_clear(Constr* c);
void CONSTR_LBOUND_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_LBOUND_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_LBOUND_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_LBOUND_free(Constr* c);

#endif
