/** @file constr_NBOUND.h
 *  @brief This file lists the constants and routines associated with the constraint of type NBOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_NBOUND_HEADER__
#define __CONSTR_NBOUND_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_NBOUND_PARAM 1e-4

// Function prototypes
Constr* CONSTR_NBOUND_new(Net* net);
void CONSTR_NBOUND_init(Constr* c);
void CONSTR_NBOUND_count_step(Constr* c, Branch* br, int t);
void CONSTR_NBOUND_allocate(Constr* c);
void CONSTR_NBOUND_clear(Constr* c);
void CONSTR_NBOUND_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_NBOUND_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_NBOUND_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_NBOUND_free(Constr* c);

#endif
