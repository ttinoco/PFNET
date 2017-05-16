/** @file constr_LOAD_PF.h
 *  @brief This file lists the constants and routines associated with the constraint of type LOAD_PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_LOAD_PF_HEADER__
#define __CONSTR_LOAD_PF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_LOAD_PF_new(Net* net);
void CONSTR_LOAD_PF_init(Constr* c);
void CONSTR_LOAD_PF_count_step(Constr* c, Branch* br, int t);
void CONSTR_LOAD_PF_allocate(Constr* c);
void CONSTR_LOAD_PF_clear(Constr* c);
void CONSTR_LOAD_PF_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_LOAD_PF_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_LOAD_PF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_LOAD_PF_free(Constr* c);

#endif
