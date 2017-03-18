/** @file constr_LINPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type LINPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_LINPF_HEADER__
#define __CONSTR_LINPF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_LINPF_new(Net* net);
void CONSTR_LINPF_init(Constr* c);
void CONSTR_LINPF_count_step(Constr* c, Branch* br, int t);
void CONSTR_LINPF_allocate(Constr* c);
void CONSTR_LINPF_clear(Constr* c);
void CONSTR_LINPF_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_LINPF_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_LINPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_LINPF_free(Constr* c);

#endif
