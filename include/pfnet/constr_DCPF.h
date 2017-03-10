/** @file constr_DCPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_DCPF_HEADER__
#define __CONSTR_DCPF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
void CONSTR_DCPF_init(Constr* c);
void CONSTR_DCPF_count_step(Constr* c, Branch* br, int t);
void CONSTR_DCPF_allocate(Constr* c);
void CONSTR_DCPF_clear(Constr* c);
void CONSTR_DCPF_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_DCPF_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_DCPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_DCPF_free(Constr* c);

#endif
