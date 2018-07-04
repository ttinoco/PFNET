/** @file constr_ACPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type ACPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_ACPF_HEADER__
#define __CONSTR_ACPF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_ACPF_new(Net* net);
void CONSTR_ACPF_count_step(Constr* c, Bus* bus, int t);
void CONSTR_ACPF_analyze_step(Constr* c, Bus* bus, int t);
void CONSTR_ACPF_eval_step(Constr* c, Bus* bus, int t, Vec* v, Vec* ve);
void CONSTR_ACPF_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
