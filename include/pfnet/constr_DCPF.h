/** @file constr_DCPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_DCPF_HEADER__
#define __CONSTR_DCPF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_DCPF_new(Net* net);
void CONSTR_DCPF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_DCPF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_DCPF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_DCPF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
