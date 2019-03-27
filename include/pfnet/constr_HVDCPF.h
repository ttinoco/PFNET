/** @file constr_HVDCPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type HVDCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_HVDCPF_HEADER__
#define __CONSTR_HVDCPF_HEADER__

#include <math.h>
#include "constr.h"

#define CONSTR_HVDCPF_MINR 1e-5 // p.u.

// Function prototypes
Constr* CONSTR_HVDCPF_new(Net* net);
void CONSTR_HVDCPF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_HVDCPF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_HVDCPF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_HVDCPF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
