/** @file constr_AC_FLOW_LIM.h
 *  @brief This file lists the constants and routines associated with the constraint of type AC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_AC_FLOW_LIM_HEADER__
#define __CONSTR_AC_FLOW_LIM_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_AC_FLOW_LIM_PARAM 1e-6

// Function prototypes
Constr* CONSTR_AC_FLOW_LIM_new(Net* net);
void CONSTR_AC_FLOW_LIM_count_step(Constr* c, Bus* bus, int t);
void CONSTR_AC_FLOW_LIM_analyze_step(Constr* c, Bus* bus, int t);
void CONSTR_AC_FLOW_LIM_eval_step(Constr* c, Bus* bus, int t, Vec* v, Vec* ve);
void CONSTR_AC_FLOW_LIM_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
