/** @file constr_REG_VSET.h
 *  @brief This file lists the constants and routines associated with the constraint of type REG_VSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_REG_VSET_HEADER__
#define __CONSTR_REG_VSET_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_REG_VSET_PARAM 1e-8
#define CONSTR_REG_VSET_MAX_YZ 1e8

// Function prototypes
Constr* CONSTR_REG_VSET_new(Net* net);
void CONSTR_REG_VSET_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_REG_VSET_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_REG_VSET_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_REG_VSET_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);

#endif
