/** @file constr_BAT_DYN.h
 *  @brief This file lists the constants and routines associated with the constraint of type BAT_DYN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_BAT_DYN_HEADER__
#define __CONSTR_BAT_DYN_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
Constr* CONSTR_BAT_DYN_new(Net* net);
void CONSTR_BAT_DYN_init(Constr* c);
void CONSTR_BAT_DYN_count_step(Constr* c, Branch* br, int t);
void CONSTR_BAT_DYN_allocate(Constr* c);
void CONSTR_BAT_DYN_clear(Constr* c);
void CONSTR_BAT_DYN_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_BAT_DYN_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_BAT_DYN_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_BAT_DYN_free(Constr* c);

#endif
