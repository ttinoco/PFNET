/** @file constr_GEN_RAMP.h
 *  @brief This file lists the constants and routines associated with the constraint of type GEN_RAMP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_GEN_RAMP_HEADER__
#define __CONSTR_GEN_RAMP_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
void CONSTR_GEN_RAMP_init(Constr* c);
void CONSTR_GEN_RAMP_count_step(Constr* c, Branch* br, int t);
void CONSTR_GEN_RAMP_allocate(Constr* c);
void CONSTR_GEN_RAMP_clear(Constr* c);
void CONSTR_GEN_RAMP_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_GEN_RAMP_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ev);
void CONSTR_GEN_RAMP_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_GEN_RAMP_free(Constr* c);

#endif
