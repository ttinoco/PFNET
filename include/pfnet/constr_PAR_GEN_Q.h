/** @file constr_PAR_GEN_Q.h
 *  @brief This file lists the constants and routines associated with the constraint of type PAR_GEN_Q.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_PAR_GEN_Q_HEADER__
#define __CONSTR_PAR_GEN_Q_HEADER__

#include <math.h>
#include "constr.h"

#define CONSTR_PAR_GEN_Q_PARAM 1e-4

// Function prototypes
Constr* CONSTR_PAR_GEN_Q_new(Net* net);
void CONSTR_PAR_GEN_Q_init(Constr* c);
void CONSTR_PAR_GEN_Q_count_step(Constr* c, Branch* br, int t);
void CONSTR_PAR_GEN_Q_allocate(Constr* c);
void CONSTR_PAR_GEN_Q_clear(Constr* c);
void CONSTR_PAR_GEN_Q_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_PAR_GEN_Q_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_PAR_GEN_Q_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_PAR_GEN_Q_free(Constr* c);

#endif
