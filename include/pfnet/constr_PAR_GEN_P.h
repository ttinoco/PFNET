/** @file constr_PAR_GEN_P.h
 *  @brief This file lists the constants and routines associated with the constraint of type PAR_GEN_P.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_PAR_GEN_P_HEADER__
#define __CONSTR_PAR_GEN_P_HEADER__

#include <math.h>
#include "constr.h"

#define CONSTR_PAR_GEN_P_PARAM 1e-4

// Function prototypes
Constr* CONSTR_PAR_GEN_P_new(Net* net);
void CONSTR_PAR_GEN_P_init(Constr* c);
void CONSTR_PAR_GEN_P_count_step(Constr* c, Branch* br, int t);
void CONSTR_PAR_GEN_P_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_PAR_GEN_P_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_PAR_GEN_P_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_PAR_GEN_P_free(Constr* c);

#endif
