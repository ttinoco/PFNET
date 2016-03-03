/** @file constr_PAR_GEN_Q.h
 *  @brief This file lists the constants and routines associated with the constraint of type PAR_GEN_Q.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_PAR_GEN_Q_HEADER__
#define __CONSTR_PAR_GEN_Q_HEADER__

#include <math.h>
#include "constr.h"

#define CONSTR_PAR_GEN_Q_PARAM 1e-4

// Function prototypes
void CONSTR_PAR_GEN_Q_init(Constr* c);
void CONSTR_PAR_GEN_Q_count_branch(Constr* c, Branch* b);
void CONSTR_PAR_GEN_Q_allocate(Constr* c);
void CONSTR_PAR_GEN_Q_clear(Constr* c);
void CONSTR_PAR_GEN_Q_analyze_branch(Constr* c, Branch* b);
void CONSTR_PAR_GEN_Q_eval_branch(Constr* c, Branch* b, Vec* var_values);
void CONSTR_PAR_GEN_Q_store_sens_branch(Constr* c, Branch* b, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_PAR_GEN_Q_free(Constr* c);

#endif
