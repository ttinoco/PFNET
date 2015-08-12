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
void CONSTR_PAR_GEN_P_init(Constr* c);
void CONSTR_PAR_GEN_P_count_branch(Constr* c, Branch* b);
void CONSTR_PAR_GEN_P_allocate(Constr* c);
void CONSTR_PAR_GEN_P_clear(Constr* c);
void CONSTR_PAR_GEN_P_analyze_branch(Constr* c, Branch* b);
void CONSTR_PAR_GEN_P_eval_branch(Constr* c, Branch *b, Vec* var_values);
void CONSTR_PAR_GEN_P_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_PAR_GEN_P_free(Constr* c);

#endif
