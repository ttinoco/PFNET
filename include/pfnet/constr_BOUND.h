/** @file constr_BOUND.h
 *  @brief This file lists the constants and routines associated with the constraint of type BOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_BOUND_HEADER__
#define __CONSTR_BOUND_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_BOUND_PARAM 1e-4

// Function prototypes
void CONSTR_BOUND_init(Constr* c);
void CONSTR_BOUND_count_branch(Constr* c, Branch* b);
void CONSTR_BOUND_allocate(Constr* c);
void CONSTR_BOUND_clear(Constr* c);
void CONSTR_BOUND_analyze_branch(Constr* c, Branch* b);
void CONSTR_BOUND_eval_branch(Constr* c, Branch *b, Vec* var_values);
void CONSTR_BOUND_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_BOUND_free(Constr* c);

#endif
