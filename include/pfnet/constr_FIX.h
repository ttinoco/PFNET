/** @file constr_FIX.h
 *  @brief This file lists the constants and routines associated with the constraint of type FIX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_FIX_HEADER__
#define __CONSTR_FIX_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
void CONSTR_FIX_init(Constr* c);
void CONSTR_FIX_count_branch(Constr* c, Branch* b);
void CONSTR_FIX_allocate(Constr* c);
void CONSTR_FIX_clear(Constr* c);
void CONSTR_FIX_analyze_branch(Constr* c, Branch* b);
void CONSTR_FIX_eval_branch(Constr* c, Branch* b, Vec* var_values);
void CONSTR_FIX_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_FIX_free(Constr* c);

#endif
