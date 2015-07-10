/** @file constr_REG_GEN.h
 *  @brief This file lists the constants and routines associated with the constraint of type REG_GEN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_REG_GEN_HEADER__
#define __CONSTR_REG_GEN_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_REG_GEN_PARAM 1e-8

// Function prototypes
void CONSTR_REG_GEN_init(Constr* c);
void CONSTR_REG_GEN_count_branch(Constr* c, Branch* b);
void CONSTR_REG_GEN_allocate(Constr* c);
void CONSTR_REG_GEN_clear(Constr* c);
void CONSTR_REG_GEN_analyze_branch(Constr* c, Branch* b);
void CONSTR_REG_GEN_eval_branch(Constr* c, Branch *b, Vec* var_values);
void CONSTR_REG_GEN_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_REG_GEN_free(Constr* c);

#endif
