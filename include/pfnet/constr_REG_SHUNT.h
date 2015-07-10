/** @file constr_REG_SHUNT.h
 *  @brief This file lists the constants and routines associated with the constraint of type REG_SHUNT.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_REG_SHUNT_HEADER__
#define __CONSTR_REG_SHUNT_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_REG_SHUNT_PARAM 1e-8
#define CONSTR_REG_SHUNT_NORM 1e0

// Function prototypes
void CONSTR_REG_SHUNT_init(Constr* c);
void CONSTR_REG_SHUNT_count_branch(Constr* c, Branch* b);
void CONSTR_REG_SHUNT_allocate(Constr* c);
void CONSTR_REG_SHUNT_clear(Constr* c);
void CONSTR_REG_SHUNT_analyze_branch(Constr* c, Branch* b);
void CONSTR_REG_SHUNT_eval_branch(Constr* c, Branch *b, Vec* var_values);
void CONSTR_REG_SHUNT_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_REG_SHUNT_free(Constr* c);

#endif
