/** @file constr_REG_TRAN.h
 *  @brief This file lists the constants and routines associated with the constraint of type REG_TRAN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_REG_TRAN_HEADER__
#define __CONSTR_REG_TRAN_HEADER__

#include <math.h>
#include "constr.h"

// Parameters
#define CONSTR_REG_TRAN_PARAM 1e-8
#define CONSTR_REG_TRAN_NORM 1e0

// Function prototypes
Constr* CONSTR_REG_TRAN_new(Net* net);
void CONSTR_REG_TRAN_init(Constr* c);
void CONSTR_REG_TRAN_count_step(Constr* c, Branch* br, int t);
void CONSTR_REG_TRAN_allocate(Constr* c);
void CONSTR_REG_TRAN_clear(Constr* c);
void CONSTR_REG_TRAN_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_REG_TRAN_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_REG_TRAN_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_REG_TRAN_free(Constr* c);

#endif
