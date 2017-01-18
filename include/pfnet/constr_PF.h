/** @file constr_PF.h
 *  @brief This file lists the constants and routines associated with the constraint of type PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_PF_HEADER__
#define __CONSTR_PF_HEADER__

#include <math.h>
#include "constr.h"

// Data
typedef struct Constr_PF_Data Constr_PF_Data;

// Function prototypes
void CONSTR_PF_init(Constr* c);
void CONSTR_PF_count_step(Constr* c, Branch* br, int t);
void CONSTR_PF_allocate(Constr* c);
void CONSTR_PF_clear(Constr* c);
void CONSTR_PF_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_PF_eval_step(Constr* c, Branch* br, int t, Vec* var_values);
void CONSTR_PF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_PF_free(Constr* c);

#endif
