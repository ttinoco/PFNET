/** @file constr_PVPQ_SWITCHING.h
 *  @brief This file lists the constants and routines associated with the constraint of type PVPQ_SWITCHING.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_PVPQ_SWITCHING_HEADER__
#define __CONSTR_PVPQ_SWITCHING_HEADER__

#include <math.h>
#include "constr.h"

#define CONSTR_PVPQ_SWITCHING_PARAM 1e-4

// Data
typedef struct Constr_PVPQ_SWITCHING_Data Constr_PVPQ_SWITCHING_Data;

// Function prototypes
Constr* CONSTR_PVPQ_SWITCHING_new(Net* net);
void CONSTR_PVPQ_SWITCHING_count_step(Constr* c, Branch* br, int t);
void CONSTR_PVPQ_SWITCHING_allocate(Constr* c);
void CONSTR_PVPQ_SWITCHING_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_PVPQ_SWITCHING_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_PVPQ_SWITCHING_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_PVPQ_SWITCHING_free(Constr* c);

char* CONSTR_PVPQ_SWITCHING_get_flags(Constr* c);

#endif
