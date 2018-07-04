/** @file constr_CFUNC.h
 *  @brief This file lists the constants and routines associated with the constraint of type CFUNC.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_CFUNC_HEADER__
#define __CONSTR_CFUNC_HEADER__

#include <math.h>
#include "func.h"
#include "constr.h"

// Constants
#define CONSTR_CFUNC_BUFFER_SIZE 10
#define CONSTR_CFUNC_EXTRA_VAR_INF 1e8

// Data
typedef struct Constr_CFUNC_Data Constr_CFUNC_Data;

// Function prototypes
Constr* CONSTR_CFUNC_new(Net* net);
void CONSTR_CFUNC_init(Constr* c);
void CONSTR_CFUNC_count_step(Constr* c, Bus* bus, int t);
void CONSTR_CFUNC_analyze_step(Constr* c, Bus* bus, int t);
void CONSTR_CFUNC_clear(Constr* c);
void CONSTR_CFUNC_allocate(Constr* c);
void CONSTR_CFUNC_eval_step(Constr* c, Bus* bus, int t, Vec* v, Vec* ve);
void CONSTR_CFUNC_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_CFUNC_set_parameter(Constr* c, char* key, void* value);
void CONSTR_CFUNC_free(Constr* c);

#endif
