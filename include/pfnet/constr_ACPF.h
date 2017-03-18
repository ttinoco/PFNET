/** @file constr_ACPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type ACPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_ACPF_HEADER__
#define __CONSTR_ACPF_HEADER__

#include <math.h>
#include "constr.h"

// Data
typedef struct Constr_ACPF_Data Constr_ACPF_Data;

// Function prototypes
Constr* CONSTR_ACPF_new(Net* net);
void CONSTR_ACPF_init(Constr* c);
void CONSTR_ACPF_count_step(Constr* c, Branch* br, int t);
void CONSTR_ACPF_allocate(Constr* c);
void CONSTR_ACPF_clear(Constr* c);
void CONSTR_ACPF_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_ACPF_eval_step(Constr* c, Branch* br, int t, Vec* v);
void CONSTR_ACPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_ACPF_free(Constr* c);

#endif
