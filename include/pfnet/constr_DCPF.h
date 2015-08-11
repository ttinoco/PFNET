/** @file constr_DCPF.h
 *  @brief This file lists the constants and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_DCPF_HEADER__
#define __CONSTR_DCPF_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
void CONSTR_DCPF_init(Constr* c);
void CONSTR_DCPF_count_branch(Constr* c, Branch* b);
void CONSTR_DCPF_allocate(Constr* c);
void CONSTR_DCPF_clear(Constr* c);
void CONSTR_DCPF_analyze_branch(Constr* c, Branch* b);
void CONSTR_DCPF_eval_branch(Constr* c, Branch *b, Vec* var_values);
void CONSTR_DCPF_store_sens_branch(Constr* c, Branch *b, Vec* sens);
void CONSTR_DCPF_free(Constr* c);

#endif
