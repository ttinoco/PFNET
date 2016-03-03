/** @file constr_DC_FLOW_LIM.h
 *  @brief This file lists the constants and routines associated with the constraint of type DC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_DC_FLOW_LIM_HEADER__
#define __CONSTR_DC_FLOW_LIM_HEADER__

#include <math.h>
#include "constr.h"

// Function prototypes
void CONSTR_DC_FLOW_LIM_init(Constr* c);
void CONSTR_DC_FLOW_LIM_count_branch(Constr* c, Branch* b);
void CONSTR_DC_FLOW_LIM_allocate(Constr* c);
void CONSTR_DC_FLOW_LIM_clear(Constr* c);
void CONSTR_DC_FLOW_LIM_analyze_branch(Constr* c, Branch* b);
void CONSTR_DC_FLOW_LIM_eval_branch(Constr* c, Branch* b, Vec* var_values);
void CONSTR_DC_FLOW_LIM_store_sens_branch(Constr* c, Branch* b, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_DC_FLOW_LIM_free(Constr* c);

#endif
