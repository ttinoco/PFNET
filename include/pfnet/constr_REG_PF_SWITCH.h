/** @file constr_REG_PF_SWITCH.h
 *  @brief This file lists the constants and routines associated with the constraint of type REG_PF_SWITCH.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_REG_PF_SWITCH_HEADER__
#define __CONSTR_REG_PF_SWITCH_HEADER__

#include <math.h>
#include "constr.h"

// Data
typedef struct Constr_REG_PF_SWITCH_Data Constr_REG_PF_SWITCH_Data;

// Function prototypes
Constr* CONSTR_REG_PF_SWITCH_new(Net* net);
void CONSTR_REG_PF_SWITCH_count_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_REG_PF_SWITCH_allocate(Constr* c);
void CONSTR_REG_PF_SWITCH_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t);
void CONSTR_REG_PF_SWITCH_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* v, Vec* ve);
void CONSTR_REG_PF_SWITCH_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_REG_PF_SWITCH_free(Constr* c);

char* CONSTR_REG_PF_SWITCH_get_flags(Constr* c);

#endif
