/** @file heur_REG_PF_SWITCH.h
 *  @brief This file lists the constants and routines associated with the heuristic of type REG_PF_SWITCH.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __HEUR_REG_PF_SWITCH_HEADER__
#define __HEUR_REG_PF_SWITCH_HEADER__

#include <math.h>
#include "net.h"
#include "heur.h"

#define HEUR_REG_PF_SWITCH_DEBUG FALSE

// Function prototypes
Heur* HEUR_REG_PF_SWITCH_new(Net* net);
void HEUR_REG_PF_SWITCH_apply_step(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values);

#endif
