/** @file heur_PVPQ.h
 *  @brief This file lists the constants and routines associated with the heuristic of type PVPQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __HEUR_PVPQ_HEADER__
#define __HEUR_PVPQ_HEADER__

#include <math.h>
#include "heur.h"

#define HEUR_PVPQ_DEBUG FALSE

// Function prototypes
Heur* HEUR_PVPQ_new(Net* net);
void HEUR_PVPQ_init(Heur* h);
void HEUR_PVPQ_clear(Heur* h);
void HEUR_PVPQ_apply_step(Heur* h, Constr* clist, Branch* br, int t, Vec* var_values);
void HEUR_PVPQ_free(Heur* h);

#endif
