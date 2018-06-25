/** @file heur_PVPQ_SWITCHING.h
 *  @brief This file lists the constants and routines associated with the heuristic of type PVPQ_SWITCHING.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __HEUR_PVPQ_SWITCHING_HEADER__
#define __HEUR_PVPQ_SWITCHING_HEADER__

#include <math.h>
#include "heur.h"

#define HEUR_PVPQ_SWITCHING_DEBUG FALSE

// Function prototypes
Heur* HEUR_PVPQ_SWITCHING_new(Net* net);
void HEUR_PVPQ_SWITCHING_init(Heur* h);
void HEUR_PVPQ_SWITCHING_clear(Heur* h);
void HEUR_PVPQ_SWITCHING_apply_step(Heur* h, Constr** cptrs, int cnum, Branch* br, int t, Vec* var_values);
void HEUR_PVPQ_SWITCHING_free(Heur* h);

#endif
