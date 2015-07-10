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
#include "net.h"
#include "heur.h"

// Data
typedef struct Heur_PVPQ_Data Heur_PVPQ_Data;

// Function prototypes
void HEUR_PVPQ_init(Heur* h, Net* net);
void HEUR_PVPQ_clear(Heur* h, Net* net);
void HEUR_PVPQ_apply_to_branch(Heur* h, Constr* clist, Net* net, Branch* br, Vec* var_values);
void HEUR_PVPQ_free(Heur* h);

#endif
