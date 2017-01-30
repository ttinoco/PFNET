/** @file heur.h
 *  @brief This file lists the constants and routines associated with the Heur data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __HEUR_HEADER__
#define __HEUR_HEADER__

#include "net.h"
#include "types.h"
#include "list.h"
#include "vector.h"
#include "matrix.h"
#include "constr.h"

// Types
#define HEUR_TYPE_PVPQ 0      // PV-PQ switching

// Heuristic
typedef struct Heur Heur;

// Prototypes
void HEUR_clear_bus_counted(Heur* h, int num);
void HEUR_del(Heur* h);
int HEUR_get_type(Heur* h);
char* HEUR_get_bus_counted(Heur *h);
void* HEUR_get_data(Heur* h);
Heur* HEUR_get_next(Heur* h);
Heur* HEUR_list_add(Heur* hlist, Heur* nh);
void HEUR_list_apply_step(Heur* hlist, Constr* clist, Net* net, Branch* br, int t, Vec* var_values);
void HEUR_list_clear(Heur* hlist, Net* net);
void HEUR_list_del(Heur* hlist);
int HEUR_list_len(Heur* hlist);
Heur* HEUR_new(int type, Net* net);
void HEUR_set_bus_counted(Heur* h, char* counted);
void HEUR_set_data(Heur* h, void* data);
void HEUR_clear(Heur* h, Net* net);
void HEUR_apply_step(Heur* h, Constr* clist, Net* net, Branch* br, int t, Vec* var_values);

#endif
