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

// Buffer
#define HEUR_BUFFER_SIZE 1024 /**< @brief Default heuristic buffer size for strings */

// Heuristic
typedef struct Heur Heur;

// Prototypes
void HEUR_clear_bus_counted(Heur* h);
void HEUR_del(Heur* h);
char* HEUR_get_name(Heur* h);
char* HEUR_get_bus_counted(Heur *h);
void* HEUR_get_data(Heur* h);
Heur* HEUR_get_next(Heur* h);
Net* HEUR_get_network(Heur* h);

void HEUR_list_clear_error(Heur* hlist);
BOOL HEUR_list_has_error(Heur* hlist);
char* HEUR_list_get_error_string(Heur* hlist);

Heur* HEUR_list_add(Heur* hlist, Heur* nh);
void HEUR_list_apply_step(Heur* hlist, Constr** cptrs, int cnum, Branch* br, int t, Vec* var_values);
void HEUR_list_clear(Heur* hlist);
void HEUR_list_del(Heur* hlist);
int HEUR_list_len(Heur* hlist);
Heur* HEUR_new(Net* net);
void HEUR_init(Heur* h);
void HEUR_clear(Heur* h);
void HEUR_apply(Heur* h, Constr** cptrs, int cnum, Vec* var_values);
void HEUR_apply_step(Heur* h, Constr** cptrs, int cnum, Branch* br, int t, Vec* var_values);

void HEUR_clear_error(Heur * h);
BOOL HEUR_has_error(Heur* h);
char* HEUR_get_error_string(Heur* h);
void HEUR_set_error(Heur* h, char* error_string);

void HEUR_set_bus_counted(Heur* h, char* counted, int size);
void HEUR_set_data(Heur* h, void* data);
void HEUR_set_name(Heur* h, char* name);
void HEUR_set_func_init(Heur* h, void (*func)(Heur* h));
void HEUR_set_func_clear(Heur* h, void (*func)(Heur* h));
void HEUR_set_func_apply_step(Heur* h, void (*func)(Heur* h, Constr** cptrs, int cnum, Branch* br, int t, Vec* var_values));
void HEUR_set_func_free(Heur* h, void (*func)(Heur* h));

void HEUR_update_network(Heur* h);

#endif
