/** @file contingency.h
 *  @brief This file lists the constants and routines associated with the Cont data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONT_HEADER__
#define __CONT_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "net.h"

// Buffer
#define CONT_BUFFER_SIZE 100      /**< @brief Default contingency buffer size for strings */
#define CONT_JSON_BUFFER_SIZE 200 /**< @brief Default contingency buffer size for json strings */

// Contingency
typedef struct Cont Cont;

void CONT_add_branch_outage(Cont* cont, int br_index);
void CONT_add_gen_outage(Cont* cont, int gen_index);
void CONT_add_load_outage(Cont* cont, int load_index);
void CONT_apply(Cont* cont, Net* net);
void CONT_clear(Cont* cont, Net* net);
void CONT_del(Cont* cont);
char* CONT_get_show_str(Cont* cont);
char* CONT_get_json_string(Cont* cont);
char* CONT_get_name(Cont* cont);
int CONT_get_num_gen_outages(Cont* cont);
int CONT_get_num_branch_outages(Cont* cont);
int CONT_get_num_load_outages(Cont* cont);
int* CONT_get_branch_outages(Cont* cont);
int* CONT_get_gen_outages(Cont* cont);
int* CONT_get_load_outages(Cont* cont);
BOOL CONT_has_gen_outage(Cont* cont, int gen_index);
BOOL CONT_has_branch_outage(Cont* cont, int br_index);
BOOL CONT_has_load_outage(Cont* cont, int load_index);
void CONT_init(Cont* cont);
Cont* CONT_new(void);
void CONT_show(Cont* cont);
void CONT_set_name(Cont* cont, char* name);

#endif
