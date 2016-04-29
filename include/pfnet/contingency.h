/** @file contingency.h
 *  @brief This file lists the constants and routines associated with the Cont data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONT_HEADER__
#define __CONT_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"

// Buffer
#define CONT_BUFFER_SIZE 64 /**< @brief Default contingency buffer size for strings */

// Contingency
typedef struct Cont Cont;

// Other
typedef struct Gen Gen;
typedef struct Branch Branch;

void CONT_add_branch_outage(Cont* cont, Branch* br);
void CONT_add_gen_outage(Cont* cont, Gen* gen);
void CONT_apply(Cont* cont);
void CONT_clear(Cont* cont);
void CONT_del(Cont* cont);
int CONT_get_num_gen_outages(Cont* cont);
int CONT_get_num_branch_outages(Cont* cont);
BOOL CONT_has_gen_outage(Cont* cont, Gen* gen);
BOOL CONT_has_branch_outage(Cont* cont, Branch* br);
void CONT_init(Cont* cont);
Cont* CONT_new(void);
void CONT_show(Cont* cont);
char* CONT_get_show_str(Cont* cont);

#endif
