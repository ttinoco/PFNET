/** @file load.h
 *  @brief This file lists the constants and routines associated with the Load data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __LOAD_HEADER__
#define __LOAD_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"

// Load
typedef struct Load Load;

// Other
typedef struct Bus Bus;

void* LOAD_array_get(void* load, int index);
Load* LOAD_array_new(int num);
void LOAD_array_show(Load* load, int num);
char LOAD_get_obj_type(Load* load);
Bus* LOAD_get_bus(Load* load);
int LOAD_get_index(Load* load);
Load* LOAD_get_next(Load* load);
REAL LOAD_get_P(Load* load);
REAL LOAD_get_Q(Load* load);
void LOAD_init(Load* load);
Load* LOAD_list_add(Load *load_list, Load* load);
int LOAD_list_len(Load* load_list);
Load* LOAD_new(void);
void LOAD_set_bus(Load* load, Bus* bus);
void LOAD_set_index(Load* load, int index);
void LOAD_set_P(Load* load, REAL P);
void LOAD_set_Q(Load* load, REAL Q);
void LOAD_show(Load* load);

#endif
