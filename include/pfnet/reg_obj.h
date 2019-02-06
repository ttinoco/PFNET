/** @file reg_obj.h
 *  @brief This file lists the routines for handling regulating objects.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __REG_OBJ_HEADER__
#define __REG_OBJ_HEADER__

#include "stdio.h"
#include "types.h"
#include "gen.h"
#include "bus.h"
#include "conv_vsc.h"
#include "facts.h"

void REG_OBJ_next(char* obj_type, void** obj, Bus* bus);
void REG_OBJ_init(char* obj_type, void** obj, Bus* bus);
void REG_OBJ_set_Q(char obj_type, void* obj, REAL Q, int t);
void REG_OBJ_set_Q_par(char obj_type, void* obj, REAL Q_par);
void REG_OBJ_show(char obj_type, void* obj);
Bus* REG_OBJ_get_bus(char obj_type, void* obj);
int REG_OBJ_get_index_Q(char obj_type, void* obj, int t);
REAL REG_OBJ_get_Q(char obj_type, void* obj, int t);
REAL REG_OBJ_get_Q_max(char obj_type, void* obj);
REAL REG_OBJ_get_Q_min(char obj_type, void* obj);
REAL REG_OBJ_get_Q_par(char obj_type, void* obj);
BOOL REG_OBJ_is_candidate(char obj_type, void* obj);
int REG_OBJ_count_candidates(Bus* bus);

#endif
