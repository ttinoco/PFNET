/** @file vector.h
 *  @brief This file lists the constants and routines associated with the Vec data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __VEC_HEADER__
#define __VEC_HEADER__

#include <stdio.h>
#include <math.h>
#include "types.h"

// Vector
typedef struct Vec Vec;

// Function prototypes
void VEC_add_to_entry(Vec* v, int index, REAL value);
void VEC_del(Vec* v);
REAL VEC_get(Vec* v, int index);
REAL* VEC_get_data(Vec* v);
REAL VEC_get_max(Vec* v);
REAL VEC_get_min(Vec* v);
int VEC_get_size(Vec* v);
Vec* VEC_new(int size);
Vec* VEC_new_from_array(REAL* data, int size);
void VEC_set(Vec* v, int index, REAL value);
void VEC_set_zero(Vec* v);
void VEC_show(Vec* v);
void VEC_sub_inplace(Vec* v,Vec* w);

#endif
