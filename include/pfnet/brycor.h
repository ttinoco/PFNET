/** @file brycor.h
 *  @brief This file lists the constants and routines associated with the BrYCor data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BRYCOR_HEADER__
#define __BRYCOR_HEADER__

#include <stdio.h>
#include "types.h"

#define BRYCOR_MAX_NUM_VALUES 20
#define BRYCOR_BUFFER_SIZE 100
#define BRYCOR_TYPE_RATIO 0
#define BRYCOR_TYPE_PHASE 1

// Structure
typedef struct BrYCor BrYCor;

// Function prototypes
void BRYCOR_array_del(BrYCor* b_array, int size);
char* BRYCOR_get_name(BrYCor* b);
//BrYCor* BRYCOR_get_copy(BrYCor* b);
int BRYCOR_get_type(BrYCor* b);
int BRYCOR_get_num_values(BrYCor* b);
int BRYCOR_get_max_num_values(BrYCor* b);
REAL BRYCOR_get_value(BrYCor* b, int i);
REAL* BRYCOR_get_values(BrYCor* b);
REAL BRYCOR_get_correction(BrYCor* b, int i);
REAL* BRYCOR_get_corrections(BrYCor* b);
//char* BRYCOR_get_json_string(BrYCor* b, char* output);
BOOL BRYCOR_is_based_on_tap_ratio(BrYCor* b);
BOOL BRYCOR_is_based_on_phase_shift(BrYCor* b);
BrYCor* BRYCOR_new(void);
void BRYCOR_set_type(BrYCor* b, int type);
void BRYCOR_set_name(BrYCor* b, char* name);
void BRYCOR_set_num_values(BrYCor* b, int num);
void BRYCOR_set_value(BrYCor* b, REAL val, int i);
void BRYCOR_set_correction(BrYCor* b, REAL corr, int i);

#endif
