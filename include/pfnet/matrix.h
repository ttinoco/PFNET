/** @file matrix.h
 *  @brief This file lists the constants and routines associated with the Mat data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __MAT_HEADER__
#define __MAT_HEADER__

#include <stdio.h>
#include "types.h"

// Types
typedef struct Mat Mat;
typedef struct Vec Vec;

// Function prototypes
void MAT_add_to_dentry(Mat* m, int index, REAL value);
void MAT_array_del(Mat* m, int size);
Mat* MAT_array_new(int size);
Mat* MAT_array_get(Mat* m, int index);
void MAT_array_set_zero_d(Mat* m, int size);
Mat* MAT_copy(Mat* m);
void MAT_del(Mat* m);
int MAT_get_i(Mat* m, int index);
int MAT_get_j(Mat* m, int index);
REAL MAT_get_d(Mat* m, int index);
int MAT_get_nnz(Mat* m);
int MAT_get_size1(Mat* m);
int MAT_get_size2(Mat* m);
int* MAT_get_row_array(Mat* m);
int* MAT_get_col_array(Mat* m);
REAL* MAT_get_data_array(Mat* m);
void MAT_init(Mat* m);
Mat* MAT_new(int size1, int size2, int nnz);
Vec* MAT_rmul_by_vec(Mat* m, Vec* v);
void MAT_set_i(Mat* m, int index, int value);
void MAT_set_j(Mat* m, int index, int value);
void MAT_set_d(Mat* m, int index, REAL value);
void MAT_set_size1(Mat* m, int size1);
void MAT_set_size2(Mat* m, int size2);
void MAT_set_zero_d(Mat* m);
void MAT_set_row_array(Mat* m, int* array);
void MAT_set_col_array(Mat* m, int* array);
void MAT_set_data_array(Mat* m, REAL* array);
void MAT_set_nnz(Mat* m, int nnz);
void MAT_set_owns_rowcol(Mat* m, BOOL flag);
void MAT_show(Mat* m);

#endif
