/** @file matrix.c
 *  @brief This file defines the Mat data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/matrix.h>

struct Mat {

  int size1;
  int size2;

  int* row;
  int* col;
  REAL* data;
  int nnz;
  
  BOOL owns_rowcol;
};

void MAT_add_to_dentry(Mat* m, int index, REAL value) {
  if (m)
    m->data[index] += value;
}

void MAT_array_del(Mat *m, int size) {
  int i;
  if (m) {
    for (i = 0; i < size; i++) {
      if (m[i].owns_rowcol) {
	free(m[i].row);
	free(m[i].col);
      }
      free(m[i].data);
    }
    free(m);
  }
}

Mat* MAT_array_get(Mat* m, int index) {
  if (m)
    return &(m[index]);
  else
    return NULL;
}

Mat* MAT_array_new(int size) {
  int i;
  Mat* m = (Mat*)malloc(sizeof(Mat)*size);
  for (i = 0; i < size; i++) 
    MAT_init(&(m[i]));
  return m;  
}

void MAT_array_set_zero_d(Mat* m, int size) {
  int i;
  if (m) {
    for (i = 0; i < size; i++) 
      MAT_set_zero_d(&(m[i]));
  }
}

void MAT_del(Mat* m) {
  if (m) {
    if (m->owns_rowcol) {
      free(m->row);
      free(m->col);
    }
    free(m->data);
    free(m);
  }
}

int MAT_get_i(Mat* m, int index) {
  return m->row[index];
}

int MAT_get_j(Mat* m, int index) {
  return m->col[index];
}

REAL MAT_get_d(Mat* m, int index) {
  return m->data[index];
}

int MAT_get_nnz(Mat* m) {
  if (m)
    return m->nnz;
  else
    return 0;
}

int MAT_get_size1(Mat* m) {
  if (m)
    return m->size1;
  else
    return 0;
}

int MAT_get_size2(Mat* m) {
  if (m)
    return m->size2;
  else
    return 0;
}

int* MAT_get_row_array(Mat* m) {
  if (m)
    return m->row;
  else
    return NULL;  
}

int* MAT_get_col_array(Mat* m) {
  if (m)
    return m->col;
  else
    return NULL;  
}

REAL* MAT_get_data_array(Mat* m) {
  if (m)
    return m->data;
  else
    return NULL;  
}

void MAT_init(Mat* m) {
  if (m) {
    m->size1 = 0;
    m->size2 = 0;
    m->row = NULL;
    m->col = NULL;
    m->data = NULL;
    m->nnz = 0;
    m->owns_rowcol = TRUE;
  }
}

Mat* MAT_new(int size1, int size2, int nnz) {
  Mat* m = (Mat*)malloc(sizeof(Mat));
  m->size1 = size1;
  m->size2 = size2;
  m->row = (int*)calloc(nnz,sizeof(int));
  m->col = (int*)calloc(nnz,sizeof(int));
  m->data = (REAL*)calloc(nnz,sizeof(REAL));
  m->nnz = nnz;
  return m;
}

void MAT_set_i(Mat* m, int index, int value) {
  if (m)
    m->row[index] = value;
}

void MAT_set_j(Mat* m, int index, int value) {
  if (m)
    m->col[index] = value;
}

void MAT_set_d(Mat* m, int index, REAL value) {
  if (m)
    m->data[index] = value;
}

void MAT_set_size1(Mat* m, int size1) {
  if (m)
    m->size1 = size1;
}

void MAT_set_size2(Mat* m, int size2) {
  if (m)
    m->size2 = size2;
}

void MAT_set_zero_d(Mat* m) {
  int i;
  if (m) {
    for (i = 0; i < m->nnz; i++) 
      m->data[i] = 0.;
  }
}

void MAT_set_row_array(Mat* m, int* array) {
  if (m)
    m->row = array;
}

void MAT_set_col_array(Mat* m, int* array) {
  if (m)
    m->col = array;
}

void MAT_set_data_array(Mat* m, REAL* array) {
  if (m)
    m->data = array;
}

void MAT_set_nnz(Mat* m, int nnz) {
  if (m)
    m->nnz = nnz;
}

void MAT_set_owns_rowcol(Mat* m, BOOL flag) {
  if (m)
    m->owns_rowcol = flag;
}

void MAT_show(Mat* m) {
  int i;
  if (m) {
    printf("\nMatrix\n");
    printf("size1 : %d\n",m->size1);
    printf("size2 : %d\n",m->size2);
    printf("nnz   : %d\n",m->nnz);
    //for (i = 0; i < m->nnz; i++) 
    //  printf("(%d,%d,%.2e)\n",m->row[i],m->col[i],m->data[i]);
  }
}
