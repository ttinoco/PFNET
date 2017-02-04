/** @file matrix.c
 *  @brief This file defines the Mat data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/matrix.h>
#include <pfnet/vector.h>

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

void MAT_array_del(Mat* m, int size) {
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

Mat* MAT_copy(Mat* m) {

  Mat* newm;
  int k;

  if (!m)
    return NULL;

  newm = MAT_new(m->size1,m->size2,m->nnz);
  for (k = 0; k < m->nnz; k++) {
    newm->row[k] = m->row[k];
    newm->col[k] = m->col[k];
    newm->data[k] = m->data[k];
  }

  return newm; 
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
  MAT_init(m);
  m->size1 = size1;
  m->size2 = size2;
  ARRAY_zalloc(m->row,int,nnz);
  ARRAY_zalloc(m->col,int,nnz);
  ARRAY_zalloc(m->data,REAL,nnz);
  m->nnz = nnz;
  return m;
}

Vec* MAT_rmul_by_vec(Mat* m, Vec* v) {
  
  int k;
  int i;
  int j;
  Vec* w;
  REAL* wdata;
  REAL* vdata;

  if (!m || !v)
    return NULL;

  if (VEC_get_size(v) != m->size2) {
    printf("MAT_rmul_by_vec error: incompatible dimensions\n");
    return NULL;
  }

  w = VEC_new(m->size1);   // all zeros
  wdata = VEC_get_data(w);
  vdata = VEC_get_data(v);
  for (k = 0; k < m->nnz; k++) {
    i = m->row[k];
    j = m->col[k];
    if (i < m->size1 && j < m->size2)
      wdata[i] += m->data[k]*vdata[j];
    else
      printf("MAT_rmul_by_vec warning: invalid index\n");
  }
  return w;
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
  if (m)
    ARRAY_clear(m->data,REAL,m->nnz);
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
  if (m) {
    printf("\nMatrix\n");
    printf("size1 : %d\n",m->size1);
    printf("size2 : %d\n",m->size2);
    printf("nnz   : %d\n",m->nnz);
  }
}
