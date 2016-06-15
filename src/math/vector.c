/** @file vector.c
 *  @brief This file defines the Vec data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/vector.h>

struct Vec {
  int size;
  REAL* data;
};

void VEC_add_to_entry(Vec* v, int index, REAL value) {
  if (v) 
    v->data[index] += value;
}

void VEC_del(Vec* v) {
  if (v)
    free(v->data);
}

REAL VEC_get(Vec* v, int index) {
  if (v)
    return v->data[index];
  else
    return 0;
}

REAL* VEC_get_data(Vec* v) {
  if (v)
    return v->data;
  else
    return NULL;
}

REAL VEC_get_max(Vec* v) {
  int i;
  REAL val = 0;
  if (!v)
    return val;
  else {
    for (i = 0; i < v->size; i++) {
      if (i == 0 || v->data[i] > val)
	val = v->data[i];
    }
    return val;
  }  
}

REAL VEC_get_min(Vec* v) {
  int i;
  REAL val = 0;
  if (!v)
    return val;
  else {
    for (i = 0; i < v->size; i++) {
      if (i == 0 || v->data[i] < val)
	val = v->data[i];
    }
    return val;
  }  
}

int VEC_get_size(Vec* v) {
  if (v)
    return v->size;
  else
    return 0;
}

Vec* VEC_new(int size) {
  Vec* v = (Vec*)malloc(sizeof(Vec));
  v->size = size;
  v->data = (REAL*)calloc(size,sizeof(REAL));
  return v;
}

Vec* VEC_new_from_array(REAL* data, int size) {
  Vec* v = (Vec*)malloc(sizeof(Vec));
  v->size = size;
  v->data = data;
  return v;
}

void VEC_set(Vec* v, int index, REAL value) {
  if (v)
    v->data[index] = value;
}

void VEC_set_zero(Vec* v) {
  int i;
  if (v) {
    for (i = 0; i < v->size; i++) 
      v->data[i] = 0;
  }
}

void VEC_show(Vec* v) {
  if (v) {
    printf("\nVector\n");
    printf("size : %d\n",(int)(v->size));
    printf("max  : %.5e\n",VEC_get_max(v));
    printf("min  : %.5e\n",VEC_get_min(v));
  }
}

void VEC_sub_inplace(Vec* v,Vec* w) {
  
  int k;

  if (!v || !w)
    return;
  
  if (v->size != w->size) {
    printf("VEC_sub_inplace error: incompatible dimensions\n");
    return;
  }
  
  for (k = 0; k < v->size; k++)
    v->data[k] -= w->data[k];  
}
