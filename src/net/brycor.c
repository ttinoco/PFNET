/** @file brycor.c
 *  @brief This file defines the BrYCor data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/brycor.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct BrYCor {

  // Properties
  char name[BRYCOR_BUFFER_SIZE]; /**< @brief Name of branch admittance correction table */
  int type;                      /**< @brief Type of correction table: based on tap ratio or phase shift */
  
  // Values and corrections
  int num_values;                         /**< @brief Number of valid values in table */
  REAL value[BRYCOR_MAX_NUM_VALUES];      /**< @brief Array of values (tap ratios in p.u. or phase shifts in radians) */
  REAL correction[BRYCOR_MAX_NUM_VALUES]; /**< @brief Array of corrections (scaling factor of branch series admittance) */
};

void BRYCOR_array_del(BrYCor* b_array, int size) {
  if (b_array)
    free(b_array);
}

char* BRYCOR_get_name(BrYCor* b) {
  if (b)
    return b->name;
  else
    return NULL;
}

int BRYCOR_get_type(BrYCor* b) {
  if (b)
    return b->type;
  else
    return BRYCOR_TYPE_RATIO;
}

int BRYCOR_get_num_values(BrYCor* b) {
  if (b)
    return b->num_values;
  else
    return 0;
}

int BRYCOR_get_max_num_values(BrYCor* b) {
  if (b)
    return BRYCOR_MAX_NUM_VALUES;
  else
    return 0;
}

REAL BRYCOR_get_value(BrYCor* b, int i) {
  if (b && 0 <= i && i < b->num_values)
    return b->value[i];
  else
    return 0;
}

REAL* BRYCOR_get_values(BrYCor* b) {
  if (b)
    return b->value;
  else
    return NULL;
}

REAL BRYCOR_get_correction(BrYCor* b, int i) {
  if (b && 0 <= i && i < b->num_values)
    return b->correction[i];
  else
    return 0;
}

REAL* BRYCOR_get_corrections(BrYCor* b) {
  if (b)
    return b->correction;
  else
    return NULL;
}

BOOL BRYCOR_is_based_on_tap_ratio(BrYCor* b) {
  return BRYCOR_get_type(b) == BRYCOR_TYPE_RATIO;
}

BOOL BRYCOR_is_based_on_phase_shift(BrYCor* b) {
  return BRYCOR_get_type(b) == BRYCOR_TYPE_PHASE;
}

BrYCor* BRYCOR_new(void) {
  int i;
  BrYCor* b = (BrYCor*)malloc(sizeof(BrYCor));
  strcpy(b->name,"");
  b->type = BRYCOR_TYPE_RATIO;
  b->num_values = 0;
  for (i = 0; i < BRYCOR_MAX_NUM_VALUES; i++) {
    b->value[i] = 0;
    b->correction[i] = 0;
  }
  return b;
}

void BRYCOR_set_type(BrYCor* b, int type) {
  if (b)
    b->type = type;
}

void BRYCOR_set_name(BrYCor* b, char* name) {
  if (b)
    strncpy(b->name,name,(size_t)(BRYCOR_BUFFER_SIZE-1));
}

void BRYCOR_set_num_values(BrYCor* b, int num) {
  if (b)
    b->num_values = (num <= BRYCOR_MAX_NUM_VALUES) ? num : BRYCOR_MAX_NUM_VALUES;
}

void BRYCOR_set_value(BrYCor* b, REAL val, int i) {
  if (b && 0 <= i && i < b->num_values)
    b->value[i] = val;
}

void BRYCOR_set_correction(BrYCor* b, REAL corr, int i) {
  if (b && 0 <= i && i < b->num_values)
    b->correction[i] = corr;
}
