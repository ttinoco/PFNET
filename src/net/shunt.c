/** @file shunt.c
 *  @brief This file defines the Shunt data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/shunt.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>

struct Shunt {
  
  // Bus
  Bus* bus;       /**< @brief Bus where the shunt is connected */
  Bus* reg_bus;   /**< @brief Bus regulated by this shunt */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Conductance
  REAL g;         /**< @brief Conductance (p.u) */

  // Susceptance
  REAL* b;         /**< @brief Susceptance (p.u.) */
  REAL b_max;     /**< @brief Maximum susceptance (p.u.) */
  REAL b_min;     /**< @brief Minimum susceptance (p.u.) */
  REAL* b_values; /**< @brief Array of valid susceptances (p.u.) */
  char num_b;     /**< @brief Number of valid susceptances (p.u.) */
 
  // Flags
  char vars;      /**< @brief Flags for indicating which quantities are treated as variables **/
  char fixed;     /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;   /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;    /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Indices
  int index;      /**< @brief Shunt index */
  int* index_b;    /**< @brief Susceptance index */
  int* index_y;    /**< @brief Susceptance positive deviation index */
  int* index_z;    /**< @brief Susceptance negative deviation index */

  // List
  Shunt* next;     /**< @brief List of shunts connceted to a bus */
  Shunt* reg_next; /**< @brief List of shunts regulated the same bus */
};

void* SHUNT_array_get(void* shunt_array, int index) { 
  if (shunt_array)
    return (void*)&(((Shunt*)shunt_array)[index]);
  else
    return NULL;
}

void SHUNT_array_del(Shunt* shunt_array, int size) {
  int i;
  Shunt* shunt;
  if (shunt_array) {
    for (i = 0; i < size; i++) {
      shunt = &(shunt_array[i]);
      free(shunt->b_values);
      free(shunt->b);
      free(shunt->index_b);
      free(shunt->index_y);
      free(shunt->index_z);
    }
    free(shunt_array);
  }
}

Shunt* SHUNT_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    Shunt* shunt_array = (Shunt*)malloc(sizeof(Shunt)*size);
    for (i = 0; i < size; i++) {
      SHUNT_init(&(shunt_array[i]),num_periods);
      SHUNT_set_index(&(shunt_array[i]),i);
    }
    return shunt_array;
  }
  else
    return NULL;
}

void SHUNT_array_show(Shunt* shunt_array, int size, int t) { 
  int i;
  if (shunt_array) {
    for (i = 0; i < size; i++) 
      SHUNT_show(&(shunt_array[i]),t);
  }
}

void SHUNT_clear_flags(Shunt* shunt, char flag_type) {
  if (shunt) {
    if (flag_type == FLAG_VARS)
      shunt->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      shunt->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      shunt->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      shunt->sparse = 0x00;
  }
}

int SHUNT_get_num_periods(Shunt* shunt) {
  if (shunt)
    return shunt->num_periods;
  else
    return 0;
}

char SHUNT_get_obj_type(void* shunt) {
  if (shunt)
    return OBJ_SHUNT;
  else
    return OBJ_UNKNOWN;
}

int SHUNT_get_index(Shunt* shunt) {
  if (shunt)
    return shunt->index;
  else
    return 0;
}

int SHUNT_get_index_b(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->index_b[t];
  else
    return 0;
}

int SHUNT_get_index_y(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->index_y[t];
  else
    return 0;
}

int SHUNT_get_index_z(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->index_z[t];
  else
    return 0;
}

Bus* SHUNT_get_bus(Shunt* shunt) {
  if (shunt)
    return shunt->bus;
  else
    return NULL;
}
Bus* SHUNT_get_reg_bus(Shunt* shunt) {
  if (shunt)
    return shunt->reg_bus;
  else
    return NULL;
}

REAL SHUNT_get_g(Shunt* shunt) {
  if (shunt)
    return shunt->g;
  else
    return 0;
}

REAL SHUNT_get_b(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->b[t];
  else
    return 0;
}

REAL SHUNT_get_b_max(Shunt* shunt) {
  if (shunt)
    return shunt->b_max;
  else
    return 0;
}

REAL SHUNT_get_b_min(Shunt* shunt) {
  if (shunt)
    return shunt->b_min;
  else
    return 0;
}

Shunt* SHUNT_get_next(Shunt* shunt) {
  if (shunt)
    return shunt->next;
  else
    return NULL;
}

Shunt* SHUNT_get_reg_next(Shunt* shunt) {
  if (shunt)
    return shunt->reg_next;
  else
    return NULL;
}

void SHUNT_get_var_values(Shunt* shunt, Vec* values, int code) {

  // Local vars
  int t;

  // No shunt
  if (!shunt)
    return;

  // Time loop
  for (t = 0; t < shunt->num_periods; t++) {
    
    if (shunt->vars & SHUNT_VAR_SUSC) { // susceptance
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,shunt->index_b[t],shunt->b_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,shunt->index_b[t],shunt->b_min);
	break;
      default:
	VEC_set(values,shunt->index_b[t],shunt->b[t]);
      }
    }
    if (shunt->vars & SHUNT_VAR_SUSC_DEV) { // susceptance deviations
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,shunt->index_y[t],SHUNT_INF_SUSC);
	VEC_set(values,shunt->index_z[t],SHUNT_INF_SUSC);
	break;
      case LOWER_LIMITS:
	VEC_set(values,shunt->index_y[t],0.);
	VEC_set(values,shunt->index_z[t],0.);
	break;
      default:
	VEC_set(values,shunt->index_y[t],0.);
	VEC_set(values,shunt->index_z[t],0.);
      }
    }    
  }
}

Vec* SHUNT_get_var_indices(void* vshunt, char var) {
  Shunt* shunt = (Shunt*)vshunt;
  Vec* indices;
  int t;
  if (!shunt)
    return NULL;
  if (var == SHUNT_VAR_SUSC) {
    indices = VEC_new(shunt->num_periods);
    for (t = 0; t < shunt->num_periods; t++)
      VEC_set(indices,t,shunt->index_b[t]);
    return indices;
  }
  if (var == SHUNT_VAR_SUSC_DEV) {
    indices = VEC_new(2*shunt->num_periods);
    for (t = 0; t < shunt->num_periods; t++) {
      VEC_set(indices,2*t,shunt->index_y[t]);
      VEC_set(indices,2*t+1,shunt->index_z[t]);
    }
    return indices;
  }
  return NULL;
}

BOOL SHUNT_has_flags(void* vshunt, char flag_type, char mask) {
  Shunt* shunt = (Shunt*)vshunt;
  if (shunt) {
    if (flag_type == FLAG_VARS)
      return (shunt->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (shunt->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (shunt->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (shunt->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL SHUNT_has_properties(void* vshunt, char prop) {
  Shunt* shunt = (Shunt*)vshunt;
  if (!shunt)
    return FALSE;
  if ((prop & SHUNT_PROP_SWITCHED_V) && !SHUNT_is_switched_v(shunt))
    return FALSE;
  return TRUE;
}

void SHUNT_init(Shunt* shunt, int num_periods) {

  // Local vars
  int T;
  
  // No gen
  if (!shunt)
    return;

  T = num_periods;
  shunt->num_periods = num_periods;
 
  shunt->bus = NULL;
  shunt->reg_bus = NULL;
  shunt->g = 0;
  shunt->b_max = 0;
  shunt->b_min = 0;
  shunt->b_values = NULL;
  shunt->num_b = 0;
  shunt->vars = 0x00;
  shunt->fixed = 0x00;
  shunt->bounded = 0x00;
  shunt->sparse = 0x00;
  shunt->index = 0;

  ARRAY_zalloc(shunt->b,REAL,T);
  ARRAY_zalloc(shunt->index_b,int,T);
  ARRAY_zalloc(shunt->index_y,int,T);
  ARRAY_zalloc(shunt->index_z,int,T);

  shunt->next = NULL;
  shunt->reg_next = NULL;
}

BOOL SHUNT_is_fixed(Shunt* shunt) {
  if (!shunt)
    return FALSE;
  else
    return !SHUNT_is_switched(shunt);
}

BOOL SHUNT_is_switched(Shunt* shunt) {
  return SHUNT_is_switched_v(shunt);
}

BOOL SHUNT_is_switched_v(Shunt* shunt) {
  if (!shunt)
    return FALSE;
  else
    return shunt->reg_bus != NULL;
}

Shunt* SHUNT_list_add(Shunt* shunt_list, Shunt* shunt) {
  LIST_add(Shunt,shunt_list,shunt,next);
  return shunt_list;
}

int SHUNT_list_len(Shunt* shunt_list) {
  int len;
  LIST_len(Shunt,shunt_list,next,len);
  return len;
}

Shunt* SHUNT_list_reg_add(Shunt* reg_shunt_list, Shunt* reg_shunt) {
  LIST_add(Shunt,reg_shunt_list,reg_shunt,reg_next);
  return reg_shunt_list;
}

int SHUNT_list_reg_len(Shunt* reg_shunt_list) {
  int len;
  LIST_len(Shunt,reg_shunt_list,reg_next,len);
  return len;
}

Shunt* SHUNT_new(int num_periods) { 
  if (num_periods > 0) {
    Shunt* shunt = (Shunt*)malloc(sizeof(Shunt));
    SHUNT_init(shunt,num_periods);
    return shunt;
  }
  else
    return NULL;
}

void SHUNT_set_bus(Shunt* shunt, Bus* bus) { 
  if (shunt)
    shunt->bus = (Bus*)bus;
}

void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus) { 
  if (shunt)
    shunt->reg_bus = (Bus*)reg_bus;
}

void SHUNT_set_index(Shunt* shunt, int index) { 
  if (shunt)
    shunt->index = index;
}

void SHUNT_set_g(Shunt* shunt, REAL g) { 
  if (shunt)
    shunt->g = g;
}

void SHUNT_set_b(Shunt* shunt, REAL b, int t) { 
  if (shunt && t >= 0 && t < shunt->num_periods)
    shunt->b[t] = b;
}

void SHUNT_set_b_max(Shunt* shunt, REAL b_max) {
  if (shunt)
    shunt->b_max = b_max;
}

void SHUNT_set_b_min(Shunt* shunt, REAL b_min) {
  if (shunt)
    shunt->b_min = b_min;
}

void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num, REAL norm) {
  int i;
  if (shunt) {
    shunt->b_values = (REAL*)malloc(sizeof(REAL)*num);
    shunt->num_b = num;
    for (i = 0; i < shunt->num_b; i++) 
      shunt->b_values[i] = values[i]/norm; // note normalization
  }
}

int SHUNT_set_flags(void* vshunt, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Shunt* shunt = (Shunt*)vshunt;
  int t;

  // No shunt
  if (!shunt)
    return index;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(shunt->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(shunt->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(shunt->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(shunt->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & SHUNT_VAR_SUSC) && (mask & SHUNT_VAR_SUSC)) { // shunt susceptance
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < shunt->num_periods; t++)
	shunt->index_b[t] = index+t;
    }
    (*flags_ptr) |= SHUNT_VAR_SUSC;
    index += shunt->num_periods;
  }
  if (!((*flags_ptr) & SHUNT_VAR_SUSC_DEV) && (mask & SHUNT_VAR_SUSC_DEV)) { // shunt susceptance deviations
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < shunt->num_periods; t++) {
	shunt->index_y[t] = index+2*t;
	shunt->index_z[t] = index+2*t+1;
      }
    }
    (*flags_ptr) |= SHUNT_VAR_SUSC_DEV;
    index += 2*shunt->num_periods;
  }
  return index;
}

void SHUNT_set_var_values(Shunt* shunt, Vec* values) {

  // Local vars
  int t;

  // No shunt
  if (!shunt)
    return;

  // Time loop
  for (t = 0; t < shunt->num_periods; t++) {

    if (shunt->vars & SHUNT_VAR_SUSC) // shunt susceptance (p.u.)
      shunt->b[t] = VEC_get(values,shunt->index_b[t]); 
  }
}

void SHUNT_show(Shunt* shunt, int t) { 
  if (shunt)
    printf("shunt %d\t%d\n",
	   BUS_get_number(shunt->bus),
	   shunt->index);
}

void SHUNT_propagate_data_in_time(Shunt* shunt) {
  int t;
  if (shunt) {
    for (t = 1; t < shunt->num_periods; t++)
      shunt->b[t] = shunt->b[0];
  }
}
