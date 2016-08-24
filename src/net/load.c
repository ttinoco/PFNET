/** @file load.c
 *  @brief This file defines the Load data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/load.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>

struct Load {

  // Bus
  Bus* bus;            /**< @brief Bus to which the load is connected */

  // Times
  int num_periods;   /**< @brief Number of time periods. */
  
  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Active power
  REAL* P;              /**< @brief Load active power consumption (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum load active power consumption (p.u.) */
  REAL P_min;          /**< @brief Minimum load active power consumption (p.u.) */

  // Reactive power
  REAL* Q;              /**< @brief Load reactive power (p.u. system base power) */

  // Utility
  REAL util_coeff_Q0;  /**< @brief Load utility coefficient (constant term, units of $/hr ) */
  REAL util_coeff_Q1;  /**< @brief Load utility coefficient (linear term, units of $/(hr p.u.) ) */
  REAL util_coeff_Q2;  /**< @brief Load utility coefficient (quadratic term, units of $/(hr p.u.^2) ) */

  // Indices
  int index;           /**< @brief Load index */
  int* index_P;         /**< @brief Active power index */

  // Sensitivities
  REAL* sens_P_u_bound;  /**< @brief Sensitivity of active power upper bound */
  REAL* sens_P_l_bound;  /**< @brief Sensitivity of active power lower bound */

  // List
  Load* next;           /**< @brief List of loads connected to a bus */
};

void* LOAD_array_get(void* load_array, int index) { 
  if (load_array) 
    return (void*)&(((Load*)load_array)[index]);
  else
    return NULL;
}

void LOAD_array_del(Load* load_array, int size) {
  int i;
  Load* load;
  if (load_array) {
    for (i = 0; i < size; i++) {
      load = &(load_array[i]);
      free(load->P);
      free(load->Q);
      free(load->index_P);
      free(load->sens_P_u_bound);
      free(load->sens_P_l_bound);
    }
    free(load_array);
  }  
}

Load* LOAD_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    Load* load_array = (Load*)malloc(sizeof(Load)*size);
    for (i = 0; i < size; i++) {
      LOAD_init(&(load_array[i]),num_periods);
      LOAD_set_index(&(load_array[i]),i);
    }
    return load_array;
  }
  else
    return NULL;
}

void LOAD_array_show(Load* load_array, int size, int t) { 
  int i;
  if (load_array) {
    for (i = 0; i < size; i++) 
      LOAD_show(&(load_array[i]),t);
  }
}

void LOAD_clear_sensitivities(Load* load) {
  int t;
  if (load) {
    for (t = 0; t < load->num_periods; t++) {
      load->sens_P_u_bound[t] = 0;
      load->sens_P_l_bound[t] = 0;
    }
  }
}

void LOAD_clear_flags(Load* load, char flag_type) {
  if (load) {
    if (flag_type == FLAG_VARS)
      load->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      load->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      load->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      load->sparse = 0x00;
  }
}

REAL LOAD_get_sens_P_u_bound(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->sens_P_u_bound[t];
  else
    return 0;
}

REAL LOAD_get_sens_P_l_bound(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->sens_P_l_bound[t];
  else
    return 0;
}

char LOAD_get_obj_type(void* load) {
  if (load)
    return OBJ_LOAD;
  else
    return OBJ_UNKNOWN;
}

Bus* LOAD_get_bus(Load* load) {
  if (load)
    return load->bus;
  else
    return NULL;
}

REAL LOAD_get_P_util(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return LOAD_get_P_util_for(load,load->P[t]); // $/hr
  else
    return 0;
}

REAL LOAD_get_P_util_for(Load* load, REAL P) {
  if (load)
    return (load->util_coeff_Q0 + 
	    load->util_coeff_Q1*P +
	    load->util_coeff_Q2*pow(P,2.)); // $/hr
  else
    return 0;
}

REAL LOAD_get_util_coeff_Q0(Load* load) {
  if (load)
    return load->util_coeff_Q0;
  else
    return 0;
}

REAL LOAD_get_util_coeff_Q1(Load* load) {
  if (load)
    return load->util_coeff_Q1;
  else
    return 0;
}

REAL LOAD_get_util_coeff_Q2(Load* load) {
  if (load)
    return load->util_coeff_Q2;
  else
    return 0;
}

int LOAD_get_index(Load* load) {
  if (load)
    return load->index;
  else
    return 0;
}

int LOAD_get_index_P(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->index_P[t];
  else
    return 0;
}

Load* LOAD_get_next(Load* load) {
  if (load)
    return load->next;
  else
    return NULL;
}

REAL LOAD_get_P(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->P[t];
  else
    return 0;
}

REAL LOAD_get_P_max(Load* load) {
  if (load)
    return load->P_max;
  else 
    return 0;
}

REAL LOAD_get_P_min(Load* load) {
  if (load)
    return load->P_min;
  else 
    return 0;
}

REAL LOAD_get_Q(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->Q[t];
  else
    return 0;
}

void LOAD_get_var_values(Load* load, Vec* values, int code) {
 
  // Local vars
  int t;
 
  // No laod
  if (!load)
    return;

  // Time loop
  for (t = 0; t < load->num_periods; t++) {

    if (load->vars & LOAD_VAR_P) { // active power
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,load->index_P[t],load->P_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,load->index_P[t],load->P_min);
	break;
      default:
	VEC_set(values,load->index_P[t],load->P[t]);
      }
    }
  }
}

Vec* LOAD_get_var_indices(void* vload, char var) {
  Load* load = (Load*)vload;
  Vec* indices;
  int t;
  if (!load)
    return NULL;
  if (var == LOAD_VAR_P) {
    indices = VEC_new(load->num_periods);
    for (t = 0; t < load->num_periods; t++)
      VEC_set(indices,t,load->index_P[t]);
    return indices;
  }
  return NULL;
}

BOOL LOAD_has_flags(void* vload, char flag_type, char mask) {
  Load* load = (Load*)vload;
  if (load) {
    if (flag_type == FLAG_VARS)
      return (load->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (load->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (load->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (load->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL LOAD_has_properties(void* vload, char prop) {
  Load* load = (Load*)vload;
  if (!load)
    return FALSE;
  if ((prop & LOAD_PROP_P_ADJUST) && !LOAD_is_P_adjustable(load))
    return FALSE;
  return TRUE;
}

void LOAD_init(Load* load, int num_periods) {

  // Local vars
  int T;

  // No load
  if (!load)
    return;

  T = num_periods;
  load->num_periods = num_periods;
    
  load->bus = NULL;
  
  load->fixed = 0x00;
  load->bounded = 0x00;
  load->sparse = 0x00;
  load->vars = 0x00;
  
  load->P_max = 0;
  load->P_min = 0;
  
  load->util_coeff_Q0 = 0;
  load->util_coeff_Q1 = 20000.;
  load->util_coeff_Q2 = -100.;
  
  load->index = 0;

  ARRAY_zalloc(load->P,REAL,T);
  ARRAY_zalloc(load->Q,REAL,T);
  ARRAY_zalloc(load->index_P,int,T);
  ARRAY_zalloc(load->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(load->sens_P_l_bound,REAL,T);
  
  load->next = NULL;
}

BOOL LOAD_is_P_adjustable(Load* load) {
  if (load)
    return load->P_min < load->P_max;
  else
    return FALSE;
}

Load* LOAD_list_add(Load* load_list, Load* load) {
  LIST_add(Load,load_list,load,next);
  return load_list;
}

int LOAD_list_len(Load* load_list) {
  int len;
  LIST_len(Load,load_list,next,len);
  return len;
}

Load* LOAD_new(int num_periods) {
  if (num_periods > 0) {
    Load* load = (Load*)malloc(sizeof(Load));
    LOAD_init(load,num_periods);
    return load;
  }
  else
    return NULL;
}

void LOAD_set_sens_P_u_bound(Load* load, REAL value, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->sens_P_u_bound[t] = value;
}

void LOAD_set_sens_P_l_bound(Load* load, REAL value, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->sens_P_l_bound[t] = value;
}

void LOAD_set_util_coeff_Q0(Load* load, REAL q) {
  if (load)
    load->util_coeff_Q0 = q;
}

void LOAD_set_util_coeff_Q1(Load* load, REAL q) {
  if (load)
    load->util_coeff_Q1 = q;
}

void LOAD_set_util_coeff_Q2(Load* load, REAL q) {
  if (load)
    load->util_coeff_Q2 = q;
}

void LOAD_set_bus(Load* load, Bus* bus) { 
  if (load)
    load->bus = (Bus*)bus;
}

void LOAD_set_index(Load* load, int index) { 
  if (load)
    load->index = index;
}

void LOAD_set_P(Load* load, REAL P, int t) { 
  if (load && t >= 0 && t < load->num_periods)
    load->P[t] = P;
}

void LOAD_set_P_max(Load* load, REAL P_max) {
  if (load)
    load->P_max = P_max;
}

void LOAD_set_P_min(Load* load, REAL P_min) {
  if (load)
    load->P_min = P_min;
}

void LOAD_set_Q(Load* load, REAL Q, int t) { 
  if (load && t >= 0 && t < load->num_periods)
    load->Q[t] = Q;
}

int LOAD_set_flags(void* vload, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Load* load = (Load*)vload;
  int t;

  // Check load
  if (!load)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(load->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(load->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(load->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(load->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & LOAD_VAR_P) && (mask & LOAD_VAR_P)) {
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < load->num_periods; t++)
	load->index_P[t] = index+t;
    }
    (*flags_ptr) |= LOAD_VAR_P;
    index += load->num_periods;
  }
  return index;  
}

void LOAD_set_var_values(Load* load, Vec* values) {

  // Local vars
  int t;

  // No load
  if (!load)
    return;

  // Time loop
  for (t = 0; t < load->num_periods; t++) {

    if (load->vars & LOAD_VAR_P) // active power (p.u.)
      load->P[t] = VEC_get(values,load->index_P[t]);
  }
}

void LOAD_show(Load* load, int t) { 
  if (load)
    printf("load %d\t%d\n",
	   BUS_get_number(load->bus),
	   load->index);
}

