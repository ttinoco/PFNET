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

struct Load {

  // Bus
  Bus* bus;            /**< @brief Bus to which the load is connected */

  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Active power
  REAL P;              /**< @brief Load active power consumption (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum load active power consumption (p.u.) */
  REAL P_min;          /**< @brief Minimum load active power consumption (p.u.) */

  // Reactive power
  REAL Q;              /**< @brief Load reactive power (p.u. system base power) */

  // Utility
  REAL util_coeff_Q0;  /**< @brief Load utility coefficient (constant term, units of $/hr ) */
  REAL util_coeff_Q1;  /**< @brief Load utility coefficient (linear term, units of $/(hr p.u.) ) */
  REAL util_coeff_Q2;  /**< @brief Load utility coefficient (quadratic term, units of $/(hr p.u.^2) ) */

  // Indices
  int index;           /**< @brief Load index */
  int index_P;         /**< @brief Active power index */

  // Sensitivities
  REAL sens_P_u_bound;  /**< @brief Sensitivity of active power upper bound */
  REAL sens_P_l_bound;  /**< @brief Sensitivity of active power lower bound */

  // List
  Load* next;           /**< @brief List of loads connected to a bus */
};

void* LOAD_array_get(void* load, int index) { 
  if (load) 
    return (void*)&(((Load*)load)[index]);
  else
    return NULL;
}

Load* LOAD_array_new(int num) { 
  int i;
  Load* load = (Load*)malloc(sizeof(Load)*num);
  for (i = 0; i < num; i++) {
    LOAD_init(&(load[i]));
    LOAD_set_index(&(load[i]),i);
  }
  return load;
}

void LOAD_array_show(Load* load, int num) { 
  int i;
  for (i = 0; i < num; i++) 
    LOAD_show(&(load[i]));
}

void LOAD_clear_sensitivities(Load* load) {
  if (load) {
    load->sens_P_u_bound = 0;
    load->sens_P_l_bound = 0;
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

REAL LOAD_get_sens_P_u_bound(Load* load) {
  if (load)
    return load->sens_P_u_bound;
  else
    return 0;
}

REAL LOAD_get_sens_P_l_bound(Load* load) {
  if (load)
    return load->sens_P_l_bound;
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

REAL LOAD_get_P_util(Load* load) {
  if (load)
    return LOAD_get_P_util_at(load,load->P); // $/hr
  else
    return 0;
}

REAL LOAD_get_P_util_at(Load* load, REAL P) {
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

int LOAD_get_index_P(Load* load) {
  if (load)
    return load->index_P;
  else
    return 0;
}

Load* LOAD_get_next(Load* load) {
  if (load)
    return load->next;
  else
    return NULL;
}

REAL LOAD_get_P(Load* load) {
  if (load)
    return load->P;
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

REAL LOAD_get_Q(Load* load) {
  if (load)
    return load->Q;
  else
    return 0;
}

void LOAD_get_var_values(Load* load, Vec* values, int code) {
  
  if (!load)
    return;

  if (load->vars & LOAD_VAR_P) { // active power
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,load->index_P,load->P_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,load->index_P,load->P_min);
      break;
    default:
      VEC_set(values,load->index_P,load->P);
    }
  }
}

Vec* LOAD_get_var_indices(void* vload, char var) {
  Load* load = (Load*)vload;
  Vec* indices;
  if (!load)
    return NULL;
  if (var == LOAD_VAR_P) {
    indices = VEC_new(1);
    VEC_set(indices,0,load->index_P);
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

void LOAD_init(Load* load) { 
  if (load) {
    
    load->bus = NULL;
    
    load->fixed = 0x00;
    load->bounded = 0x00;
    load->sparse = 0x00;
    load->vars = 0x00;
    
    load->P = 0;
    load->P_max = 0;
    load->P_min = 0;

    load->Q = 0;

    load->util_coeff_Q0 = 0;
    load->util_coeff_Q1 = 20000.;
    load->util_coeff_Q2 = -100.;

    load->index = 0;
    load->index_P = 0;

    load->sens_P_u_bound = 0;
    load->sens_P_l_bound = 0;

    load->next = NULL;
  }
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

Load* LOAD_new(void) { 
  Load* load = (Load*)malloc(sizeof(Load));
  LOAD_init(load);
  return load;
}

void LOAD_set_sens_P_u_bound(Load* load, REAL value) {
  if (load)
    load->sens_P_u_bound = value;
}

void LOAD_set_sens_P_l_bound(Load* load, REAL value) {
  if (load)
    load->sens_P_l_bound = value;
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

void LOAD_set_P(Load* load, REAL P) { 
  if (load)
    load->P = P;
}

void LOAD_set_P_max(Load* load, REAL P_max) {
  if (load)
    load->P_max = P_max;
}

void LOAD_set_P_min(Load* load, REAL P_min) {
  if (load)
    load->P_min = P_min;
}

void LOAD_set_Q(Load* load, REAL Q) { 
  if (load)
    load->Q = Q;
}

int LOAD_set_flags(void* vload, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Load* load = (Load*)vload;

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
    if (flag_type == FLAG_VARS)
      load->index_P = index;
    (*flags_ptr) |= LOAD_VAR_P;
    index++;
  }
  return index;  
}

void LOAD_set_var_values(Load* load, Vec* values) {
  
  if (!load)
    return;
  if (load->vars & LOAD_VAR_P) // active power (p.u.)
    load->P = VEC_get(values,load->index_P);
}

void LOAD_show(Load* load) { 
  if (load)
    printf("load %d\t%d\n",BUS_get_number(load->bus),load->index);
}

