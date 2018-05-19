/** @file load.c
 *  @brief This file defines the Load data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/load.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct Load {

  // Bus
  Bus* bus;            /**< @brief Bus to which the load is connected */

  // Times
  int num_periods;     /**< @brief Number of time periods. */

  // Properties
  char name[LOAD_BUFFER_SIZE]; /**< @brief Load name */
  
  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Active power
  REAL* P;             /**< @brief Load active power consumption (p.u. system base power) */
  REAL* P_max;         /**< @brief Maximum load active power consumption (p.u.) */
  REAL* P_min;         /**< @brief Minimum load active power consumption (p.u.) */

  // Reactive power
  REAL* Q;              /**< @brief Load reactive power (p.u. system base power) */
  REAL* Q_max;         /**< @brief Maximum load reactive power consumption (p.u.) */
  REAL* Q_min;         /**< @brief Minimum load reactive power consumption (p.u.) */

  // Power factor
  REAL target_power_factor; /**< @brief Target power factor in (0,1] */

  // Utility
  REAL util_coeff_Q0;  /**< @brief Load utility coefficient (constant term, units of $/hr ) */
  REAL util_coeff_Q1;  /**< @brief Load utility coefficient (linear term, units of $/(hr p.u.) ) */
  REAL util_coeff_Q2;  /**< @brief Load utility coefficient (quadratic term, units of $/(hr p.u.^2) ) */

  // Indices
  int index;           /**< @brief Load index */
  int* index_P;        /**< @brief Active power index */
  int* index_Q;        /**< @brief Reactive power index */

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
      free(load->P_max);
      free(load->P_min);
      free(load->Q);
      free(load->Q_max);
      free(load->Q_min);
      free(load->index_P);
      free(load->index_Q);
      free(load->sens_P_u_bound);
      free(load->sens_P_l_bound);
      LOAD_set_bus(load,NULL);
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
      snprintf(load_array[i].name,(size_t)(LOAD_BUFFER_SIZE-1),"%d",i);
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

void LOAD_copy_from_load(Load* load, Load* other) {

  // Local variables
  int num_periods;

  // Check
  if (!load || !other)
    return;

  // Min num periods
  if (load->num_periods < other->num_periods)
    num_periods = load->num_periods;
  else
    num_periods = other->num_periods;

  // Bus
  // skip bue

  // Times
  // skip num periods

  // Properties
  strcpy(load->name,other->name);

  // Flags
  load->fixed = other->fixed;
  load->bounded = other->bounded;
  load->sparse = other->sparse;
  load->vars = other->vars;

  // Active power
  memcpy(load->P,other->P,num_periods*sizeof(REAL));
  memcpy(load->P_max,other->P_max,num_periods*sizeof(REAL));
  memcpy(load->P_min,other->P_min,num_periods*sizeof(REAL));

  // Reactive power
  memcpy(load->Q,other->Q,num_periods*sizeof(REAL));
  memcpy(load->Q_max,other->Q_max,num_periods*sizeof(REAL));
  memcpy(load->Q_min,other->Q_min,num_periods*sizeof(REAL));

  // Power factor
  load->target_power_factor = other->target_power_factor;

  // Utility
  load->util_coeff_Q0 = other->util_coeff_Q0;
  load->util_coeff_Q1 = other->util_coeff_Q1;
  load->util_coeff_Q2 = other->util_coeff_Q2;

  // Indices
  // skip index
  memcpy(load->index_P,other->index_P,num_periods*sizeof(int));
  memcpy(load->index_Q,other->index_Q,num_periods*sizeof(int));

  // Sensitivities
  memcpy(load->sens_P_u_bound,other->sens_P_u_bound,num_periods*sizeof(REAL));
  memcpy(load->sens_P_l_bound,other->sens_P_l_bound,num_periods*sizeof(REAL));

  // List
  // skip next
}

char LOAD_get_flags_vars(Load* load) {
  if (load)
    return load->vars;
  else
    return 0;
}

char LOAD_get_flags_fixed(Load* load) {
  if (load)
    return load->fixed;
  else
    return 0;
}

char LOAD_get_flags_bounded(Load* load) {
  if (load)
    return load->bounded;
  else
    return 0;
}

char LOAD_get_flags_sparse(Load* load) {
  if (load)
    return load->sparse;
  else
    return 0;
}

char* LOAD_get_name(Load* load) {
  if (load)
    return load->name;
  else
    return NULL;
}

int LOAD_get_num_periods(Load* load) {
  if (load)
    return load->num_periods;
  else
    return 0;
}

REAL LOAD_get_power_factor(Load* load, int t) {
  REAL S;
  if (load && t >= 0 && t < load->num_periods) {
    S = sqrt((load->P[t])*(load->P[t]) + (load->Q[t])*(load->Q[t]));
    if (S != 0.)
      return load->P[t]/S;
    else
      return 1.;
  }
  else
    return 1.;
}

REAL LOAD_get_target_power_factor(Load* load) {
  if (load)
    return load->target_power_factor;
  else
    return 1.;
}

REAL LOAD_get_sens_P_u_bound(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->sens_P_u_bound[t];
  else
    return 0;
}

REAL* LOAD_get_sens_P_u_bound_array(Load* load) {
  if (load)
    return load->sens_P_u_bound;
  else
    return NULL;
}

REAL LOAD_get_sens_P_l_bound(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->sens_P_l_bound[t];
  else
    return 0;
}

REAL* LOAD_get_sens_P_l_bound_array(Load* load) {
  if (load)
    return load->sens_P_l_bound;
  else
    return NULL;
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
    return -1;
}

int LOAD_get_index_P(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->index_P[t];
  else
    return -1;
}

int LOAD_get_index_Q(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->index_Q[t];
  else
    return -1;
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

REAL LOAD_get_P_max(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->P_max[t];
  else 
    return 0;
}

REAL LOAD_get_P_min(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->P_min[t];
  else 
    return 0;
}

REAL LOAD_get_Q(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->Q[t];
  else
    return 0;
}

REAL LOAD_get_Q_max(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->Q_max[t];
  else
    return 0;
}

REAL LOAD_get_Q_min(Load* load, int t) {
  if (load && t >= 0 && t < load->num_periods)
    return load->Q_min[t];
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
	if (load->bounded & LOAD_VAR_P)
	  VEC_set(values,load->index_P[t],load->P_max[t]);
	else
	  VEC_set(values,load->index_P[t],LOAD_INF_P);
	break;

      case LOWER_LIMITS:
	if (load->bounded & LOAD_VAR_P)
	  VEC_set(values,load->index_P[t],load->P_min[t]);
	else
	  VEC_set(values,load->index_P[t],-LOAD_INF_P);
	break;

      default:
	VEC_set(values,load->index_P[t],load->P[t]);
      }
    }

    if (load->vars & LOAD_VAR_Q) { // reactive power
      switch(code) {

      case UPPER_LIMITS:
	if (load->bounded & LOAD_VAR_Q)
	  VEC_set(values,load->index_Q[t],load->Q_max[t]);
	else
	  VEC_set(values,load->index_Q[t],LOAD_INF_Q);
	break;

      case LOWER_LIMITS:
	if (load->bounded & LOAD_VAR_Q)
	  VEC_set(values,load->index_Q[t],load->Q_min[t]);
	else
	  VEC_set(values,load->index_Q[t],-LOAD_INF_Q);
	break;

      default:
	VEC_set(values,load->index_Q[t],load->Q[t]);
      }
    }
  }
}

char* LOAD_get_var_info_string(Load* load, int index) {

  // Local variables
  char* info;

  //Check
  if (!load)
    return NULL;

  // Active power
  if ((load->vars & LOAD_VAR_P) &&
      index >= load->index_P[0] &&
      index <= load->index_P[load->num_periods-1]) {
    info = (char*)malloc(LOAD_BUFFER_SIZE*sizeof(char));
    snprintf(info,LOAD_BUFFER_SIZE*sizeof(char),
	     "load:%d:active power:%d",load->index,index-load->index_P[0]);
    return info;
  }

  // Reactive power
  if ((load->vars & LOAD_VAR_Q) &&
      index >= load->index_Q[0] &&
      index <= load->index_Q[load->num_periods-1]) {
    info = (char*)malloc(LOAD_BUFFER_SIZE*sizeof(char));
    snprintf(info,LOAD_BUFFER_SIZE*sizeof(char),
	     "load:%d:reactive power:%d",load->index,index-load->index_Q[0]);
    return info;
  }

  // Return
  return NULL;
}

int LOAD_get_num_vars(void* vload, unsigned char var, int t_start, int t_end) {

  // Local vars
  Load* load = (Load*)vload;
  int num_vars = 0;
  int dt;

  // Checks
  if (!load)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > load->num_periods-1)
    t_end = load->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & LOAD_VAR_P) && (load->vars & LOAD_VAR_P)) // active power
    num_vars += dt;
  if ((var & LOAD_VAR_Q) && (load->vars & LOAD_VAR_Q)) // reactive power
    num_vars += dt;
  return num_vars;
}

Vec* LOAD_get_var_indices(void* vload, unsigned char var, int t_start, int t_end) {

  // Local vars
  Load* load = (Load*)vload;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!load)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > load->num_periods-1)
    t_end = load->num_periods-1;

  // Indices
  indices = VEC_new(LOAD_get_num_vars(vload,var,t_start,t_end));
  if ((var & LOAD_VAR_P) && (load->vars & LOAD_VAR_P)) { // active power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,load->index_P[t]);
      offset++;
    }
  }
  if ((var & LOAD_VAR_Q) && (load->vars & LOAD_VAR_Q)) { // reactive power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,load->index_Q[t]);
      offset++;
    }
  }
  return indices;
}

char* LOAD_get_json_string(Load* load, char* output) {

  // Local variables
  char temp[LOAD_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No load
  if (!load)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*LOAD_JSON_BUFFER_SIZE*LOAD_NUM_JSON_FIELDS*load->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",load->index,FALSE);
  JSON_obj(temp,output,"bus",load->bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"num_periods",load->num_periods,FALSE);
  JSON_str(temp,output,"name",load->name,FALSE);
  JSON_array_float(temp,output,"P",load->P,load->num_periods,FALSE);
  JSON_array_float(temp,output,"P_max",load->P_max,load->num_periods,FALSE);
  JSON_array_float(temp,output,"P_min",load->P_min,load->num_periods,FALSE);
  JSON_array_float(temp,output,"Q",load->Q,load->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_max",load->Q_max,load->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_min",load->Q_min,load->num_periods,FALSE);
  JSON_float(temp,output,"target_power_factor",load->target_power_factor,FALSE);
  JSON_float(temp,output,"util_coeff_Q0",load->util_coeff_Q0,FALSE);
  JSON_float(temp,output,"util_coeff_Q1",load->util_coeff_Q1,FALSE);
  JSON_float(temp,output,"util_coeff_Q2",load->util_coeff_Q2,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL LOAD_has_flags(void* vload, char flag_type, unsigned char mask) {
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
  ARRAY_clear(load->name,char,LOAD_BUFFER_SIZE);
  
  load->bus = NULL;
  
  load->fixed = 0x00;
  load->bounded = 0x00;
  load->sparse = 0x00;
  load->vars = 0x00;
    
  load->util_coeff_Q0 = 0;
  load->util_coeff_Q1 = 20000.;
  load->util_coeff_Q2 = -100.;
  
  load->index = -1;

  ARRAY_zalloc(load->P,REAL,T);
  ARRAY_zalloc(load->P_max,REAL,T);
  ARRAY_zalloc(load->P_min,REAL,T);
  ARRAY_zalloc(load->Q,REAL,T);
  ARRAY_zalloc(load->Q_max,REAL,T);
  ARRAY_zalloc(load->Q_min,REAL,T);
  ARRAY_zalloc(load->index_P,int,T);
  ARRAY_zalloc(load->index_Q,int,T);
  ARRAY_zalloc(load->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(load->sens_P_l_bound,REAL,T);

  load->target_power_factor = 1.;
  
  load->next = NULL;
}

BOOL LOAD_is_equal(Load* load, Load* other) {
  return load == other;
}

BOOL LOAD_is_P_adjustable(Load* load) {
  int t;
  if (load) {
    for (t = 0; t < load->num_periods; t++) {
      if (load->P_min[t] < load->P_max[t])
	return TRUE;
    }
    return FALSE;
  }    
  else
    return FALSE;
}

Load* LOAD_list_add(Load* load_list, Load* load) {
  LIST_add(Load,load_list,load,next);
  return load_list;
}

Load* LOAD_list_del(Load* load_list, Load* load) {
  LIST_del(Load,load_list,load,next);
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

void LOAD_set_name(Load* load, char* name) {
  if (load)
    strncpy(load->name,name,(size_t)(LOAD_BUFFER_SIZE-1));
}

void LOAD_set_target_power_factor(Load* load, REAL pf) {
  if (load) {
    if (pf > 1.)
      load->target_power_factor = 1.;
    else if (pf < -1.)
      load->target_power_factor = -1.;
    else if (fabs(pf) < LOAD_MIN_TARGET_PF)
      load->target_power_factor = LOAD_MIN_TARGET_PF;
    else
      load->target_power_factor = pf;
  }
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
  Bus* old_bus;
  if (load) {
    old_bus = load->bus;
    load->bus = NULL;
    BUS_del_load(old_bus,load);
    load->bus = bus;
    BUS_add_load(load->bus,load);
  }
}

void LOAD_set_index(Load* load, int index) { 
  if (load)
    load->index = index;
}

void LOAD_set_P(Load* load, REAL P, int t) { 
  if (load && t >= 0 && t < load->num_periods)
    load->P[t] = P;
}

void LOAD_set_P_max(Load* load, REAL P_max, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->P_max[t] = P_max;
}

void LOAD_set_P_min(Load* load, REAL P_min, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->P_min[t] = P_min;
}

void LOAD_set_Q(Load* load, REAL Q, int t) { 
  if (load && t >= 0 && t < load->num_periods)
    load->Q[t] = Q;
}

void LOAD_set_Q_max(Load* load, REAL Q_max, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->Q_max[t] = Q_max;
}

void LOAD_set_Q_min(Load* load, REAL Q_min, int t) {
  if (load && t >= 0 && t < load->num_periods)
    load->Q_min[t] = Q_min;
}

int LOAD_set_flags(void* vload, char flag_type, unsigned char mask, int index) {

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
  if (!((*flags_ptr) & LOAD_VAR_P) && (mask & LOAD_VAR_P)) { // active power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < load->num_periods; t++)
	load->index_P[t] = index+t;
    }
    (*flags_ptr) |= LOAD_VAR_P;
    index += load->num_periods;
  }
  if (!((*flags_ptr) & LOAD_VAR_Q) && (mask & LOAD_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < load->num_periods; t++)
	load->index_Q[t] = index+t;
    }
    (*flags_ptr) |= LOAD_VAR_Q;
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

    if (load->vars & LOAD_VAR_Q) // reactive power (p.u.)
      load->Q[t] = VEC_get(values,load->index_Q[t]);
  }
}

void LOAD_show(Load* load, int t) { 
  if (load)
    printf("load %d\t%d\n",
	   BUS_get_number(load->bus),
	   load->index);
}

void LOAD_propagate_data_in_time(Load* load, int start, int end) {
  int t;
  if (load) {
    if (start < 0)
      start = 0;
    if (end > load->num_periods)
      end = load->num_periods;
    for (t = start+1; t < end; t++) {
      load->P[t] = load->P[start];
      load->P_max[t] = load->P_max[start];
      load->P_min[t] = load->P_min[start];
      load->Q[t] = load->Q[start];
      load->Q_max[t] = load->Q_max[start];
      load->Q_min[t] = load->Q_min[start];
    }
  }
}

