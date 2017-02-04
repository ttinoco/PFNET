/** @file bat.c
 *  @brief This file defines the Bat data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/bat.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>

struct Bat {

// Bus
  Bus* bus;            /**< @brief Bus to which the battery is connected */
 
  // Times
  int num_periods;   /**< @brief Number of time periods. */ 

 // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Charging power
  REAL* P;             /**< @brief Battery charging power during a time period (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum charging power (p.u.) */
  REAL P_min;          /**< @brief Minimum charging power (p.u.) */

  // Efficiencies
  REAL eta_c;          /**< @brief Battery charging efficiency (unitless) */
  REAL eta_d;          /**< @brief Battery discharging efficiency (unitless) */

  // Energy level
  REAL* E;             /**< @brief Battery energy level at the beginning of a period (p.u. system base power times time unit) */
  REAL E_init;         /**< @brief Initial battery energy level (p.u. system base power times time unit) */
  REAL E_final;        /**< @brief Battery energy level at the end of last period (p.u. system base power times time unit) */
  REAL E_max;          /**< @brief Maximum energy level (p.u. times time unit) */
  
  // Indices
  int index;           /**< @brief Battery index */
  int* index_Pc;       /**< @brief charging power index */
  int* index_Pd;       /**< @brief discharging power index */
  int* index_E;        /**< @brief energy level index */

  // List
  Bat* next;           /**< @brief List of batteries connected to a bus */
};

void* BAT_array_get(void* bat_array, int index) { 
  if (bat_array) 
    return (void*)&(((Bat*)bat_array)[index]);
  else
    return NULL;
}

void BAT_array_del(Bat* bat_array, int size) {
  int i;
  Bat* bat;
  if (bat_array) {
    for (i = 0; i < size; i++) {
      bat = &(bat_array[i]);
      free(bat->P);
      free(bat->E);
      free(bat->index_Pc);
      free(bat->index_Pd);
      free(bat->index_E);
    }
    free(bat_array);
  }  
}

Bat* BAT_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    Bat* bat_array = (Bat*)malloc(sizeof(Bat)*size);
    for (i = 0; i < size; i++) {
      BAT_init(&(bat_array[i]),num_periods);
      BAT_set_index(&(bat_array[i]),i);
    }
    return bat_array;
  }
  else
    return NULL;
}

void BAT_array_show(Bat* bat_array, int size, int t) { 
  int i;
  if (bat_array) {
    for (i = 0; i < size; i++) 
      BAT_show(&(bat_array[i]),t);
  }
}

void BAT_clear_flags(Bat* bat, char flag_type) {
  if (bat) {
    if (flag_type == FLAG_VARS)
      bat->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      bat->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      bat->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      bat->sparse = 0x00;
  }
}

int BAT_get_num_periods(Bat* bat) {
  if (bat)
    return bat->num_periods;
  else
    return 0;
}

char BAT_get_obj_type(void* bat) {
  if (bat)
    return OBJ_BAT;
  else
    return OBJ_UNKNOWN;
}

Bus* BAT_get_bus(Bat* bat) {
  if (bat)
    return bat->bus;
  else
    return NULL;
}

int BAT_get_index(Bat* bat) {
  if (bat)
    return bat->index;
  else
    return 0;
}

int BAT_get_index_Pc(Bat* bat, int t) {
  if (bat && t >= 0 && t < bat->num_periods)
    return bat->index_Pc[t];
  else
    return 0;
}

int BAT_get_index_Pd(Bat* bat, int t) {
  if (bat && t >= 0 && t < bat->num_periods)
    return bat->index_Pd[t];
  else
    return 0;
}

int BAT_get_index_E(Bat* bat, int t) {
  if (bat && t >= 0 && t < bat->num_periods)
    return bat->index_E[t];
  else
    return 0;
}

Bat* BAT_get_next(Bat* bat) {
  if (bat)
    return bat->next;
  else
    return NULL;
}

REAL BAT_get_P(Bat* bat, int t) {
  if (bat && t >= 0 && t < bat->num_periods)
    return bat->P[t];
  else
    return 0;
}

REAL BAT_get_P_max(Bat* bat) {
  if (bat)
    return bat->P_max;
  else 
    return 0;
}

REAL BAT_get_P_min(Bat* bat) {
  if (bat)
    return bat->P_min;
  else 
    return 0;
}

REAL BAT_get_E(Bat* bat, int t) {
  if (bat && t >= 0 && t < bat->num_periods)
    return bat->E[t];
  else
    return 0;
}

REAL BAT_get_E_init(Bat* bat) {
  if (bat)
    return bat->E_init;
  else 
    return 0;
}

REAL BAT_get_E_final(Bat* bat) {
  if (bat)
    return bat->E_final;
  else 
    return 0;
}

REAL BAT_get_E_max(Bat* bat) {
  if (bat)
    return bat->E_max;
  else 
    return 0;
}

REAL BAT_get_eta_c(Bat* bat) {
  if (bat)
    return bat->eta_c;
  else
    return 0;
}

REAL BAT_get_eta_d(Bat* bat) {
  if (bat)
    return bat->eta_d;
  else
    return 0;
}

void BAT_get_var_values(Bat* bat, Vec* values, int code) {
 
  // Local vars
  int t;
 
  // No bat
  if (!bat)
    return;

  // Time loop
  for (t = 0; t < bat->num_periods; t++) {

    // Charging power
    if (bat->vars & BAT_VAR_P) {
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,bat->index_Pc[t],bat->P_max);
	VEC_set(values,bat->index_Pd[t],-bat->P_min);
	break;
      case LOWER_LIMITS:
	VEC_set(values,bat->index_Pc[t],0.);
	VEC_set(values,bat->index_Pd[t],0.);
	break;
      default:
	if (bat->P[t] >= 0) {
	  VEC_set(values,bat->index_Pc[t],bat->P[t]);
	  VEC_set(values,bat->index_Pd[t],0.);
	}
	else {	  
	  VEC_set(values,bat->index_Pc[t],0.);
	  VEC_set(values,bat->index_Pd[t],-bat->P[t]);
	}
      }
    }

    // Energy level
    if (bat->vars & BAT_VAR_E) {
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,bat->index_E[t],bat->E_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,bat->index_E[t],0.);
	break;
      default:
	VEC_set(values,bat->index_E[t],bat->E[t]);
      }
    }
  }
}

int BAT_get_num_vars(void* vbat, unsigned char var, int t_start, int t_end) {

  // Local vars
  Bat* bat = (Bat*)vbat;
  int num_vars = 0;
  int dt;

  // Checks
  if (!bat)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > bat->num_periods-1)
    t_end = bat->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & BAT_VAR_P) && (bat->vars & BAT_VAR_P))
    num_vars += 2*dt;
  if ((var & BAT_VAR_E) && (bat->vars & BAT_VAR_E)) 
    num_vars += dt;
  return num_vars;
}

Vec* BAT_get_var_indices(void* vbat, unsigned char var, int t_start, int t_end) {

  // Local vars
  Bat* bat = (Bat*)vbat;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!bat)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > bat->num_periods-1)
    t_end = bat->num_periods-1;

  // Indices
  indices = VEC_new(BAT_get_num_vars(vbat,var,t_start,t_end));
  if ((var & BAT_VAR_P) && (bat->vars & BAT_VAR_P)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,bat->index_Pc[t]);
      VEC_set(indices,offset+1,bat->index_Pd[t]);
      offset += 2;
    }
  }
  if ((var & BAT_VAR_E) && (bat->vars & BAT_VAR_E)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,bat->index_E[t]);
      offset++;
    }
  }
  return indices;
}

BOOL BAT_has_flags(void* vbat, char flag_type, unsigned char mask) {
  Bat* bat = (Bat*)vbat;
  if (bat) {
    if (flag_type == FLAG_VARS)
      return (bat->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (bat->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (bat->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (bat->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL BAT_has_properties(void* vbat, char prop) {
  Bat* bat = (Bat*)vbat;
  if (!bat)
    return FALSE;
  return TRUE;
}

void BAT_init(Bat* bat, int num_periods) { 

  // Local vars
  int T;
  
  // No gen
  if (!bat)
    return;
  
  T = num_periods;
  bat->num_periods = num_periods;
  
  bat->bus = NULL;
  
  bat->fixed = 0x00;
  bat->bounded = 0x00;
  bat->sparse = 0x00;
  bat->vars = 0x00;
   
  bat->P_max = 0;
  bat->P_min = 0;
  
  bat->E_init = 0;
  bat->E_final = 0;
  bat->E_max = 0;
  
  bat->eta_c = 1.;
  bat->eta_d = 1.;
  
  bat->index = 0;

  ARRAY_zalloc(bat->P,REAL,T);
  ARRAY_zalloc(bat->E,REAL,T);
  ARRAY_zalloc(bat->index_Pc,int,T);
  ARRAY_zalloc(bat->index_Pd,int,T);
  ARRAY_zalloc(bat->index_E,int,T);
  
  bat->next = NULL;
}

Bat* BAT_list_add(Bat* bat_list, Bat* bat) {
  LIST_add(Bat,bat_list,bat,next);
  return bat_list;
}

int BAT_list_len(Bat* bat_list) {
  int len;
  LIST_len(Bat,bat_list,next,len);
  return len;
}

Bat* BAT_new(int num_periods) {
  if (num_periods > 0) {
    Bat* bat = (Bat*)malloc(sizeof(Bat));
    BAT_init(bat,num_periods);
    return bat;
  }
  else
    return NULL;
}

void BAT_set_bus(Bat* bat, Bus* bus) { 
  if (bat)
    bat->bus = (Bus*)bus;
}

void BAT_set_index(Bat* bat, int index) { 
  if (bat)
    bat->index = index;
}

void BAT_set_P(Bat* bat, REAL P, int t) { 
  if (bat && t >= 0 && t < bat->num_periods)
    bat->P[t] = P;
}

void BAT_set_P_max(Bat* bat, REAL P_max) {
  if (bat)
    bat->P_max = P_max;
}

void BAT_set_P_min(Bat* bat, REAL P_min) {
  if (bat)
    bat->P_min = P_min;
}

void BAT_set_E(Bat* bat, REAL E, int t) { 
  if (bat && t >= 0 && t < bat->num_periods)
    bat->E[t] = E;
}

void BAT_set_E_init(Bat* bat, REAL E) {
  if (bat)
    bat->E_init = E;
}

void BAT_set_E_final(Bat* bat, REAL E) {
  if (bat)
    bat->E_final = E;
}

void BAT_set_E_max(Bat* bat, REAL E) {
  if (bat)
    bat->E_max = E;
}

void BAT_set_eta_c(Bat* bat, REAL eta_c) { 
  if (bat)
    bat->eta_c = eta_c;
}

void BAT_set_eta_d(Bat* bat, REAL eta_d) { 
  if (bat)
    bat->eta_d = eta_d;
}

int BAT_set_flags(void* vbat, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Bat* bat = (Bat*)vbat;
  int t;

  // No bat
  if (!bat)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(bat->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(bat->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(bat->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(bat->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & BAT_VAR_P) && (mask & BAT_VAR_P)) { // charging/discharging power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < bat->num_periods; t++) {
	bat->index_Pc[t] = index+2*t;
	bat->index_Pd[t] = index+2*t+1;
      }
    }
    (*flags_ptr) |= BAT_VAR_P;
    index += 2*bat->num_periods;
  }
  if (!((*flags_ptr) & BAT_VAR_E) && (mask & BAT_VAR_E)) { // energy level
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < bat->num_periods; t++)
	bat->index_E[t] = index+t;
    }
    (*flags_ptr) |= BAT_VAR_E;
    index += bat->num_periods;
  }
  return index;  
}

void BAT_set_var_values(Bat* bat, Vec* values) {

  // Local vars
  int t;
  
  // No bat
  if (!bat)
    return;
  
  // Time loop
  for (t = 0; t < bat->num_periods; t++) {

    if (bat->vars & BAT_VAR_P) // charging/discharging power
      bat->P[t] = VEC_get(values,bat->index_Pc[t])-VEC_get(values,bat->index_Pd[t]);
    if (bat->vars & BAT_VAR_E) // energy level
      bat->E[t] = VEC_get(values,bat->index_E[t]);
  }
}

void BAT_show(Bat* bat, int t) { 
  if (bat)
    printf("bat %d\t%d\n",
	   BUS_get_number(bat->bus),
	   bat->index);
}

void BAT_propagate_data_in_time(Bat* bat) {
  int t;
  if (bat) {
    for (t = 1; t < bat->num_periods; t++) {
      bat->P[t] = bat->P[0];
      bat->E[t] = bat->E[0];
    }
  }
}

