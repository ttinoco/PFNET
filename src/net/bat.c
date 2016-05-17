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

struct Bat {

  // Bus
  Bus* bus;            /**< @brief Bus to which the battery is connected */
  
  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Charging power
  REAL P;              /**< @brief Battery charging power (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum charging power (p.u.) */
  REAL P_min;          /**< @brief Minimum charging power (p.u.) */

  // Efficiencies
  REAL eta_c;          /**< @brief Battery charging efficiency (unitless) */
  REAL eta_d;          /**< @brief Battery discharging efficiency (unitless) */

  // Energy level
  REAL E;              /**< @brief Battery energy level (p.u. system base power times time unit) */
  REAL E_max;          /**< @brief Maximum energy level (p.u. times time unit) */
  
  // Indices
  int index;           /**< @brief Battery index */
  int index_Pc;        /**< @brief charging power index */
  int index_Pd;        /**< @brief discharging power index */
  int index_E;         /**< @brief energy level index */

  // List
  Bat* next;           /**< @brief List of batteries connected to a bus */
};

void* BAT_array_get(void* bat, int index) { 
  if (bat) 
    return (void*)&(((Bat*)bat)[index]);
  else
    return NULL;
}

Bat* BAT_array_new(int num) { 
  int i;
  Bat* bat = (Bat*)malloc(sizeof(Bat)*num);
  for (i = 0; i < num; i++) {
    BAT_init(&(bat[i]));
    BAT_set_index(&(bat[i]),i);
  }
  return bat;
}

void BAT_array_show(Bat* bat, int num) { 
  int i;
  for (i = 0; i < num; i++) 
    BAT_show(&(bat[i]));
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

int BAT_get_index_Pc(Bat* bat) {
  if (bat)
    return bat->index_Pc;
  else
    return 0;
}

int BAT_get_index_Pd(Bat* bat) {
  if (bat)
    return bat->index_Pd;
  else
    return 0;
}

int BAT_get_index_E(Bat* bat) {
  if (bat)
    return bat->index_E;
  else
    return 0;
}

Bat* BAT_get_next(Bat* bat) {
  if (bat)
    return bat->next;
  else
    return NULL;
}

REAL BAT_get_P(Bat* bat) {
  if (bat)
    return bat->P;
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

REAL BAT_get_E(Bat* bat) {
  if (bat)
    return bat->E;
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
  
  if (!bat)
    return;

  // Charging power
  if (bat->vars & BAT_VAR_P) {
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,bat->index_Pc,bat->P_max);
      VEC_set(values,bat->index_Pd,-bat->P_min);
      break;
    case LOWER_LIMITS:
      VEC_set(values,bat->index_Pc,0.);
      VEC_set(values,bat->index_Pd,0.);
      break;
    default:
      if (bat->P >= 0) {
	VEC_set(values,bat->index_Pc,bat->P);
	VEC_set(values,bat->index_Pd,0.);
      }
      else {
	VEC_set(values,bat->index_Pc,0.);
	VEC_set(values,bat->index_Pd,-bat->P);
      }
    }
  }

  // Energy level
  if (bat->vars & BAT_VAR_E) {
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,bat->index_E,bat->E_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,bat->index_E,0.);
      break;
    default:
      VEC_set(values,bat->index_E,bat->E);
    }
  }
}

Vec* BAT_get_var_indices(void* vbat, char var) {
  Bat* bat = (Bat*)vbat;
  Vec* indices;
  if (!bat)
    return NULL;
  if (var == BAT_VAR_P) {
    indices = VEC_new(2);
    VEC_set(indices,0,bat->index_Pc);
    VEC_set(indices,1,bat->index_Pd);
    return indices;
  }
  if (var == BAT_VAR_E) {
    indices = VEC_new(1);
    VEC_set(indices,0,bat->index_E);
    return indices;
  }
  return NULL;
}

BOOL BAT_has_flags(void* vbat, char flag_type, char mask) {
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

void BAT_init(Bat* bat) { 
  if (bat) {
    
    bat->bus = NULL;
    
    bat->fixed = 0x00;
    bat->bounded = 0x00;
    bat->sparse = 0x00;
    bat->vars = 0x00;
    
    bat->P = 0;
    bat->P_max = 0;
    bat->P_min = 0;

    bat->E = 0;
    bat->E_max = 0;

    bat->eta_c = 1.;
    bat->eta_d = 1.;
    
    bat->index = 0;
    bat->index_Pc = 0;
    bat->index_Pd = 0;
    bat->index_E = 0;

    bat->next = NULL;
  }
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

Bat* BAT_new(void) { 
  Bat* bat = (Bat*)malloc(sizeof(Bat));
  BAT_init(bat);
  return bat;
}

void BAT_set_bus(Bat* bat, Bus* bus) { 
  if (bat)
    bat->bus = (Bus*)bus;
}

void BAT_set_index(Bat* bat, int index) { 
  if (bat)
    bat->index = index;
}

void BAT_set_P(Bat* bat, REAL P) { 
  if (bat)
    bat->P = P;
}

void BAT_set_P_max(Bat* bat, REAL P_max) {
  if (bat)
    bat->P_max = P_max;
}

void BAT_set_P_min(Bat* bat, REAL P_min) {
  if (bat)
    bat->P_min = P_min;
}

void BAT_set_E(Bat* bat, REAL E) { 
  if (bat)
    bat->E = E;
}

void BAT_set_E_max(Bat* bat, REAL E_max) {
  if (bat)
    bat->E_max = E_max;
}

void BAT_set_eta_c(Bat* bat, REAL eta_c) { 
  if (bat)
    bat->eta_c = eta_c;
}

void BAT_set_eta_d(Bat* bat, REAL eta_d) { 
  if (bat)
    bat->eta_d = eta_d;
}

int BAT_set_flags(void* vbat, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Bat* bat = (Bat*)vbat;

  // Check bat
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
      bat->index_Pc = index;
      bat->index_Pd = index+1;
    }
    (*flags_ptr) |= BAT_VAR_P;
    index += 2;
  }
  if (!((*flags_ptr) & BAT_VAR_E) && (mask & BAT_VAR_E)) { // energy level
    if (flag_type == FLAG_VARS)
      bat->index_E = index;
    (*flags_ptr) |= BAT_VAR_E;
    index++;
  }
  return index;  
}

void BAT_set_var_values(Bat* bat, Vec* values) {
  
  if (!bat)
    return;
  if (bat->vars & BAT_VAR_P) // charging/discharging power
    bat->P = VEC_get(values,bat->index_Pc)-VEC_get(values,bat->index_Pd);
  if (bat->vars & BAT_VAR_E) // energy level
    bat->E = VEC_get(values,bat->index_E);
}

void BAT_show(Bat* bat) { 
  if (bat)
    printf("bat %d\t%d\n",BUS_get_number(bat->bus),bat->index);
}

