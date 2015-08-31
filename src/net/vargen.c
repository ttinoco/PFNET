/** @file vargen.c
 *  @brief This file defines the Vargen data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/vargen.h>
#include <pfnet/bus.h>

struct Vargen {

  // Bus
  Bus* bus;            /**< @brief Bus to which variable generator is connected */

  // Properties
  char type;         /**< @brief Variable generator type */
  
  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Active power
  REAL P;              /**< @brief Variable generator active power (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum variable generator active power (p.u.) */
  REAL P_std;          /**< @brief Standard deviation of active power (p.u. system base power) */

  // Indices
  int index;           /**< @brief Generator index */
  int index_P;         /**< @brief Active power index */

  // List
  Vargen* next;        /**< @brief List of variable generators connected to a bus */
};

void* VARGEN_array_get(void* gen, int index) {
  if (gen)
    return (void*)&(((Vargen*)gen)[index]);
  else 
    return NULL;
}

Vargen* VARGEN_array_new(int num) {
  int i;
  Vargen* gen = (Vargen*)malloc(sizeof(Vargen)*num);
  for (i = 0; i < num; i++) {
    VARGEN_init(&(gen[i]));
    VARGEN_set_index(&(gen[i]),i);
  }
  return gen;
}

void VARGEN_array_show(Vargen* gen, int num) {
  int i;
  if (gen) {
    for (i = 0; i < num; i++) 
      VARGEN_show(&(gen[i]));
  }
}

void VARGEN_clear_flags(Vargen* gen, char flag_type) {
  if (gen) {
    if (flag_type == FLAG_VARS)
      gen->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      gen->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      gen->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      gen->sparse = 0x00;
  }
}

void* VARGEN_get_bus(Vargen* gen) {
  if (gen)
    return (void*)gen->bus;
  else
    return NULL;
}

int VARGEN_get_index(Vargen* gen) {
  if (gen)
    return gen->index;
  else
    return 0;
}

int VARGEN_get_index_P(Vargen* gen) {
  if (gen)
    return gen->index_P;
  else
    return 0;
}

Vargen* VARGEN_get_next(Vargen* gen) {
  if (gen)
    return gen->next;
  else
    return NULL;
}

REAL VARGEN_get_P(Vargen* gen) {
  if (gen)
    return gen->P;
  else 
    return 0;
}

REAL VARGEN_get_P_max(Vargen* gen) {
  if (gen)
    return gen->P_max;
  else 
    return 0;
}

REAL VARGEN_get_P_std(Vargen* gen) {
  if (gen)
    return gen->P_std;
  else 
    return 0;
}

void VARGEN_get_var_values(Vargen* gen, Vec* values, int code) {  

  if (!gen)
    return;

  if (gen->vars & VARGEN_VAR_P) { // active power
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,gen->index_P,gen->P_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,gen->index_P,0.);
      break;
    default:
      VEC_set(values,gen->index_P,gen->P);
    }
  }
}

int VARGEN_get_var_index(void* vgen, char var) {
  Vargen* gen = (Vargen*)vgen;
  if (!gen)
    return 0;
  if (var == VARGEN_VAR_P)
    return gen->index_P;
  return 0;
}

BOOL VARGEN_has_flags(void* vgen, char flag_type, char mask) {
  Vargen* gen = (Vargen*)vgen;
  if (gen) {
    if (flag_type == FLAG_VARS)
      return (gen->vars & mask);
    else if (flag_type == FLAG_BOUNDED)
      return (gen->bounded & mask);
    else if (flag_type == FLAG_FIXED)
      return (gen->fixed & mask);
    else if (flag_type == FLAG_SPARSE)
      return (gen->sparse & mask);
    return FALSE;
  }
  else
    return FALSE;
}

BOOL VARGEN_has_properties(void* vgen, char prop) {
  Vargen* gen = (Vargen*)vgen;
  if (!gen)
    return FALSE;
  return TRUE;
}

void VARGEN_init(Vargen* gen) {
  gen->bus = NULL;
  gen->type = VARGEN_TYPE_WIND;
  gen->fixed = 0x00;
  gen->bounded = 0x00;
  gen->sparse = 0x00;
  gen->vars = 0x00;
  gen->P = 0;
  gen->P_max = 0;
  gen->P_std = 0;
  gen->index = 0;
  gen->index_P = 0;
  gen->next = NULL;
}

BOOL VARGEN_is_wind_farm(Vargen* gen) {
  if (gen)
    return gen->type == VARGEN_TYPE_WIND;
  else
    return FALSE;
}

BOOL VARGEN_is_solar_plant(Vargen* gen) {
  if (gen)
    return gen->type == VARGEN_TYPE_SOLAR;
  else
    return FALSE;
}

Vargen* VARGEN_list_add(Vargen *gen_list, Vargen* gen) {
  LIST_add(gen_list,gen,next);
  return gen_list;
}

int VARGEN_list_len(Vargen* gen_list) {
  int len;
  LIST_len(Vargen,gen_list,next,len);
  return len;
}

Vargen* VARGEN_new(void) {
  Vargen* gen = (Vargen*)malloc(sizeof(Vargen));
  VARGEN_init(gen);
  return gen;
}

void VARGEN_set_bus(Vargen* gen, void* bus) {
  gen->bus = (Bus*)bus;
}

void VARGEN_set_index(Vargen* gen, int index) {
  gen->index = index;
}

void VARGEN_set_P(Vargen* gen, REAL P) {
  gen->P = P;
}

void VARGEN_set_P_max(Vargen* gen, REAL P_max) {
  gen->P_max = P_max;
}

void VARGEN_set_P_std(Vargen* gen, REAL P_std) {
  gen->P_std = P_std;
}

int VARGEN_set_flags(void* vgen, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Vargen* gen = (Vargen*)vgen;

  // Check gen
  if (!gen)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(gen->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(gen->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(gen->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(gen->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & VARGEN_VAR_P) && (mask & VARGEN_VAR_P)) {
    if (flag_type == FLAG_VARS)
      gen->index_P = index;
    (*flags_ptr) |= VARGEN_VAR_P;
    index++;
  }
  return index;  
}

void VARGEN_set_var_values(Vargen* gen, Vec* values) {
  
  if (!gen)
    return;
  if (gen->vars & VARGEN_VAR_P) // active power (p.u.)
    gen->P = VEC_get(values,gen->index_P);
}

void VARGEN_show(Vargen* gen) {
  printf("vargen %d\t%d\n",
	 BUS_get_number(VARGEN_get_bus(gen)),
	 VARGEN_is_wind_farm(gen));
}

