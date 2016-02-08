/** @file vargen.c
 *  @brief This file defines the Vargen data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
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

  // Reactive power
  REAL Q;              /**< @brief Variable generator reactive power (p.u. system base power) */
  REAL Q_max;          /**< @brief Maximum variable generator reactive power (p.u.) */
  REAL Q_min;          /**< @brief Minimum variable generator reactive power (p.u.) */ 
  
  // Indices
  int index;           /**< @brief Generator index */
  int index_P;         /**< @brief Active power index */
  int index_Q;         /**< @brief Reactive power index */

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
    if (flag_type == FLAG_VARS)         // variables
      gen->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED) // bounded
      gen->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)   // fixed
      gen->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)  // sparse
      gen->sparse = 0x00;
  }
}

Bus* VARGEN_get_bus(Vargen* gen) {
  if (gen)
    return gen->bus;
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

int VARGEN_get_index_Q(Vargen* gen) {
  if (gen)
    return gen->index_Q;
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

REAL VARGEN_get_Q(Vargen* gen) {
  if (gen)
    return gen->Q;
  else 
    return 0;
}

REAL VARGEN_get_Q_max(Vargen* gen) {
  if (gen)
    return gen->Q_max;
  else 
    return 0;
}

REAL VARGEN_get_Q_min(Vargen* gen) {
  if (gen)
    return gen->Q_min;
  else 
    return 0;
}

void VARGEN_get_var_values(Vargen* gen, Vec* values, int code) {  

  if (!gen)
    return;

  // active power
  if (gen->vars & VARGEN_VAR_P) {
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

  // reactive power
  if (gen->vars & VARGEN_VAR_Q) {
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,gen->index_Q,gen->Q_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,gen->index_Q,gen->Q_min);
      break;
    default:
      VEC_set(values,gen->index_Q,gen->Q);
    }
  }
}

int VARGEN_get_var_index(void* vgen, char var) {
  Vargen* gen = (Vargen*)vgen;
  if (!gen)
    return 0;
  if (var == VARGEN_VAR_P)
    return gen->index_P;
  if (var == VARGEN_VAR_Q)
    return gen->index_Q;
  return 0;
}

BOOL VARGEN_has_flags(void* vgen, char flag_type, char mask) {
  Vargen* gen = (Vargen*)vgen;
  if (gen) {
    if (flag_type == FLAG_VARS)         // variables
      return (gen->vars & mask);
    else if (flag_type == FLAG_BOUNDED) // bounded
      return (gen->bounded & mask);
    else if (flag_type == FLAG_FIXED)   // fixed
      return (gen->fixed & mask);
    else if (flag_type == FLAG_SPARSE)  // sparse
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
  gen->Q = 0;
  gen->Q_max = 0;
  gen->Q_min = 0;
  gen->index = 0;
  gen->index_P = 0;
  gen->index_Q = 0;
  gen->next = NULL;
}

BOOL VARGEN_is_wind(Vargen* gen) {
  if (gen)
    return gen->type == VARGEN_TYPE_WIND;
  else
    return FALSE;
}

BOOL VARGEN_is_solar(Vargen* gen) {
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

void VARGEN_set_bus(Vargen* gen, Bus* bus) {
  if (gen)
    gen->bus = (Bus*)bus;
}

void VARGEN_set_index(Vargen* gen, int index) {
  if (gen)
    gen->index = index;
}

void VARGEN_set_P(Vargen* gen, REAL P) {
  if (gen)
    gen->P = P;
}

void VARGEN_set_P_max(Vargen* gen, REAL P_max) {
  if (gen)
    gen->P_max = P_max;
}

void VARGEN_set_P_std(Vargen* gen, REAL P_std) {
  if (gen)
    gen->P_std = P_std;
}

void VARGEN_set_Q(Vargen* gen, REAL Q) {
  if (gen)
    gen->Q = Q;
}

void VARGEN_set_Q_max(Vargen* gen, REAL Q) {
  if (gen)
    gen->Q_max = Q;
}

void VARGEN_set_Q_min(Vargen* gen, REAL Q) {
  if (gen)
    gen->Q_min = Q;
}

int VARGEN_set_flags(void* vgen, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Vargen* gen = (Vargen*)vgen;

  // Check gen
  if (!gen)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)         // variables
    flags_ptr = &(gen->vars);
  else if (flag_type == FLAG_FIXED)   // fixed
    flags_ptr = &(gen->fixed);
  else if (flag_type == FLAG_BOUNDED) // bounded
    flags_ptr = &(gen->bounded);
  else if (flag_type == FLAG_SPARSE)  // sparse
    flags_ptr = &(gen->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & VARGEN_VAR_P) && (mask & VARGEN_VAR_P)) { // active power
    if (flag_type == FLAG_VARS)
      gen->index_P = index;
    (*flags_ptr) |= VARGEN_VAR_P;
    index++;
  }
  if (!((*flags_ptr) & VARGEN_VAR_Q) && (mask & VARGEN_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS)
      gen->index_Q = index;
    (*flags_ptr) |= VARGEN_VAR_Q;
    index++;
  }
  return index;  
}

void VARGEN_set_var_values(Vargen* gen, Vec* values) {
  
  if (!gen)
    return;
  if (gen->vars & VARGEN_VAR_P) // active power (p.u.)
    gen->P = VEC_get(values,gen->index_P);
  if (gen->vars & VARGEN_VAR_Q) // reactive power (p.u.)
    gen->Q = VEC_get(values,gen->index_Q);
}

void VARGEN_show(Vargen* gen) {
  printf("vargen %d\t%d\n",
	 BUS_get_number(VARGEN_get_bus(gen)),
	 VARGEN_is_wind(gen));
}

