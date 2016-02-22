/** @file gen.c
 *  @brief This file defines the Gen data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/gen.h>
#include <pfnet/bus.h>

struct Gen {

  // Properties
  char obj_type;       /**< @brief Object type */

  // Bus
  Bus* bus;            /**< @brief Bus to which generator is connected */
  Bus* reg_bus;        /**< @brief Bus regulated by this generator */
  
  // Flags
  BOOL regulator;      /**< @brief Flag for indicating that this generator provides voltage regulation */
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Active power
  REAL P;              /**< @brief Generator active power (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum generator active power (p.u.) */
  REAL P_min;          /**< @brief Minimum generator active power (p.u.) */

  // Reactive power
  REAL Q;              /**< @brief Generator reactive power (p.u. system base power) */
  REAL Q_max;          /**< @brief Maximum generator reactive power (p.u.) */
  REAL Q_min;          /**< @brief Minimum generator reactive power (p.u.) */

  // Cost
  REAL cost_coeff_Q0;  /**< @brief Generator cost coefficient (constant term) */
  REAL cost_coeff_Q1;  /**< @brief Generator cost coefficient (linear term) */
  REAL cost_coeff_Q2;  /**< @brief Generator cost coefficient (quadratic term) */

  // Indices
  int index;           /**< @brief Generator index */
  int index_P;         /**< @brief Active power index */
  int index_Q;         /**< @brief Reactive power index */

  // List
  Gen* next;     /**< @brief List of generators connected to a bus */
  Gen* reg_next; /**< @brief List of generators regulating a bus */
};

void* GEN_array_get(void* gen, int index) {
  if (gen)
    return (void*)&(((Gen*)gen)[index]);
  else 
    return NULL;
}

Gen* GEN_array_new(int num) {
  int i;
  Gen* gen = (Gen*)malloc(sizeof(Gen)*num);
  for (i = 0; i < num; i++) {
    GEN_init(&(gen[i]));
    GEN_set_index(&(gen[i]),i);
  }
  return gen;
}

void GEN_array_show(Gen* gen, int num) {
  int i;
  if (gen) {
    for (i = 0; i < num; i++) 
      GEN_show(&(gen[i]));
  }
}

void GEN_clear_flags(Gen* gen, char flag_type) {
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

char GEN_get_obj_type(Gen* gen) {
  if (gen)
    return gen->obj_type;
  else
    return OBJ_UNKNOWN;
}

Bus* GEN_get_bus(Gen* gen) {
  if (gen)
    return gen->bus;
  else
    return NULL;
}

Bus* GEN_get_reg_bus(Gen* gen) {
  if (gen)
    return gen->reg_bus;
  else
    return NULL;
}

REAL GEN_get_cost_coeff_Q0(Gen* gen) {
  if (gen)
    return gen->cost_coeff_Q0;
  else
    return 0;
}

REAL GEN_get_cost_coeff_Q1(Gen* gen) {
  if (gen)
    return gen->cost_coeff_Q1;
  else
    return 0;
}

REAL GEN_get_cost_coeff_Q2(Gen* gen) {
  if (gen)
    return gen->cost_coeff_Q2;
  else
    return 0;
}

int GEN_get_index(Gen* gen) {
  if (gen)
    return gen->index;
  else
    return 0;
}

int GEN_get_index_P(Gen* gen) {
  if (gen)
    return gen->index_P;
  else
    return 0;
}

int GEN_get_index_Q(Gen* gen) {
  if (gen)
    return gen->index_Q;
  else
    return 0;
}

Gen* GEN_get_next(Gen* gen) {
  if (gen)
    return gen->next;
  else
    return NULL;
}

Gen* GEN_get_reg_next(Gen* gen) {
  if (gen)
    return gen->reg_next;
  else
    return NULL;
}

REAL GEN_get_P(Gen* gen) {
  if (gen)
    return gen->P;
  else 
    return 0;
}

REAL GEN_get_P_max(Gen* gen) {
  if (gen)
    return gen->P_max;
  else 
    return 0;
}

REAL GEN_get_P_min(Gen* gen) {
  return gen->P_min;
}

REAL GEN_get_Q(Gen* gen) {
  return gen->Q;
}

REAL GEN_get_Q_max(Gen* gen) {
  return gen->Q_max;
}

REAL GEN_get_Q_min(Gen* gen) {
  return gen->Q_min;
}

void GEN_get_var_values(Gen* gen, Vec* values, int code) {
  
  if (!gen)
    return;

  if (gen->vars & GEN_VAR_P) { // active power
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,gen->index_P,gen->P_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,gen->index_P,gen->P_min);
      break;
    default:
      VEC_set(values,gen->index_P,gen->P);
    }
  }
  if (gen->vars & GEN_VAR_Q) { // reactive power
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

int GEN_get_var_index(void* vgen, char var) {
  Gen* gen = (Gen*)vgen;
  if (!gen)
    return 0;
  if (var == GEN_VAR_P)
    return gen->index_P;
  if (var == GEN_VAR_Q)
    return gen->index_Q;
  return 0;
}

BOOL GEN_has_flags(void* vgen, char flag_type, char mask) {
  Gen* gen = (Gen*)vgen;
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

BOOL GEN_has_properties(void* vgen, char prop) {
  Gen* gen = (Gen*)vgen;
  if (!gen)
    return FALSE;
  if ((prop & GEN_PROP_SLACK) && !GEN_is_slack(gen))
    return FALSE;
  if ((prop & GEN_PROP_REG) && !GEN_is_regulator(gen))
    return FALSE;    
  if ((prop & GEN_PROP_NOT_REG) && GEN_is_regulator(gen))
    return FALSE;
  if ((prop & GEN_PROP_NOT_SLACK) && GEN_is_slack(gen))
    return FALSE;
  if ((prop & GEN_PROP_P_ADJUST) && !GEN_is_P_adjustable(gen))
    return FALSE;
  return TRUE;
}

void GEN_init(Gen* gen) {
  gen->obj_type = OBJ_GEN;
  gen->bus = NULL;
  gen->reg_bus = NULL;
  gen->regulator = FALSE;
  gen->fixed = 0x00;
  gen->bounded = 0x00;
  gen->sparse = 0x00;
  gen->vars = 0x00;
  gen->P = 0;
  gen->P_max = 0;
  gen->P_min = 0;
  gen->Q = 0;
  gen->Q_max = 0;
  gen->Q_min = 0;
  gen->cost_coeff_Q0 = 0/100.;
  gen->cost_coeff_Q1 = 20.;
  gen->cost_coeff_Q2 = 0.01*100.;
  gen->index = 0;
  gen->index_P = 0;
  gen->index_Q = 0;
  gen->next = NULL;
  gen->reg_next = NULL;
}

BOOL GEN_is_P_adjustable(Gen* gen) {
  if (gen)
    return gen->P_min < gen->P_max;
  else
    return FALSE;
}

BOOL GEN_is_regulator(Gen* gen) {
  if (gen)
    return gen->regulator;
  else
    return FALSE;
}

BOOL GEN_is_slack(Gen* gen) {
  if (gen)
    return BUS_is_slack(gen->bus);
  else
    return FALSE;
}

Gen* GEN_list_add(Gen *gen_list, Gen* gen) {
  LIST_add(gen_list,gen,next);
  return gen_list;
}

int GEN_list_len(Gen* gen_list) {
  int len;
  LIST_len(Gen,gen_list,next,len);
  return len;
}

Gen* GEN_list_reg_add(Gen* reg_gen_list, Gen* reg_gen) {
  LIST_add(reg_gen_list,reg_gen,reg_next);
  return reg_gen_list;
}

int GEN_list_reg_len(Gen* reg_gen_list) {
  int len;
  LIST_len(Gen,reg_gen_list,reg_next,len);
  return len;
}

Gen* GEN_new(void) {
  Gen* gen = (Gen*)malloc(sizeof(Gen));
  GEN_init(gen);
  return gen;
}

void GEN_set_cost_coeff_Q0(Gen* gen, REAL q) {
  if (gen)
    gen->cost_coeff_Q0 = q;
}

void GEN_set_cost_coeff_Q1(Gen* gen, REAL q) {
  if (gen)
    gen->cost_coeff_Q1 = q;
}

void GEN_set_cost_coeff_Q2(Gen* gen, REAL q) {
  if (gen)
    gen->cost_coeff_Q2 = q;
}

void GEN_set_bus(Gen* gen, Bus* bus) {
  gen->bus = (Bus*)bus;
}

void GEN_set_reg_bus(Gen* gen, Bus* reg_bus) {
  gen->reg_bus = (Bus*)reg_bus;
}

void GEN_set_regulator(Gen* gen, BOOL regulator) {
  gen->regulator = regulator;
}

void GEN_set_index(Gen* gen, int index) {
  gen->index = index;
}

void GEN_set_P(Gen* gen, REAL P) {
  gen->P = P;
}

void GEN_set_P_max(Gen* gen, REAL P_max) {
  gen->P_max = P_max;
}

void GEN_set_P_min(Gen* gen, REAL P_min) {
  gen->P_min = P_min;
}

void GEN_set_Q(Gen* gen, REAL Q) {
  gen->Q = Q;
}

void GEN_set_Q_max(Gen* gen, REAL Q_max) {
  gen->Q_max = Q_max;
}

void GEN_set_Q_min(Gen* gen, REAL Q_min) {
  gen->Q_min = Q_min;
}

int GEN_set_flags(void* vgen, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Gen* gen = (Gen*)vgen;

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
  if (!((*flags_ptr) & GEN_VAR_P) && (mask & GEN_VAR_P)) {
    if (flag_type == FLAG_VARS)
      gen->index_P = index;
    (*flags_ptr) |= GEN_VAR_P;
    index++;
  }
  if (!((*flags_ptr) & GEN_VAR_Q) && (mask & GEN_VAR_Q)) {
    if (flag_type == FLAG_VARS)
      gen->index_Q = index;
    (*flags_ptr) |= GEN_VAR_Q;
    index++;
  }
  return index;  
}

void GEN_set_var_values(Gen* gen, Vec* values) {
  
  if (!gen)
    return;
  if (gen->vars & GEN_VAR_P) // active power (p.u.)
    gen->P = VEC_get(values,gen->index_P);
  if (gen->vars & GEN_VAR_Q) // reactive power (p.u.)
    gen->Q = VEC_get(values,gen->index_Q);
}

void GEN_show(Gen* gen) {
  printf("gen %d\t%d\n",BUS_get_number(gen->bus),BUS_get_number(gen->reg_bus));
}

