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
#include <pfnet/array.h>

struct Gen {

  // Bus
  Bus* bus;            /**< @brief Bus to which generator is connected */
  Bus* reg_bus;        /**< @brief Bus regulated by this generator */

  // Times
  int num_periods;   /**< @brief Number of time periods. */
  
  // Flags
  BOOL outage;         /**< @brief Flag for indicating that generator in on outage */
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Active power
  REAL* P;              /**< @brief Generator active power production (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum generator active power production (p.u.) */
  REAL P_min;          /**< @brief Minimum generator active power production (p.u.) */

  // Reactive power
  REAL* Q;              /**< @brief Generator reactive power production (p.u. system base power) */
  REAL Q_max;          /**< @brief Maximum generator reactive power production (p.u.) */
  REAL Q_min;          /**< @brief Minimum generator reactive power production (p.u.) */

  // Cost
  REAL cost_coeff_Q0;  /**< @brief Generator cost coefficient (constant term, units of $/hr ) */
  REAL cost_coeff_Q1;  /**< @brief Generator cost coefficient (linear term, units of $/(hr p.u.) ) */
  REAL cost_coeff_Q2;  /**< @brief Generator cost coefficient (quadratic term, units of $/(hr p.u.^2) ) */

  // Indices
  int index;           /**< @brief Generator index */
  int* index_P;         /**< @brief Active power index */
  int* index_Q;         /**< @brief Reactive power index */

  // Sensitivities
  REAL* sens_P_u_bound;  /**< @brief Sensitivity of active power upper bound */
  REAL* sens_P_l_bound;  /**< @brief Sensitivity of active power lower bound */

  // List
  Gen* next;     /**< @brief List of generators connected to a bus */
  Gen* reg_next; /**< @brief List of generators regulating a bus */
};

void* GEN_array_get(void* gen_array, int index) {
  if (gen_array)
    return (void*)&(((Gen*)gen_array)[index]);
  else 
    return NULL;
}

void GEN_array_del(Gen* gen_array, int size) {
  int i;
  Gen* gen;
  if (gen_array) {
    for (i = 0; i < size; i++) {
      gen = &(gen_array[i]);
      free(gen->P);
      free(gen->Q);
      free(gen->index_P);
      free(gen->index_Q);
      free(gen->sens_P_u_bound);
      free(gen->sens_P_l_bound);
    }
    free(gen_array);
  }  
}

Gen* GEN_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    Gen* gen_array = (Gen*)malloc(sizeof(Gen)*size);
    for (i = 0; i < size; i++) {
      GEN_init(&(gen_array[i]),num_periods);
      GEN_set_index(&(gen_array[i]),i);
    }
    return gen_array;
  }
  else
    return NULL;
}

void GEN_array_show(Gen* gen_array, int size, int t) {
  int i;
  if (gen_array) {
    for (i = 0; i < size; i++) 
      GEN_show(&(gen_array[i]),t);
  }
}

void GEN_clear_sensitivities(Gen* gen) {
  int t;
  if (gen) {
    for (t = 0; t < gen->num_periods; t++) {
      gen->sens_P_u_bound[t] = 0;
      gen->sens_P_l_bound[t] = 0;
    }
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

int GEN_get_num_periods(Gen* gen) {
  if (gen)
    return gen->num_periods;
  else
    return 0;
}

REAL GEN_get_sens_P_u_bound(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->sens_P_u_bound[t];
  else
    return 0;
}

REAL GEN_get_sens_P_l_bound(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->sens_P_l_bound[t];
  else
    return 0;
}

char GEN_get_obj_type(void* gen) {
  if (gen)
    return OBJ_GEN;
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

REAL GEN_get_P_cost(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return GEN_get_P_cost_for(gen,gen->P[t]); // $/hr
  else
    return 0;
}

REAL GEN_get_P_cost_for(Gen* gen, REAL P) {
  if (gen)
    return (gen->cost_coeff_Q0 + 
	    gen->cost_coeff_Q1*P +
	    gen->cost_coeff_Q2*pow(P,2.)); // $/hr
  else
    return 0;
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

int GEN_get_index_P(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_P[t];
  else
    return 0;
}

int GEN_get_index_Q(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_Q[t];
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

REAL GEN_get_P(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->P[t];
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
  if (gen)
    return gen->P_min;
  else 
    return 0;
}

REAL GEN_get_Q(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->Q[t];
  else
    return 0;
}

REAL GEN_get_Q_max(Gen* gen) {
  if (gen)
    return gen->Q_max;
  else
    return 0;
}

REAL GEN_get_Q_min(Gen* gen) {
  if (gen)
    return gen->Q_min;
  else
    return 0;
}

void GEN_get_var_values(Gen* gen, Vec* values, int code) {

  // Local vars
  int t;
  
  // No gen
  if (!gen)
    return;

  // Time loop
  for (t = 0; t < gen->num_periods; t++) {
    
    if (gen->vars & GEN_VAR_P) { // active power
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,gen->index_P[t],gen->P_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,gen->index_P[t],gen->P_min);
	break;
      default:
	VEC_set(values,gen->index_P[t],gen->P[t]);
      }
    }
    if (gen->vars & GEN_VAR_Q) { // reactive power
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,gen->index_Q[t],gen->Q_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,gen->index_Q[t],gen->Q_min);
	break;
      default:
	VEC_set(values,gen->index_Q[t],gen->Q[t]);
      }
    }
  }
}

Vec* GEN_get_var_indices(void* vgen, char var) {
  Gen* gen = (Gen*)vgen;
  Vec* indices;
  int t;
  if (!gen)
    return NULL;
  if (var == GEN_VAR_P) {
    indices = VEC_new(gen->num_periods);
    for (t = 0; t < gen->num_periods; t++)
      VEC_set(indices,t,gen->index_P[t]);
    return indices;
  }
  if (var == GEN_VAR_Q) {
    indices = VEC_new(gen->num_periods);
    for (t = 0; t < gen->num_periods; t++)
      VEC_set(indices,t,gen->index_Q[t]);
    return indices;
  }
  return NULL;
}

BOOL GEN_has_flags(void* vgen, char flag_type, char mask) {
  Gen* gen = (Gen*)vgen;
  if (gen) {
    if (flag_type == FLAG_VARS)
      return (gen->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (gen->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (gen->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (gen->sparse & mask) == mask;
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
  if ((prop & GEN_PROP_NOT_OUT) && GEN_is_on_outage(gen))
    return FALSE;
  if ((prop & GEN_PROP_P_ADJUST) && !GEN_is_P_adjustable(gen))
    return FALSE;
  return TRUE;
}

void GEN_init(Gen* gen, int num_periods) {

  // Local vars
  int T;

  // No gen
  if (!gen)
    return;

  T = num_periods;
  gen->num_periods = num_periods;
        
  gen->bus = NULL;
  gen->reg_bus = NULL;
  
  gen->outage = FALSE;
  gen->fixed = 0x00;
  gen->bounded = 0x00;
  gen->sparse = 0x00;
  gen->vars = 0x00;
  
  gen->P_max = 0;
  gen->P_min = 0;
    
  gen->Q_max = 0;
  gen->Q_min = 0;
  
  gen->cost_coeff_Q0 = 0;
  gen->cost_coeff_Q1 = 2000.;
  gen->cost_coeff_Q2 = 100.;
  
  gen->index = 0;

  ARRAY_zalloc(gen->P,REAL,T);
  ARRAY_zalloc(gen->Q,REAL,T);
  ARRAY_zalloc(gen->index_P,int,T);
  ARRAY_zalloc(gen->index_Q,int,T);
  ARRAY_zalloc(gen->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(gen->sens_P_l_bound,REAL,T);
  
  gen->next = NULL;
  gen->reg_next = NULL;
}

BOOL GEN_is_equal(Gen* gen, Gen* other) {
  return gen == other;
}

BOOL GEN_is_on_outage(Gen* gen) {
  if (gen)
    return gen->outage;
  else
    return FALSE;
}

BOOL GEN_is_P_adjustable(Gen* gen) {
  if (gen)
    return gen->P_min < gen->P_max;
  else
    return FALSE;
}

BOOL GEN_is_regulator(Gen* gen) {
  if (gen)
    return gen->reg_bus != NULL;
  else
    return FALSE;
}

BOOL GEN_is_slack(Gen* gen) {
  if (gen)
    return BUS_is_slack(gen->bus);
  else
    return FALSE;
}

Gen* GEN_list_add(Gen* gen_list, Gen* gen) {
  LIST_add(Gen,gen_list,gen,next);
  return gen_list;
}

Gen* GEN_list_del(Gen* gen_list, Gen* gen) {
  LIST_del(Gen,gen_list,gen,next);
  return gen_list;
}

int GEN_list_len(Gen* gen_list) {
  int len;
  LIST_len(Gen,gen_list,next,len);
  return len;
}

Gen* GEN_list_reg_add(Gen* reg_gen_list, Gen* reg_gen) {
  LIST_add(Gen,reg_gen_list,reg_gen,reg_next);
  return reg_gen_list;
}

Gen* GEN_list_reg_del(Gen* reg_gen_list, Gen* reg_gen) {
  LIST_del(Gen,reg_gen_list,reg_gen,reg_next);
  return reg_gen_list;
}

int GEN_list_reg_len(Gen* reg_gen_list) {
  int len;
  LIST_len(Gen,reg_gen_list,reg_next,len);
  return len;
}

Gen* GEN_new(int num_periods) {
  if (num_periods > 0) {
    Gen* gen = (Gen*)malloc(sizeof(Gen));
    GEN_init(gen,num_periods);
    return gen;
  }
  else
    return NULL;
}

void GEN_set_sens_P_u_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_P_u_bound[t] = value;
}

void GEN_set_sens_P_l_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_P_l_bound[t] = value;
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
  if (gen)
    gen->bus = bus;
}

void GEN_set_reg_bus(Gen* gen, Bus* reg_bus) {
  if (gen)
    gen->reg_bus = reg_bus;
}

void GEN_set_outage(Gen* gen, BOOL outage) {
  if (gen)
    gen->outage = outage;
}

void GEN_set_index(Gen* gen, int index) {
  if (gen)
    gen->index = index;
}

void GEN_set_P(Gen* gen, REAL P, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->P[t] = P;
}

void GEN_set_P_max(Gen* gen, REAL P_max) {
  if (gen)
    gen->P_max = P_max;
}

void GEN_set_P_min(Gen* gen, REAL P_min) {
  if (gen)
    gen->P_min = P_min;
}

void GEN_set_Q(Gen* gen, REAL Q, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->Q[t] = Q;
}

void GEN_set_Q_max(Gen* gen, REAL Q_max) {
  if (gen)
    gen->Q_max = Q_max;
}

void GEN_set_Q_min(Gen* gen, REAL Q_min) {
  if (gen)  
    gen->Q_min = Q_min;
}

int GEN_set_flags(void* vgen, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Gen* gen = (Gen*)vgen;
  int t;

  // No gen
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
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < gen->num_periods; t++)
	gen->index_P[t] = index+t;
    }
    (*flags_ptr) |= GEN_VAR_P;
    index += gen->num_periods;
  }
  if (!((*flags_ptr) & GEN_VAR_Q) && (mask & GEN_VAR_Q)) {
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < gen->num_periods; t++)
	gen->index_Q[t] = index+t;
    }
    (*flags_ptr) |= GEN_VAR_Q;
    index += gen->num_periods;
  }
  return index;  
}

void GEN_set_var_values(Gen* gen, Vec* values) {
 
  // Local vars
  int t;
 
  // No gen
  if (!gen)
    return;

  // Time loop
  for (t = 0; t < gen->num_periods; t++) {

    if (gen->vars & GEN_VAR_P) // active power (p.u.)
      gen->P[t] = VEC_get(values,gen->index_P[t]);
    if (gen->vars & GEN_VAR_Q) // reactive power (p.u.)
      gen->Q[t] = VEC_get(values,gen->index_Q[t]);
  }
}

void GEN_show(Gen* gen, int t) {
  printf("gen %d\t%d\n",
	 BUS_get_number(gen->bus),
	 BUS_get_number(gen->reg_bus));
}

void GEN_propagate_data_in_time(Gen* gen) {
  int t;
  if (gen) {
    for (t = 1; t < gen->num_periods; t++) {
      gen->P[t] = gen->P[0];
      gen->Q[t] = gen->Q[0];
    }
  }
}

