/** @file gen.c
 *  @brief This file defines the Gen data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/gen.h>
#include <pfnet/bus.h>
#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct Gen {

  // Bus
  Bus* bus;            /**< @brief Bus to which generator is connected. */
  Bus* reg_bus;        /**< @brief Bus regulated by this generator. */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Properties
  char name[GEN_BUFFER_SIZE]; /**< @brief Generator name */
  
  // Flags
  BOOL in_service;     /**< @brief Flag for indicating generator is in service. */
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value. */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded. */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables. */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse. */
  
  // Active power
  REAL* P;             /**< @brief Generator active power production (p.u. system base power). */
  REAL P_max;          /**< @brief Maximum generator active power production (p.u.). */
  REAL P_min;          /**< @brief Minimum generator active power production (p.u.). */
  REAL dP_max;         /**< @brief Maximum generator active power ramping (p.u.). */
  REAL P_prev;         /**< @brief Generator active power production during the previous time period (p.u.). */

  // Reactive power
  REAL* Q;             /**< @brief Generator reactive power production (p.u. system base power). */
  REAL Q_max;          /**< @brief Maximum generator reactive power production (p.u.). */
  REAL Q_min;          /**< @brief Minimum generator reactive power production (p.u.). */
  REAL Q_par;          /**< @brief Generator reactive power participation factor (unitless). */

  // Cost
  REAL cost_coeff_Q0;  /**< @brief Generator cost coefficient (constant term, units of $/hr ). */
  REAL cost_coeff_Q1;  /**< @brief Generator cost coefficient (linear term, units of $/(hr p.u.) ). */
  REAL cost_coeff_Q2;  /**< @brief Generator cost coefficient (quadratic term, units of $/(hr p.u.^2) ). */

  // Indices
  int index;           /**< @brief Generator index. */
  int* index_P;        /**< @brief Active power index. */
  int* index_Q;        /**< @brief Reactive power index. */

  // Sensitivities
  REAL* sens_P_u_bound;  /**< @brief Sensitivity of active power upper bound. */
  REAL* sens_P_l_bound;  /**< @brief Sensitivity of active power lower bound. */
  REAL* sens_Q_u_bound;  /**< @brief Sensitivity of reactive power upper bound. */
  REAL* sens_Q_l_bound;  /**< @brief Sensitivity of reactive power lower bound. */
  
  // Network
  Net* net; /**< @brief Network. */
  
  // List
  Gen* next;     /**< @brief List of generators connected to a bus. */
  Gen* reg_next; /**< @brief List of generators regulating a bus. */
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
      free(gen->sens_Q_u_bound);
      free(gen->sens_Q_l_bound);
      GEN_set_bus(gen,NULL);
      GEN_set_reg_bus(gen,NULL);
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
      snprintf(gen_array[i].name,(size_t)(GEN_BUFFER_SIZE-1),"%d",i);
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
      gen->sens_Q_u_bound[t] = 0;
      gen->sens_Q_l_bound[t] = 0;
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

void GEN_copy_from_gen(Gen* gen, Gen* other) {

  // Local variables
  int num_periods;

  // Check
  if (!gen || !other)
    return;

  // Min num periods
  if (gen->num_periods < other->num_periods)
    num_periods = gen->num_periods;
  else
    num_periods = other->num_periods;

  // Bus
  // skip buses

  // Times
  // skip num periods

  // Properties
  strcpy(gen->name,other->name);

  // Flags
  GEN_set_in_service(gen,GEN_is_in_service(other));
  gen->fixed = other->fixed;
  gen->bounded = other->bounded;
  gen->sparse = other->sparse;
  gen->vars = other->vars;

  // Active power
  memcpy(gen->P,other->P,num_periods*sizeof(REAL));
  gen->P_max = other->P_max;
  gen->P_min = other->P_min;
  gen->P_prev = other->P_prev;
  gen->dP_max = other->dP_max;

  // Reactive power
  memcpy(gen->Q,other->Q,num_periods*sizeof(REAL));
  gen->Q_max = other->Q_max;
  gen->Q_min = other->Q_min;
  gen->Q_par = other->Q_par;

  // Cost coefficients
  gen->cost_coeff_Q0 = other->cost_coeff_Q0;
  gen->cost_coeff_Q1 = other->cost_coeff_Q1;
  gen->cost_coeff_Q2 = other->cost_coeff_Q2;

  // Indices
  // skip index
  memcpy(gen->index_P,other->index_P,num_periods*sizeof(int));
  memcpy(gen->index_Q,other->index_Q,num_periods*sizeof(int));

  // Sensitivities
  memcpy(gen->sens_P_u_bound,other->sens_P_u_bound,num_periods*sizeof(REAL));
  memcpy(gen->sens_P_l_bound,other->sens_P_l_bound,num_periods*sizeof(REAL));
  memcpy(gen->sens_Q_u_bound,other->sens_Q_u_bound,num_periods*sizeof(REAL));
  memcpy(gen->sens_Q_l_bound,other->sens_Q_l_bound,num_periods*sizeof(REAL));

  // List
  // skip next 
}

char GEN_get_flags_vars(Gen* gen) {
  if (gen)
    return gen->vars;
  else
    return 0;
}

char GEN_get_flags_fixed(Gen* gen) {
  if (gen)
    return gen->fixed;
  else
    return 0;
}

char GEN_get_flags_bounded(Gen* gen) {
  if (gen)
    return gen->bounded;
  else
    return 0;
}

char GEN_get_flags_sparse(Gen* gen) {
  if (gen)
    return gen->sparse;
  else
    return 0;
}

char* GEN_get_name(Gen* gen) {
  if (gen)
    return gen->name;
  else
    return NULL;
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

REAL* GEN_get_sens_P_u_bound_array(Gen* gen) {
  if (gen)
    return gen->sens_P_u_bound;
  else
    return NULL;
}

REAL GEN_get_sens_P_l_bound(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->sens_P_l_bound[t];
  else
    return 0;
}

REAL* GEN_get_sens_P_l_bound_array(Gen* gen) {
  if (gen)
    return gen->sens_P_l_bound;
  else
    return NULL;
}

REAL GEN_get_sens_Q_u_bound(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->sens_Q_u_bound[t];
  else
    return 0;
}

REAL* GEN_get_sens_Q_u_bound_array(Gen* gen) {
  if (gen)
    return gen->sens_Q_u_bound;
  else
    return NULL;
}

REAL GEN_get_sens_Q_l_bound(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->sens_Q_l_bound[t];
  else
    return 0;
}

REAL* GEN_get_sens_Q_l_bound_array(Gen* gen) {
  if (gen)
    return gen->sens_Q_l_bound;
  else
    return NULL;
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
    return -1;
}

int GEN_get_index_P(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_P[t];
  else
    return -1;
}

int GEN_get_index_Q(Gen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_Q[t];
  else
    return -1;
}

int* GEN_get_index_P_array(Gen* gen) {
  if (gen)
    return gen->index_P;
  else
    return NULL;
}

int* GEN_get_index_Q_array(Gen* gen) {
  if (gen)
    return gen->index_Q;
  else
    return NULL;
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

REAL GEN_get_dP_max(Gen* gen) {
  if (gen)
    return gen->dP_max;
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

REAL GEN_get_P_prev(Gen* gen) {
  if (gen)
    return gen->P_prev;
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

REAL GEN_get_Q_par(Gen* gen) {
  if (gen)
    return gen->Q_par;
  else
    return 0;
}

REAL* GEN_get_P_array(Gen* gen) {
  if (gen)
    return gen->P;
  else
    return NULL;
}

REAL* GEN_get_Q_array(Gen* gen) {
  if (gen)
    return gen->Q;
  else
    return NULL;
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
	if (gen->bounded & GEN_VAR_P)
	  VEC_set(values,gen->index_P[t],gen->P_max);
	else
	  VEC_set(values,gen->index_P[t],GEN_INF_P);
	break;

      case LOWER_LIMITS:
	if (gen->bounded & GEN_VAR_P)
	  VEC_set(values,gen->index_P[t],gen->P_min);
	else
	  VEC_set(values,gen->index_P[t],-GEN_INF_P);
	break;

      default:
	VEC_set(values,gen->index_P[t],gen->P[t]);
      }
    }
    if (gen->vars & GEN_VAR_Q) { // reactive power
      switch(code) {

      case UPPER_LIMITS:
	if (gen->bounded & GEN_VAR_Q)
	  VEC_set(values,gen->index_Q[t],gen->Q_max);
	else
	  VEC_set(values,gen->index_Q[t],GEN_INF_Q);
	break;

      case LOWER_LIMITS:
	if (gen->bounded & GEN_VAR_Q)
	  VEC_set(values,gen->index_Q[t],gen->Q_min);
	else
	  VEC_set(values,gen->index_Q[t],-GEN_INF_Q);
	break;

      default:
	VEC_set(values,gen->index_Q[t],gen->Q[t]);
      }
    }
  }
}

char* GEN_get_var_info_string(Gen* gen, int index) {

  // Local variables
  char* info;

  //Check
  if (!gen)
    return NULL;

  // Active power
  if ((gen->vars & GEN_VAR_P) &&
      index >= gen->index_P[0] &&
      index <= gen->index_P[gen->num_periods-1]) {
    info = (char*)malloc(GEN_BUFFER_SIZE*sizeof(char));
    snprintf(info,GEN_BUFFER_SIZE*sizeof(char),
	     "generator:%d:active power:%d",gen->index,index-gen->index_P[0]);
    return info;
  }

  // Reactive power
  if ((gen->vars & GEN_VAR_Q) &&
      index >= gen->index_Q[0] &&
      index <= gen->index_Q[gen->num_periods-1]) {
    info = (char*)malloc(GEN_BUFFER_SIZE*sizeof(char));
    snprintf(info,GEN_BUFFER_SIZE*sizeof(char),
	     "generator:%d:reactive power:%d",gen->index,index-gen->index_Q[0]);
    return info;
  }

  // Return
  return NULL;
}

int GEN_get_num_vars(void* vgen, unsigned char var, int t_start, int t_end) {

  // Local vars
  Gen* gen = (Gen*)vgen;
  int num_vars = 0;
  int dt;

  // Checks
  if (!gen)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > gen->num_periods-1)
    t_end = gen->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & GEN_VAR_P) && (gen->vars & GEN_VAR_P))
    num_vars += dt;
  if ((var & GEN_VAR_Q) && (gen->vars & GEN_VAR_Q))
    num_vars += dt;
  return num_vars;
}

Vec* GEN_get_var_indices(void* vgen, unsigned char var, int t_start, int t_end) {

  // Local vars
  Gen* gen = (Gen*)vgen;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!gen)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > gen->num_periods-1)
    t_end = gen->num_periods-1;

  // Indices
  indices = VEC_new(GEN_get_num_vars(vgen,var,t_start,t_end));
  if ((var & GEN_VAR_P) && (gen->vars & GEN_VAR_P)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,gen->index_P[t]);
      offset++;
    }
  }
  if ((var & GEN_VAR_Q) && (gen->vars & GEN_VAR_Q)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,gen->index_Q[t]);
      offset++;
    }
  }
  return indices;
}

char* GEN_get_json_string(Gen* gen, char* output) {

  // Local variables
  char temp[GEN_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No gen
  if (!gen)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*GEN_BUFFER_SIZE*GEN_NUM_JSON_FIELDS*gen->num_periods);
    resize = TRUE;
  }
  output_start = output;

  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",gen->index,FALSE);
  JSON_obj(temp,output,"bus",gen->bus,BUS_get_index,FALSE);
  JSON_obj(temp,output,"reg_bus",gen->reg_bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"num_periods",gen->num_periods,FALSE);
  JSON_str(temp,output,"name",gen->name,FALSE);
  JSON_bool(temp,output,"in_service",GEN_is_in_service(gen),FALSE);
  JSON_array_float(temp,output,"P",gen->P,gen->num_periods,FALSE);
  JSON_float(temp,output,"P_max",gen->P_max,FALSE);
  JSON_float(temp,output,"P_min",gen->P_min,FALSE);
  JSON_float(temp,output,"dP_max",gen->dP_max,FALSE);
  JSON_float(temp,output,"P_prev",gen->P_prev,FALSE);
  JSON_array_float(temp,output,"Q",gen->Q,gen->num_periods,FALSE);
  JSON_float(temp,output,"Q_max",gen->Q_max,FALSE);
  JSON_float(temp,output,"Q_min",gen->Q_min,FALSE);
  JSON_float(temp,output,"Q_par",gen->Q_par,FALSE);
  JSON_float(temp,output,"cost_coeff_Q0",gen->cost_coeff_Q0,FALSE);
  JSON_float(temp,output,"cost_coeff_Q1",gen->cost_coeff_Q1,FALSE);
  JSON_float(temp,output,"cost_coeff_Q2",gen->cost_coeff_Q2,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL GEN_has_flags(void* vgen, char flag_type, unsigned char mask) {
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

  ARRAY_clear(gen->name,char,GEN_BUFFER_SIZE);
  
  gen->in_service = TRUE;
  gen->fixed = 0x00;
  gen->bounded = 0x00;
  gen->sparse = 0x00;
  gen->vars = 0x00;
  
  gen->dP_max = 0;
  gen->P_max = 0;
  gen->P_min = 0;
  gen->P_prev = 0;
    
  gen->Q_max = 0;
  gen->Q_min = 0;
  gen->Q_par = 1.;
  
  gen->cost_coeff_Q0 = 0;
  gen->cost_coeff_Q1 = 2000.;
  gen->cost_coeff_Q2 = 100.;
  
  gen->index = -1;

  ARRAY_zalloc(gen->P,REAL,T);
  ARRAY_zalloc(gen->Q,REAL,T);
  ARRAY_zalloc(gen->index_P,int,T);
  ARRAY_zalloc(gen->index_Q,int,T);
  ARRAY_zalloc(gen->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(gen->sens_P_l_bound,REAL,T);
  ARRAY_zalloc(gen->sens_Q_u_bound,REAL,T);
  ARRAY_zalloc(gen->sens_Q_l_bound,REAL,T);

  gen->net = NULL;
  
  gen->next = NULL;
  gen->reg_next = NULL;
}

BOOL GEN_is_equal(Gen* gen, Gen* other) {
  return gen == other;
}

BOOL GEN_is_in_service(void* gen) {
  if (gen)
    return ((Gen*)gen)->in_service && BUS_is_in_service(((Gen*)gen)->bus);
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

void GEN_set_name(Gen* gen, char* name) {
  if (gen)
    strncpy(gen->name,name,(size_t)(GEN_BUFFER_SIZE-1));
}

void GEN_set_sens_P_u_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_P_u_bound[t] = value;
}

void GEN_set_sens_P_l_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_P_l_bound[t] = value;
}

void GEN_set_sens_Q_u_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_Q_u_bound[t] = value;
}

void GEN_set_sens_Q_l_bound(Gen* gen, REAL value, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->sens_Q_l_bound[t] = value;
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
  Bus* old_bus;
  if (gen) {
    old_bus = gen->bus;        // save old bus
    gen->bus = NULL;           // disconnect old bus from gen
    BUS_del_gen(old_bus,gen);  // disconnect gen from old bus
    gen->bus = bus;            // connect new bus to gen
    BUS_add_gen(bus,gen);      // connect gen to new bus
  }
}

void GEN_set_reg_bus(Gen* gen, Bus* reg_bus) {
  Bus* old_reg_bus;
  if (gen) {
    old_reg_bus = gen->reg_bus;
    gen->reg_bus = NULL;
    BUS_del_reg_gen(old_reg_bus,gen);
    gen->reg_bus = reg_bus;
    BUS_add_reg_gen(gen->reg_bus,gen);
  }
}

void GEN_set_in_service(Gen* gen, BOOL in_service) {
  if (gen) {
    if (gen->in_service != in_service)
      NET_inc_state_tag(gen->net);
    gen->in_service = in_service;
  }
}

void GEN_set_index(Gen* gen, int index) {
  if (gen)
    gen->index = index;
}

void GEN_set_P(Gen* gen, REAL P, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->P[t] = P;
}

void GEN_set_dP_max(Gen* gen, REAL P) {
  if (gen)
    gen->dP_max = P;
}

void GEN_set_P_max(Gen* gen, REAL P) {
  if (gen)
    gen->P_max = P;
}

void GEN_set_P_min(Gen* gen, REAL P) {
  if (gen)
    gen->P_min = P;
}

void GEN_set_P_prev(Gen* gen, REAL P) {
  if (gen)
    gen->P_prev = P;
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

void GEN_set_Q_par(Gen* gen, REAL Q_par) {
  if (gen)  
    gen->Q_par = Q_par;
}

void GEN_set_network(Gen* gen, void* net) {
  if (gen)
    gen->net = (Net*)net;
}

int GEN_set_flags(void* vgen, char flag_type, unsigned char mask, int index) {

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

void GEN_propagate_data_in_time(Gen* gen, int start, int end) {
  int t;
  if (gen) {
    if (start < 0)
      start = 0;
    if (end > gen->num_periods)
      end = gen->num_periods;
    for (t = start+1; t < end; t++) {
      gen->P[t] = gen->P[start];
      gen->Q[t] = gen->Q[start];
    }
  }
}

