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
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct Vargen {

  // Bus
  Bus* bus;            /**< @brief Bus to which variable generator is connected */

  // Times
  int num_periods;     /**< @brief Number of time periods. */

  // Properties
  char name[VARGEN_BUFFER_SIZE]; /**< @brief Variable generator name */
  char type;                     /**< @brief Variable generator type */
  
  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Active power
  REAL* P;             /**< @brief Variable generator active power after curtailments (p.u. system base power) */
  REAL* P_ava;         /**< @brief Variable generator available active power (p.u. system base power) */
  REAL P_max;          /**< @brief Maximum variable generator active power (p.u.) */
  REAL P_min;          /**< @brief Minimum variable generator active power (p.u.) */
  REAL* P_std;         /**< @brief Standard deviation of active power (p.u. system base power) */
  
  // Reactive power
  REAL* Q;             /**< @brief Variable generator reactive power (p.u. system base power) */
  REAL Q_max;          /**< @brief Maximum variable generator reactive power (p.u.) */
  REAL Q_min;          /**< @brief Minimum variable generator reactive power (p.u.) */ 
  
  // Indices
  int index;           /**< @brief Generator index */
  int* index_P;        /**< @brief Active power index */
  int* index_Q;        /**< @brief Reactive power index */

  // List
  Vargen* next;        /**< @brief List of variable generators connected to a bus */
};

void* VARGEN_array_get(void* gen_array, int index) {
  if (gen_array)
    return (void*)&(((Vargen*)gen_array)[index]);
  else 
    return NULL;
}

void VARGEN_array_del(Vargen* gen_array, int size) {
  int i;
  Vargen* gen;
  if (gen_array) {    
    for (i = 0; i < size; i++) {
      gen = &(gen_array[i]);
      free(gen->P);
      free(gen->P_ava);
      free(gen->P_std);
      free(gen->Q);
      free(gen->index_P);
      free(gen->index_Q);
      VARGEN_set_bus(gen,NULL);
    }
    free(gen_array);
  }  
}

Vargen* VARGEN_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    Vargen* gen_array = (Vargen*)malloc(sizeof(Vargen)*size);
    for (i = 0; i < size; i++) {
      VARGEN_init(&(gen_array[i]),num_periods);
      VARGEN_set_index(&(gen_array[i]),i);
      snprintf(gen_array[i].name,(size_t)(VARGEN_BUFFER_SIZE-1),"%d",i);
    }
    return gen_array;
  }
  else
    return NULL;
}

void VARGEN_array_show(Vargen* gen_array, int size, int t) {
  int i;
  if (gen_array) {
    for (i = 0; i < size; i++) 
      VARGEN_show(&(gen_array[i]),t);
  }
}

void VARGEN_clear_sensitivities(Vargen* vargen) {

  
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

void VARGEN_copy_from_vargen(Vargen* gen, Vargen* other) {
  
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
  gen->type = other->type;

  // Flags
  gen->fixed = other->fixed;
  gen->bounded = other->bounded;
  gen->sparse = other->sparse;
  gen->vars = other->vars;

  // Active power
  memcpy(gen->P,other->P,num_periods*sizeof(REAL));
  memcpy(gen->P_ava,other->P_ava,num_periods*sizeof(REAL));
  memcpy(gen->P_std,other->P_std,num_periods*sizeof(REAL));
  gen->P_max = other->P_max;
  gen->P_min = other->P_min;

  // Reactive power
  memcpy(gen->Q,other->Q,num_periods*sizeof(REAL));
  gen->Q_max = other->Q_max;
  gen->Q_min = other->Q_min;

  // Indices
  // skip index
  memcpy(gen->index_P,other->index_P,num_periods*sizeof(int));
  memcpy(gen->index_Q,other->index_Q,num_periods*sizeof(int));

  // List
  // skip next    
}

char VARGEN_get_flags_vars(Vargen* gen) {
  if (gen)
    return gen->vars;
  else
    return 0;
}

char VARGEN_get_flags_fixed(Vargen* gen) {
  if (gen)
    return gen->fixed;
  else
    return 0;
}

char VARGEN_get_flags_bounded(Vargen* gen) {
  if (gen)
    return gen->bounded;
  else
    return 0;
}

char VARGEN_get_flags_sparse(Vargen* gen) {
  if (gen)
    return gen->sparse;
  else
    return 0;
}

int VARGEN_get_num_periods(Vargen* gen) {
  if (gen)
    return gen->num_periods;
  else
    return 0;
}

char* VARGEN_get_name(Vargen* gen) {
  if (gen)
    return gen->name;
  else
    return NULL;
}

char VARGEN_get_obj_type(void* gen) {
  if (gen)
    return OBJ_VARGEN;
  else
    return OBJ_UNKNOWN;
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
    return -1;
}

int VARGEN_get_index_P(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_P[t];
  else
    return -1;
}

int VARGEN_get_index_Q(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_Q[t];
  else
    return -1;
}

Vargen* VARGEN_get_next(Vargen* gen) {
  if (gen)
    return gen->next;
  else
    return NULL;
}

REAL VARGEN_get_P(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->P[t];
  else 
    return 0;
}

REAL VARGEN_get_P_ava(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->P_ava[t];
  else 
    return 0;
}

REAL VARGEN_get_P_max(Vargen* gen) {
  if (gen)
    return gen->P_max;
  else 
    return 0;
}

REAL VARGEN_get_P_min(Vargen* gen) {
  if (gen)
    return gen->P_min;
  else 
    return 0;
}

REAL VARGEN_get_P_std(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->P_std[t];
  else 
    return 0;
}

REAL VARGEN_get_Q(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->Q[t];
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

  // Local vars
  int t;

  // No vargen
  if (!gen)
    return;

  // Time loop
  for (t = 0; t < gen->num_periods; t++) {

    // active power
    if (gen->vars & VARGEN_VAR_P) {
      switch(code) {
	
      case UPPER_LIMITS:
	if (gen->bounded & VARGEN_VAR_P)
	  VEC_set(values,gen->index_P[t],gen->P_ava[t]);
	else
	  VEC_set(values,gen->index_P[t],VARGEN_INF_P);
	break;

      case LOWER_LIMITS:
	if (gen->bounded & VARGEN_VAR_P)
	  VEC_set(values,gen->index_P[t],gen->P_min);
	else
	  VEC_set(values,gen->index_P[t],-VARGEN_INF_P);
	break;

      default:
	VEC_set(values,gen->index_P[t],gen->P[t]);
      }
    }
    
    // reactive power
    if (gen->vars & VARGEN_VAR_Q) {
      switch(code) {

      case UPPER_LIMITS:
	if (gen->bounded & VARGEN_VAR_Q)
	  VEC_set(values,gen->index_Q[t],gen->Q_max);
	else
	  VEC_set(values,gen->index_Q[t],VARGEN_INF_Q);
	break;

      case LOWER_LIMITS:
	if (gen->bounded & VARGEN_VAR_Q)
	  VEC_set(values,gen->index_Q[t],gen->Q_min);
	else
	  VEC_set(values,gen->index_Q[t],-VARGEN_INF_Q);
	break;

      default:
	VEC_set(values,gen->index_Q[t],gen->Q[t]);
      }
    }
  }
}

char* VARGEN_get_var_info_string(Vargen* gen, int index) {

  // Local variables
  char* info;

  //Check
  if (!gen)
    return NULL;

  // Active power
  if ((gen->vars & VARGEN_VAR_P) &&
      index >= gen->index_P[0] &&
      index <= gen->index_P[gen->num_periods-1]) {
    info = (char*)malloc(VARGEN_BUFFER_SIZE*sizeof(char));
    snprintf(info,VARGEN_BUFFER_SIZE*sizeof(char),
	     "variable generator:%d:active power:%d",gen->index,index-gen->index_P[0]);
    return info;
  }
  
  // Reactive power
  if ((gen->vars & VARGEN_VAR_Q) &&
      index >= gen->index_Q[0] &&
      index <= gen->index_Q[gen->num_periods-1]) {
    info = (char*)malloc(VARGEN_BUFFER_SIZE*sizeof(char));
    snprintf(info,VARGEN_BUFFER_SIZE*sizeof(char),
	     "variable generator:%d:reactive power:%d",gen->index,index-gen->index_Q[0]);
    return info;
  }

  // Return
  return NULL;
}

int VARGEN_get_num_vars(void* vgen, unsigned char var, int t_start, int t_end) {

  // Local vars
  Vargen* gen = (Vargen*)vgen;
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
  if ((var & VARGEN_VAR_P) && (gen->vars & VARGEN_VAR_P))
    num_vars += dt;
  if ((var & VARGEN_VAR_Q) && (gen->vars & VARGEN_VAR_Q))
    num_vars += dt;
  return num_vars;
}

Vec* VARGEN_get_var_indices(void* vgen, unsigned char var, int t_start, int t_end) {

  // Local vars
  Vargen* gen = (Vargen*)vgen;
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
  indices = VEC_new(VARGEN_get_num_vars(vgen,var,t_start,t_end));
  if ((var & VARGEN_VAR_P) && (gen->vars & VARGEN_VAR_P)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,gen->index_P[t]);
      offset++;
    }
  }
  if ((var & VARGEN_VAR_Q) && (gen->vars & VARGEN_VAR_Q)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,gen->index_Q[t]);
      offset++;
    }
  }
  return indices;
}

char* VARGEN_get_json_string(Vargen* gen, char* output) {

  // Local variables
  char temp[VARGEN_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;  

  // No gen
  if (!gen)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*VARGEN_JSON_BUFFER_SIZE*VARGEN_NUM_JSON_FIELDS*gen->num_periods);
    resize = TRUE;
  }
  output_start = output;

  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",gen->index,FALSE);
  JSON_obj(temp,output,"bus",gen->bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"num_periods",gen->num_periods,FALSE);
  JSON_str(temp,output,"name",gen->name,FALSE);
  JSON_array_float(temp,output,"P",gen->P,gen->num_periods,FALSE);
  JSON_array_float(temp,output,"P_ava",gen->P_ava,gen->num_periods,FALSE)
  JSON_float(temp,output,"P_max",gen->P_max,FALSE);
  JSON_float(temp,output,"P_min",gen->P_min,FALSE);
  JSON_array_float(temp,output,"P_std",gen->P_std,gen->num_periods,FALSE);
  JSON_array_float(temp,output,"Q",gen->Q,gen->num_periods,FALSE);
  JSON_float(temp,output,"Q_max",gen->Q_max,FALSE);
  JSON_float(temp,output,"Q_min",gen->Q_min,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}


BOOL VARGEN_has_flags(void* vgen, char flag_type, unsigned char mask) {
  Vargen* gen = (Vargen*)vgen;
  if (gen) {
    if (flag_type == FLAG_VARS)         // variables
      return (gen->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED) // bounded
      return (gen->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)   // fixed
      return (gen->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)  // sparse
      return (gen->sparse & mask) == mask;
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

void VARGEN_init(Vargen* gen, int num_periods) {
  
  // Local variables
  int T;

  // No vargen
  if (!gen)
    return;

  T = num_periods;
  gen->num_periods = num_periods;

  gen->bus = NULL;
  gen->type = VARGEN_TYPE_WIND;
  ARRAY_clear(gen->name,char,VARGEN_BUFFER_SIZE);
  gen->fixed = 0x00;
  gen->bounded = 0x00;
  gen->sparse = 0x00;
  gen->vars = 0x00;
  gen->P_max = 0;
  gen->P_min = 0;
  gen->Q_max = 0;
  gen->Q_min = 0;
  
  gen->index = -1;

  gen->P = NULL;
  gen->P_ava = NULL;
  gen->P_std = NULL;
  gen->Q = NULL;
  gen->index_P = NULL;
  gen->index_Q = NULL;

  ARRAY_zalloc(gen->P,REAL,T);
  ARRAY_zalloc(gen->P_ava,REAL,T);
  ARRAY_zalloc(gen->P_std,REAL,T);
  ARRAY_zalloc(gen->Q,REAL,T);
  ARRAY_zalloc(gen->index_P,int,T);
  ARRAY_zalloc(gen->index_Q,int,T);
  
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

BOOL VARGEN_is_equal(Vargen* gen, Vargen* other) {
  return gen == other;
}

Vargen* VARGEN_list_add(Vargen* gen_list, Vargen* gen) {
  LIST_add(Vargen,gen_list,gen,next);
  return gen_list;
}

Vargen* VARGEN_list_del(Vargen* gen_list, Vargen* gen) {
  LIST_del(Vargen,gen_list,gen,next);
  return gen_list;
}

int VARGEN_list_len(Vargen* gen_list) {
  int len;
  LIST_len(Vargen,gen_list,next,len);
  return len;
}

Vargen* VARGEN_new(int num_periods) {
  if (num_periods > 0) {
    Vargen* gen = (Vargen*)malloc(sizeof(Vargen));
    VARGEN_init(gen,num_periods);
    return gen;
  }
  else
    return NULL;
}

void VARGEN_set_name(Vargen* gen, char* name) {
  if (gen)
    strncpy(gen->name,name,(size_t)(VARGEN_BUFFER_SIZE-1));
}

void VARGEN_set_bus(Vargen* gen, Bus* bus) {
  Bus* old_bus;
  if (gen) {
    old_bus = gen->bus;
    gen->bus = NULL;
    BUS_del_vargen(old_bus,gen);
    gen->bus = bus;
    BUS_add_vargen(gen->bus,gen);
  }
}

void VARGEN_set_index(Vargen* gen, int index) {
  if (gen)
    gen->index = index;
}

void VARGEN_set_P(Vargen* gen, REAL P, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->P[t] = P;
}

void VARGEN_set_P_ava(Vargen* gen, REAL P, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->P_ava[t] = P;
}

void VARGEN_set_P_max(Vargen* gen, REAL P_max) {
  if (gen)
    gen->P_max = P_max;
}

void VARGEN_set_P_min(Vargen* gen, REAL P_min) {
  if (gen)
    gen->P_min = P_min;
}

void VARGEN_set_P_std(Vargen* gen, REAL P_std, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->P_std[t] = P_std;
}

void VARGEN_set_Q(Vargen* gen, REAL Q, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    gen->Q[t] = Q;
}

void VARGEN_set_Q_max(Vargen* gen, REAL Q) {
  if (gen)
    gen->Q_max = Q;
}

void VARGEN_set_Q_min(Vargen* gen, REAL Q) {
  if (gen)
    gen->Q_min = Q;
}

int VARGEN_set_flags(void* vgen, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Vargen* gen = (Vargen*)vgen;
  int t;

  // No vargen
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
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < gen->num_periods; t++)
	gen->index_P[t] = index+t;
    }
    (*flags_ptr) |= VARGEN_VAR_P;
    index += gen->num_periods;
  }
  if (!((*flags_ptr) & VARGEN_VAR_Q) && (mask & VARGEN_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < gen->num_periods; t++)
	gen->index_Q[t] = index+t;
    }
    (*flags_ptr) |= VARGEN_VAR_Q;
    index += gen->num_periods;
  }
  return index;  
}

void VARGEN_set_var_values(Vargen* gen, Vec* values) {

  // Local vars
  int t;
  
  // No vargen
  if (!gen)
    return;

  // Time loop
  for (t = 0; t < gen->num_periods; t++) {
    
    if (gen->vars & VARGEN_VAR_P) // active power (p.u.)
      gen->P[t] = VEC_get(values,gen->index_P[t]);
    if (gen->vars & VARGEN_VAR_Q) // reactive power (p.u.)
      gen->Q[t] = VEC_get(values,gen->index_Q[t]);
  }
}

void VARGEN_show(Vargen* gen, int t) {
  printf("vargen %d\t%d\n",
	 BUS_get_number(VARGEN_get_bus(gen)),
	 VARGEN_is_wind(gen));
}

void VARGEN_propagate_data_in_time(Vargen* gen, int start, int end) {
  int t;
  if (gen) {
    if (start < 0)
      start = 0;
    if (end > gen->num_periods)
      end = gen->num_periods;
    for (t = start+1; t < end; t++) {
      gen->P[t] = gen->P[start];
      gen->P_ava[t] = gen->P_ava[start];
      gen->P_std[t] = gen->P_std[start];
      gen->Q[t] = gen->Q[start];
    }
  }
}

