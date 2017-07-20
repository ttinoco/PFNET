/** @file vargen.c
 *  @brief This file defines the Vargen data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/vargen.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>

struct Vargen {

  // Bus
  Bus* bus;            /**< @brief Bus to which variable generator is connected */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

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

  // Hash
  UT_hash_handle hh;   /**< @brief Handle for vargen hash table based on names */

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
      snprintf(gen_array[i].name,(size_t)(VARGEN_BUFFER_SIZE-1),"VARGEN %d",i+1);
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
    return 0;
}

int VARGEN_get_index_P(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_P[t];
  else
    return 0;
}

int VARGEN_get_index_Q(Vargen* gen, int t) {
  if (gen && t >= 0 && t < gen->num_periods)
    return gen->index_Q[t];
  else
    return 0;
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

char* VARGEN_get_json_string(Vargen* gen) {

  // Local variables
  char* temp;
  char* output;
  int i;

  // No gen
  if (!gen)
    return NULL;

  // Alloc
  temp = (char*)malloc(sizeof(char)*VARGEN_BUFFER_SIZE);
  output = (char*)malloc(sizeof(char)*VARGEN_BUFFER_SIZE*15*gen->num_periods);

  // Start
  strcpy(output,"{ ");

  // Bus
  if (gen->bus)
    sprintf(temp,"\"bus\" : %d", BUS_get_index(gen->bus));
  else
    sprintf(temp,"\"bus\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // Num periods
  sprintf(temp,"\"num_periods\" : %d", gen->num_periods);
  strcat(output,temp);
  strcat(output,", ");

  // Name
  sprintf(temp,"\"name\" : \"%s\"", gen->name);
  strcat(output,temp);
  strcat(output,", ");

  // Type
  sprintf(temp,"\"type\" : %d", gen->type);
  strcat(output,temp);
  strcat(output,", ");
  
  // P
  strcat(output,"\"P\" : [ ");
  for (i = 0; i < gen->num_periods; i++) {
    sprintf(temp,"%.10e", gen->P[i]);
    strcat(output,temp);
    if (i < gen->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");

  // P ava
  strcat(output,"\"P_ava\" : [ ");
  for (i = 0; i < gen->num_periods; i++) {
    sprintf(temp,"%.10e", gen->P_ava[i]);
    strcat(output,temp);
    if (i < gen->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");

  // P max
  sprintf(temp,"\"P_max\" : %.10e", gen->P_max);
  strcat(output,temp);
  strcat(output,", ");

  // P min
  sprintf(temp,"\"P_min\" : %.10e", gen->P_min);
  strcat(output,temp);
  strcat(output,", ");

  // P std
  strcat(output,"\"P_std\" : [ ");
  for (i = 0; i < gen->num_periods; i++) {
    sprintf(temp,"%.10e", gen->P_std[i]);
    strcat(output,temp);
    if (i < gen->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");

  // Q
  strcat(output,"\"Q\" : [ ");
  for (i = 0; i < gen->num_periods; i++) {
    sprintf(temp,"%.10e", gen->Q[i]);
    strcat(output,temp);
    if (i < gen->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  
  // Q max
  sprintf(temp,"\"Q_max\" : %.10e", gen->Q_max);
  strcat(output,temp);
  strcat(output,", ");

  // Q min
  sprintf(temp,"\"Q_min\" : %.10e", gen->Q_min);
  strcat(output,temp);
  strcat(output,", ");

  // Index
  sprintf(temp,"\"index\" : %d", gen->index);
  strcat(output,temp);
  strcat(output,", ");

  // Next
  if (gen->next)
    sprintf(temp,"\"next\" : %d", VARGEN_get_index(gen->next));
  else
    sprintf(temp,"\"next\" : %s", "null");
  strcat(output,temp);
  strcat(output,"");

   // End
  strcat(output," }");
  
  // Free 
  free(temp);
  output = (char*)realloc(output,sizeof(char)*(strlen(output)+1)); // +1 important!
  
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

Vargen* VARGEN_hash_name_add(Vargen* vargen_hash, Vargen* vargen) {
  HASH_ADD_STR(vargen_hash,name,vargen);
  return vargen_hash;
}

void VARGEN_hash_name_del(Vargen* vargen_hash) {
  while (vargen_hash != NULL)
    HASH_DEL(vargen_hash,vargen_hash);
}

Vargen* VARGEN_hash_name_find(Vargen* vargen_hash, char* name) {
  Vargen* vargen;
  HASH_FIND_STR(vargen_hash,name,vargen);
  return vargen;
}

int VARGEN_hash_name_len(Vargen* vargen_hash) {
  return HASH_COUNT(vargen_hash);
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
  gen->index = 0;

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

Vargen* VARGEN_list_add(Vargen* gen_list, Vargen* gen) {
  LIST_add(Vargen,gen_list,gen,next);
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
  if (gen)
    gen->bus = (Bus*)bus;
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

void VARGEN_propagate_data_in_time(Vargen* gen) {
  int t;
  if (gen) {
    for (t = 1; t < gen->num_periods; t++) {
      gen->P[t] = gen->P[0];
      gen->P_ava[t] = gen->P_ava[0];
      gen->P_std[t] = gen->P_std[0];
      gen->Q[t] = gen->Q[0];
    }
  }
}

