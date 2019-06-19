/** @file shunt.c
 *  @brief This file defines the Shunt data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/shunt.h>
#include <pfnet/bus.h>
#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct Shunt {
  
  // Bus
  Bus* bus;          /**< @brief Bus where the shunt is connected */
  Bus* reg_bus;      /**< @brief Bus regulated by this shunt */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Properties
  char type;                    /**< @brief Shunt type (fixed, switched, switched_v) */
  char mode;                    /**< @brief Shunt mode (discrete, continuous) */
  char name[SHUNT_BUFFER_SIZE]; /**< @brief Shunt name */

  // Conductance
  REAL g;            /**< @brief Conductance (p.u) */

  // Susceptance
  REAL* b;           /**< @brief Susceptance (p.u.) */
  REAL b_max;        /**< @brief Maximum susceptance (p.u.) */
  REAL b_min;        /**< @brief Minimum susceptance (p.u.) */
  REAL* b_values;    /**< @brief Array of valid susceptances (p.u.) */
  char num_b_values; /**< @brief Number of valid susceptances (p.u.) */
 
  // Flags
  BOOL in_service; /**< @brief Flag for indicating whether shunt is in service */
  char vars;       /**< @brief Flags for indicating which quantities are treated as variables **/
  char fixed;      /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;    /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;     /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Indices
  int index;      /**< @brief Shunt index */
  int* index_b;    /**< @brief Susceptance index */

  // Sensitivities
  REAL* sens_b_u_bound; /**< @brief Sensitivity of susceptance upper bound */
  REAL* sens_b_l_bound; /**< @brief Sensitivity of susceptance lower bound */

  // Transformer
  Branch* xfmr; /**< @brief Transformer associated with fixed shunt */

  // Network
  Net* net; /**< @brief Network. */

  // List
  Shunt* next;     /**< @brief List of shunts connceted to a bus */
  Shunt* reg_next; /**< @brief List of shunts regulated the same bus */
};

void* SHUNT_array_get(void* shunt_array, int index) { 
  if (shunt_array)
    return (void*)&(((Shunt*)shunt_array)[index]);
  else
    return NULL;
}

void SHUNT_array_del(Shunt* shunt_array, int size) {
  int i;
  Shunt* shunt;
  if (shunt_array) {
    for (i = 0; i < size; i++) {
      shunt = &(shunt_array[i]);
      free(shunt->b_values);
      free(shunt->b);
      free(shunt->index_b);
      free(shunt->sens_b_u_bound);
      free(shunt->sens_b_l_bound);
      SHUNT_set_bus(shunt,NULL);
      SHUNT_set_reg_bus(shunt,NULL);
    }
    free(shunt_array);
  }
}

Shunt* SHUNT_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    Shunt* shunt_array = (Shunt*)malloc(sizeof(Shunt)*size);
    for (i = 0; i < size; i++) {
      SHUNT_init(&(shunt_array[i]),num_periods);
      SHUNT_set_index(&(shunt_array[i]),i);
      snprintf(shunt_array[i].name,(size_t)(SHUNT_BUFFER_SIZE-1),"%d",i);
    }
    return shunt_array;
  }
  else
    return NULL;
}

void SHUNT_array_show(Shunt* shunt_array, int size, int t) { 
  int i;
  if (shunt_array) {
    for (i = 0; i < size; i++) 
      SHUNT_show(&(shunt_array[i]),t);
  }
}

void SHUNT_clear_sensitivities(Shunt* shunt) {
  int t;
  if (shunt) {
    for (t = 0; t < shunt->num_periods; t++) {
      shunt->sens_b_u_bound[t] = 0;
      shunt->sens_b_l_bound[t] = 0;
    }
  }
}

void SHUNT_clear_flags(Shunt* shunt, char flag_type) {
  if (shunt) {
    if (flag_type == FLAG_VARS)
      shunt->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      shunt->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      shunt->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      shunt->sparse = 0x00;
  }
}

void SHUNT_copy_from_shunt(Shunt* shunt, Shunt* other) {

  // Local variables
  int num_periods;

  // Check
  if (!shunt || !other)
    return;

  // Min num periods
  if (shunt->num_periods < other->num_periods)
    num_periods = shunt->num_periods;
  else
    num_periods = other->num_periods;

  // Bus
  // skip buses

  // Times
  // skip num periods

  // Properties
  shunt->type = other->type;
  shunt->mode = other->mode;
  strcpy(shunt->name,other->name);

  // Conductance
  shunt->g = other->g;

  // Susceptance
  memcpy(shunt->b,other->b,num_periods*sizeof(REAL));
  shunt->b_max = other->b_max;
  shunt->b_min = other->b_min;
  free(shunt->b_values);
  ARRAY_zalloc(shunt->b_values,REAL,other->num_b_values);
  memcpy(shunt->b_values,other->b_values,other->num_b_values*sizeof(REAL));
  shunt->num_b_values = other->num_b_values;

  // Flags
  SHUNT_set_in_service(shunt,SHUNT_is_in_service(other));
  shunt->fixed = other->fixed;
  shunt->bounded = other->bounded;
  shunt->sparse = other->sparse;
  shunt->vars = other->vars;

  // Indices
  // skip index
  memcpy(shunt->index_b,other->index_b,num_periods*sizeof(int));

  // Sensivitities
  memcpy(shunt->sens_b_u_bound,other->sens_b_u_bound,num_periods*sizeof(REAL));
  memcpy(shunt->sens_b_l_bound,other->sens_b_l_bound,num_periods*sizeof(REAL));

  // List
  // skip next
}

REAL SHUNT_get_sens_b_u_bound(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->sens_b_u_bound[t];
  else
    return 0;
}

REAL* SHUNT_get_sens_b_u_bound_array(Shunt* shunt) {
  if (shunt)
    return shunt->sens_b_u_bound;
  else
    return NULL;
}

REAL SHUNT_get_sens_b_l_bound(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->sens_b_l_bound[t];
  else
    return 0;
}

REAL* SHUNT_get_sens_b_l_bound_array(Shunt* shunt) {
  if (shunt)
    return shunt->sens_b_l_bound;
  else
    return NULL;
}

char SHUNT_get_flags_vars(Shunt* shunt) {
  if (shunt)
    return shunt->vars;
  else
    return 0;
}

char SHUNT_get_flags_fixed(Shunt* shunt) {
  if (shunt)
    return shunt->fixed;
  else
    return 0;
}

char SHUNT_get_flags_bounded(Shunt* shunt) {
  if (shunt)
    return shunt->bounded;
  else
    return 0;
}

char SHUNT_get_flags_sparse(Shunt* shunt) {
  if (shunt)
    return shunt->sparse;
  else
    return 0;
}

char SHUNT_get_type(Shunt* shunt) {
  if (shunt)
    return shunt->type;
  else
    return SHUNT_TYPE_FIXED;
}

char SHUNT_get_mode(Shunt* shunt) {
  if (shunt)
    return shunt->mode;
  else
    return SHUNT_MODE_CONT;
}

char* SHUNT_get_name(Shunt* shunt) {
  if (shunt)
    return shunt->name;
  else
    return NULL;
}

int SHUNT_get_num_periods(Shunt* shunt) {
  if (shunt)
    return shunt->num_periods;
  else
    return 0;
}

char SHUNT_get_obj_type(void* shunt) {
  if (shunt)
    return OBJ_SHUNT;
  else
    return OBJ_UNKNOWN;
}

int SHUNT_get_index(Shunt* shunt) {
  if (shunt)
    return shunt->index;
  else
    return -1;
}

int SHUNT_get_index_b(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->index_b[t];
  else
    return -1;
}

int* SHUNT_get_index_b_array(Shunt* shunt) {
  if (shunt)
    return shunt->index_b;
  else
    return NULL;
}

Bus* SHUNT_get_bus(Shunt* shunt) {
  if (shunt)
    return shunt->bus;
  else
    return NULL;
}
Bus* SHUNT_get_reg_bus(Shunt* shunt) {
  if (shunt)
    return shunt->reg_bus;
  else
    return NULL;
}

REAL SHUNT_get_g(Shunt* shunt) {
  if (shunt)
    return shunt->g;
  else
    return 0;
}

REAL SHUNT_get_b(Shunt* shunt, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    return shunt->b[t];
  else
    return 0;
}

REAL* SHUNT_get_b_array(Shunt* shunt) {
  if (shunt)
    return shunt->b;
  else
    return NULL;
}

REAL SHUNT_get_b_max(Shunt* shunt) {
  if (shunt)
    return shunt->b_max;
  else
    return 0;
}

REAL SHUNT_get_b_min(Shunt* shunt) {
  if (shunt)
    return shunt->b_min;
  else
    return 0;
}

REAL* SHUNT_get_b_values(Shunt* shunt) {
  if (shunt)
    return shunt->b_values;
  else
    return NULL;
}

int SHUNT_get_num_b_values(Shunt* shunt) {
  if (shunt)
    return shunt->num_b_values;
  else
    return 0;
}


Shunt* SHUNT_get_next(Shunt* shunt) {
  if (shunt)
    return shunt->next;
  else
    return NULL;
}

Shunt* SHUNT_get_reg_next(Shunt* shunt) {
  if (shunt)
    return shunt->reg_next;
  else
    return NULL;
}

void SHUNT_get_var_values(Shunt* shunt, Vec* values, int code) {

  // Local vars
  int t;

  // No shunt
  if (!shunt)
    return;

  // Time loop
  for (t = 0; t < shunt->num_periods; t++) {
    
    if (shunt->vars & SHUNT_VAR_SUSC) { // susceptance
      switch(code) {

      case UPPER_LIMITS:
        if (shunt->bounded & SHUNT_VAR_SUSC)
          VEC_set(values,shunt->index_b[t],shunt->b_max);
        else
          VEC_set(values,shunt->index_b[t],SHUNT_INF_SUSC);
        break;

      case LOWER_LIMITS:
        if (shunt->bounded & SHUNT_VAR_SUSC)
          VEC_set(values,shunt->index_b[t],shunt->b_min);
        else
          VEC_set(values,shunt->index_b[t],-SHUNT_INF_SUSC);
        break;

      default:
        VEC_set(values,shunt->index_b[t],shunt->b[t]);
      }
    }
  }
}

char* SHUNT_get_var_info_string(Shunt* shunt, int index) {

  // Local variables
  char* info;

  //Check
  if (!shunt)
    return NULL;

  // Susceptance
  if ((shunt->vars & SHUNT_VAR_SUSC) &&
      index >= shunt->index_b[0] &&
      index <= shunt->index_b[shunt->num_periods-1]) {
    info = (char*)malloc(SHUNT_BUFFER_SIZE*sizeof(char));
    snprintf(info,SHUNT_BUFFER_SIZE*sizeof(char),
             "shunt:%d:susceptance:%d",shunt->index,index-shunt->index_b[0]);
    return info;
  }
  // Return
  return NULL;
}

int SHUNT_get_num_vars(void* vshunt, unsigned char var, int t_start, int t_end) {

  // Local vars
  Shunt* shunt = (Shunt*)vshunt;
  int num_vars = 0;
  int dt;

  // Checks
  if (!shunt)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > shunt->num_periods-1)
    t_end = shunt->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & SHUNT_VAR_SUSC) && (shunt->vars & SHUNT_VAR_SUSC))
    num_vars += dt;
  return num_vars;
}

Vec* SHUNT_get_var_indices(void* vshunt, unsigned char var, int t_start, int t_end) {

  // Local vars
  Shunt* shunt = (Shunt*)vshunt;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!shunt)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > shunt->num_periods-1)
    t_end = shunt->num_periods-1;

  // Indices
  indices = VEC_new(SHUNT_get_num_vars(vshunt,var,t_start,t_end));
  if ((var & SHUNT_VAR_SUSC) && (shunt->vars & SHUNT_VAR_SUSC)) {
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,shunt->index_b[t]);
      offset++;
    }
  }
  return indices;
}

char* SHUNT_get_json_string(Shunt* shunt, char* output) {

  // Local variables
  char temp[SHUNT_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No shunt
  if (!shunt)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*SHUNT_BUFFER_SIZE*SHUNT_NUM_JSON_FIELDS*shunt->num_periods);
    resize = TRUE;
  }
  output_start = output;

  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",shunt->index,FALSE);
  JSON_obj(temp,output,"bus",shunt->bus,BUS_get_index,FALSE);
  JSON_obj(temp,output,"reg_bus",shunt->reg_bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"num_periods",shunt->num_periods,FALSE);
  JSON_str(temp,output,"name",shunt->name,FALSE);
  JSON_bool(temp,output,"in_service",SHUNT_is_in_service(shunt),FALSE);
  JSON_int(temp,output,"type",shunt->type,FALSE);
  JSON_int(temp,output,"mode",shunt->mode,FALSE);
  JSON_float(temp,output,"g",shunt->g,FALSE);
  JSON_array_float(temp,output,"b",shunt->b,shunt->num_periods,FALSE);
  JSON_float(temp,output,"b_max",shunt->b_max,FALSE);
  JSON_float(temp,output,"b_min",shunt->b_min,FALSE);
  JSON_array_float(temp,output,"b_values",shunt->b_values,shunt->num_b_values,TRUE);
  JSON_end(output);
  
  // Output
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL SHUNT_has_flags(void* vshunt, char flag_type, unsigned char mask) {
  Shunt* shunt = (Shunt*)vshunt;
  if (shunt) {
    if (flag_type == FLAG_VARS)
      return (shunt->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (shunt->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (shunt->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (shunt->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL SHUNT_has_properties(void* vshunt, char prop) {
  Shunt* shunt = (Shunt*)vshunt;
  if (!shunt)
    return FALSE;
  if ((prop & SHUNT_PROP_SWITCHED_V) && !SHUNT_is_switched_v(shunt))
    return FALSE;
  return TRUE;
}

void SHUNT_init(Shunt* shunt, int num_periods) {

  // Local vars
  int T;
  
  // No gen
  if (!shunt)
    return;

  T = num_periods;
  shunt->num_periods = num_periods;

  shunt->type = SHUNT_TYPE_FIXED;
  shunt->mode = SHUNT_MODE_CONT;
  ARRAY_clear(shunt->name,char,SHUNT_BUFFER_SIZE);
  
  shunt->bus = NULL;
  shunt->reg_bus = NULL;
  shunt->g = 0;
  shunt->b_max = 0;
  shunt->b_min = 0;
  shunt->b_values = NULL;
  shunt->num_b_values = 0;

  shunt->in_service = TRUE;
  shunt->vars = 0x00;
  shunt->fixed = 0x00;
  shunt->bounded = 0x00;
  shunt->sparse = 0x00;
  shunt->index = -1;

  ARRAY_zalloc(shunt->b,REAL,T);
  ARRAY_zalloc(shunt->index_b,int,T);
  ARRAY_zalloc(shunt->sens_b_u_bound,REAL,T);
  ARRAY_zalloc(shunt->sens_b_l_bound,REAL,T);

  shunt->net = NULL;
  shunt->xfmr = NULL;

  shunt->next = NULL;
  shunt->reg_next = NULL;
}

BOOL SHUNT_is_part_of_transformer(Shunt* shunt) {
  if (shunt)
    return shunt->xfmr != NULL;
  else
    return FALSE;
}

BOOL SHUNT_is_in_service(void* shunt) {
  if (shunt)
    return ((Shunt*)shunt)->in_service && BUS_is_in_service(((Shunt*)shunt)->bus);
  else
    return FALSE;
}

BOOL SHUNT_is_equal(Shunt* shunt, Shunt* other) {
  return shunt == other;
}

BOOL SHUNT_is_fixed(Shunt* shunt) {
  if (shunt)
    return shunt->type == SHUNT_TYPE_FIXED;
  else
    return FALSE;
}

BOOL SHUNT_is_switched(Shunt* shunt) {
  if (shunt)
    return (shunt->type == SHUNT_TYPE_SWITCHED ||
            shunt->type == SHUNT_TYPE_SWITCHED_V);
  else
    return FALSE;
}

BOOL SHUNT_is_switched_v(Shunt* shunt) {
  if (shunt)
    return shunt->type == SHUNT_TYPE_SWITCHED_V;
  else
    return FALSE;
}

BOOL SHUNT_is_switched_locked(Shunt* shunt) {
  if (shunt)
    return shunt->type == SHUNT_TYPE_SWITCHED;
  else
    return FALSE;
}

BOOL SHUNT_is_continuous(Shunt* shunt) {
  if (shunt)
    return shunt->mode == SHUNT_MODE_CONT;
  else
    return FALSE;
}

BOOL SHUNT_is_discrete(Shunt* shunt) {
  if (shunt)
    return shunt->mode == SHUNT_MODE_DIS;
  else
    return FALSE;
}

Shunt* SHUNT_list_add(Shunt* shunt_list, Shunt* shunt) {
  LIST_add(Shunt,shunt_list,shunt,next);
  return shunt_list;
}

Shunt* SHUNT_list_del(Shunt* shunt_list, Shunt* shunt) {
  LIST_del(Shunt,shunt_list,shunt,next);
  return shunt_list;
}

int SHUNT_list_len(Shunt* shunt_list) {
  int len;
  LIST_len(Shunt,shunt_list,next,len);
  return len;
}

Shunt* SHUNT_list_reg_add(Shunt* reg_shunt_list, Shunt* reg_shunt) {
  LIST_add(Shunt,reg_shunt_list,reg_shunt,reg_next);
  return reg_shunt_list;
}

Shunt* SHUNT_list_reg_del(Shunt* reg_shunt_list, Shunt* reg_shunt) {
  LIST_del(Shunt,reg_shunt_list,reg_shunt,reg_next);
  return reg_shunt_list;
}

int SHUNT_list_reg_len(Shunt* reg_shunt_list) {
  int len;
  LIST_len(Shunt,reg_shunt_list,reg_next,len);
  return len;
}

Shunt* SHUNT_new(int num_periods) { 
  if (num_periods > 0) {
    Shunt* shunt = (Shunt*)malloc(sizeof(Shunt));
    SHUNT_init(shunt,num_periods);
    return shunt;
  }
  else
    return NULL;
}

int SHUNT_round_b(Shunt* shunt, int t) {
  int i;
  int ibest;
  REAL db = 0;
  if (shunt && shunt->b_values && shunt->num_b_values > 0) {
    ibest = 0;
    db = fabs(SHUNT_get_b(shunt,t)-shunt->b_values[0]);
    for (i = 1; i < shunt->num_b_values; i++) {
      if (fabs(SHUNT_get_b(shunt,t)-shunt->b_values[i]) <= db) {
        ibest = i;
        db = fabs(SHUNT_get_b(shunt,t)-shunt->b_values[i]);
      }
    }
    SHUNT_set_b(shunt,shunt->b_values[ibest],t);
  }
  return (db > SHUNT_ROUND_THRESHOLD) ? 1 : 0;
}

void SHUNT_set_sens_b_u_bound(Shunt* shunt, REAL value, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    shunt->sens_b_u_bound[t] = value;
}

void SHUNT_set_sens_b_l_bound(Shunt* shunt, REAL value, int t) {
  if (shunt && t >= 0 && t < shunt->num_periods)
    shunt->sens_b_l_bound[t] = value;
}

void SHUNT_set_in_service(Shunt* shunt, BOOL in_service) {
  if (shunt) {
    if (shunt->in_service != in_service)
      NET_inc_state_tag(shunt->net);
    shunt->in_service = in_service;
  }
}

void SHUNT_set_network(Shunt* shunt, void* net) {
  if (shunt)
    shunt->net = (Net*)net;
}

void SHUNT_set_type(Shunt* shunt, char type) {
  if (shunt)
    shunt->type = type;
}

void SHUNT_set_mode(Shunt* shunt, char mode) {
  if (shunt)
    shunt->mode = mode;
}

void SHUNT_set_name(Shunt* shunt, char* name) {
  if (shunt)
    strncpy(shunt->name,name,(size_t)(SHUNT_BUFFER_SIZE-1));
}

void SHUNT_set_bus(Shunt* shunt, Bus* bus) {
  Bus* old_bus;
  if (shunt) {
    old_bus = shunt->bus;
    shunt->bus = NULL;
    BUS_del_shunt(old_bus,shunt);
    shunt->bus = bus;
    BUS_add_shunt(shunt->bus,shunt);
  }
}

void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus) {
  Bus* old_reg_bus;
  if (shunt) {
    old_reg_bus = shunt->reg_bus;
    shunt->reg_bus = NULL;
    BUS_del_reg_shunt(old_reg_bus,shunt);
    shunt->reg_bus = reg_bus;
    BUS_add_reg_shunt(shunt->reg_bus,shunt);
    if (reg_bus)
      shunt->type = SHUNT_TYPE_SWITCHED_V;
    else if (shunt->type == SHUNT_TYPE_SWITCHED_V)
      shunt->type = SHUNT_TYPE_SWITCHED;
  }
}

void SHUNT_set_index(Shunt* shunt, int index) { 
  if (shunt)
    shunt->index = index;
}

void SHUNT_set_g(Shunt* shunt, REAL g) { 
  if (shunt)
    shunt->g = g;
}

void SHUNT_set_b(Shunt* shunt, REAL b, int t) { 
  if (shunt && t >= 0 && t < shunt->num_periods)
    shunt->b[t] = b;
}

void SHUNT_set_b_max(Shunt* shunt, REAL b_max) {
  if (shunt)
    shunt->b_max = b_max;
}

void SHUNT_set_b_min(Shunt* shunt, REAL b_min) {
  if (shunt)
    shunt->b_min = b_min;
}

void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num) {
  if (shunt) {
    if (shunt->b_values)
      free(shunt->b_values);
    shunt->b_values = values;
    shunt->num_b_values = num;
  }
}

int SHUNT_set_flags(void* vshunt, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Shunt* shunt = (Shunt*)vshunt;
  int t;

  // No shunt
  if (!shunt)
    return index;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(shunt->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(shunt->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(shunt->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(shunt->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & SHUNT_VAR_SUSC) && (mask & SHUNT_VAR_SUSC)) { // shunt susceptance
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < shunt->num_periods; t++)
        shunt->index_b[t] = index+t;
    }
    (*flags_ptr) |= SHUNT_VAR_SUSC;
    index += shunt->num_periods;
  }
  return index;
}

void SHUNT_set_var_values(Shunt* shunt, Vec* values) {

  // Local vars
  int t;

  // No shunt
  if (!shunt)
    return;

  // Time loop
  for (t = 0; t < shunt->num_periods; t++) {

    if (shunt->vars & SHUNT_VAR_SUSC) // shunt susceptance (p.u.)
      shunt->b[t] = VEC_get(values,shunt->index_b[t]); 
  }
}

void SHUNT_set_xfmr(Shunt* shunt, Branch* xfmr) {
  if (shunt)
    shunt->xfmr = xfmr;
}

void SHUNT_show(Shunt* shunt, int t) { 
  if (shunt)
    printf("shunt %d\t%d\n",
           BUS_get_number(shunt->bus),
           shunt->index);
}

void SHUNT_propagate_data_in_time(Shunt* shunt, int start, int end) {
  int t;
  if (shunt) {
    if (start < 0)
      start = 0;
    if (end > shunt->num_periods)
      end = shunt->num_periods;
    for (t = start+1; t < end; t++)
      shunt->b[t] = shunt->b[start];
  }
}
