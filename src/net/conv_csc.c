/** @file conv_csc.c
 *  @brief This file defines the ConvCSC data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/conv_csc.h>
#include <pfnet/bus.h>
#include <pfnet/bus_dc.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct ConvCSC {
  
  // Properties
  char type;           /**< @brief Converter type (rectifier or inverter) */
  char name[CONVCSC_BUFFER_SIZE]; /**< @brief Converter name */
  
  // Times
  int num_periods;     /**< @brief Number of time periods. */

  // Buses
  Bus* ac_bus;         /**< @brief AC bus to which the converter is connected */
  BusDC* dc_bus;       /**< @brief DC bus to which the converter is connected */

  // Flags
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Control
  char mode_dc;        /**< @brief Control mode */
  REAL* P_dc_set;      /**< @brief DC power set point (p.u. system base power) */
  REAL* i_dc_set;      /**< @brief DC current set point (p.u.) */
  REAL* v_dc_set;      /**< @brief DC voltage set point (p.u.) */
  
  // AC injections
  REAL* P;             /**< @brief Active power injection into AC bus (p.u. system base power) */
  REAL* Q;             /**< @brief Rective power injection into AC bus (p.u. system base power) */

  // DC injections
  REAL* P_dc;          /**< @brief Power injection into DC bus (p.u. system base power) */

  // Bridges
  int num_bridges;    /**< @brief Number of bridges in series */
  
  // Commutating capacitor
  REAL x_cap;          /**< @brief Commutating capacitor reactance as seen by each individual bridge (p.u.) */ 
  
  // Commutating transformer
  REAL x;              /**< @brief Commutating transformer reactance as seen by each individual bridge (p.u.) */
  REAL r;              /**< @brief Commutating transformer resistance as seen by each individual bridge (p.u.) */
  REAL* ratio;         /**< @brief Commutating transformer turns ratio (p.u.) */
  REAL ratio_max;      /**< @brief Maximum commutating transformer turns ratio (p.u.) */
  REAL ratio_min;      /**< @brief Minimum commutating transformer turns ratio (p.u.) */
  
  // Angle
  REAL* angle;         /**< @brief Ignition delay angle for rectifier or extinction advance angle for inverter (radians) */
  REAL angle_max;      /**< @brief Nominal maximum angle (radians) */
  REAL angle_min;      /**< @brief Minimum steady-state angle (radians) */
  
  // Base voltages
  REAL v_base_p;       /**< @brief Primary side bus base AC voltage (kilovolts) */
  REAL v_base_s;       /**< @brief Secondary side bus base AC voltage (kilovolts) */

  // Indices
  int index;           /**< @brief Converter index */
  int* index_P;        /**< @brief Active power injection index */
  int* index_Q;        /**< @brief Reactive power injection index */
  int* index_P_dc;     /**< @brief DC power injection index */
  int* index_i_dc;     /**< @brief DC current injection index */
  int* index_ratio;    /**< @brief Commutating transformer turns ratio index */
  int* index_angle;    /**< @brief Ignition delay or extinction advance angle index */
  
  // List
  ConvCSC* next_ac;       /**< @brief List of converters connected to an AC bus */
  ConvCSC* next_dc;       /**< @brief List of converters connected to a DC bus */
};

void CONVCSC_array_del(ConvCSC* conv_array, int size) {
  int i;
  ConvCSC* conv;
  if (conv_array) {
    for (i = 0; i < size; i++) {
      conv = &(conv_array[i]);
      free(conv->P_dc_set);
      free(conv->i_dc_set);
      free(conv->v_dc_set);
      free(conv->P);
      free(conv->Q);
      free(conv->P_dc);
      free(conv->ratio);
      free(conv->angle);
      free(conv->index_P);
      free(conv->index_Q);
      free(conv->index_P_dc);
      free(conv->index_i_dc);
      free(conv->index_ratio);
      free(conv->index_angle);
      CONVCSC_set_ac_bus(conv,NULL);
      CONVCSC_set_dc_bus(conv,NULL);
    }
    free(conv_array);
  }  
}

void* CONVCSC_array_get(void* conv_array, int index) { 
  if (conv_array) 
    return (void*)&(((ConvCSC*)conv_array)[index]);
  else
    return NULL;
}

ConvCSC* CONVCSC_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    ConvCSC* conv_array = (ConvCSC*)malloc(sizeof(ConvCSC)*size);
    for (i = 0; i < size; i++) {
      CONVCSC_init(&(conv_array[i]),num_periods);
      CONVCSC_set_index(&(conv_array[i]),i);
      snprintf(conv_array[i].name,(size_t)(CONVCSC_BUFFER_SIZE-1),"%d",i);
    }
    return conv_array;
  }
  else
    return NULL;
}

void CONVCSC_clear_flags(ConvCSC* conv, char flag_mask) {
  if (conv) {
    if (flag_mask & FLAG_VARS)
      conv->vars = 0x00;
    if (flag_mask & FLAG_BOUNDED)
      conv->bounded = 0x00;
    if (flag_mask & FLAG_FIXED)
      conv->fixed = 0x00;
    if (flag_mask & FLAG_SPARSE)
      conv->sparse = 0x00;
  }
}

void CONVCSC_clear_sensitivities(ConvCSC* conv) {
  // Nothing
}

void CONVCSC_copy_from_conv(ConvCSC* conv, ConvCSC* other) {

  // Local variables
  int num_periods;

  // Check
  if (!conv || !other)
    return;

  // Min num periods
  if (conv->num_periods < other->num_periods)
    num_periods = conv->num_periods;
  else
    num_periods = other->num_periods;

  // Bus ac
  // skip bus ac

  // Bus dc
  // skip bus dc

  // Times
  // skip num periods

  // Properties
  conv->type = other->type;
  strcpy(conv->name,other->name);

  // Flags
  conv->fixed = other->fixed;
  conv->bounded = other->bounded;
  conv->sparse = other->sparse;
  conv->vars = other->vars;

  // Control
  conv->mode_dc = other->mode_dc;
  memcpy(conv->P_dc_set,other->P_dc_set,num_periods*sizeof(REAL));
  memcpy(conv->i_dc_set,other->i_dc_set,num_periods*sizeof(REAL));
  memcpy(conv->v_dc_set,other->v_dc_set,num_periods*sizeof(REAL));
  
  // AC injections
  memcpy(conv->P,other->P,num_periods*sizeof(REAL));
  memcpy(conv->Q,other->Q,num_periods*sizeof(REAL));

  // DC injection
  memcpy(conv->P_dc,other->P_dc,num_periods*sizeof(REAL));

  // Num bridges
  conv->num_bridges = other->num_bridges;

  // Commutating capacitor
  conv->x_cap = other->x_cap;

  // Commutating transformer
  conv->x = other->x;
  conv->r = other->r;
  memcpy(conv->ratio,other->ratio,num_periods*sizeof(REAL));
  conv->ratio_max = other->ratio_max;
  conv->ratio_min = other->ratio_min;

  // Angle
  memcpy(conv->angle,other->angle,num_periods*sizeof(REAL));
  conv->angle_max = other->angle_max;
  conv->angle_min = other->angle_min;

  // Base voltages
  conv->v_base_p = other->v_base_p;
  conv->v_base_s = other->v_base_s;
  
  // Indices
  // skip index
  memcpy(conv->index_P,other->index_P,num_periods*sizeof(int));
  memcpy(conv->index_Q,other->index_Q,num_periods*sizeof(int));
  memcpy(conv->index_P_dc,other->index_P_dc,num_periods*sizeof(int));
  memcpy(conv->index_i_dc,other->index_i_dc,num_periods*sizeof(int));
  memcpy(conv->index_ratio,other->index_ratio,num_periods*sizeof(int));
  memcpy(conv->index_angle,other->index_angle,num_periods*sizeof(int));
    
  // List
  // skip next
}

Bus* CONVCSC_get_ac_bus(ConvCSC* conv) {
  if (conv)
    return conv->ac_bus;
  else
    return NULL;
}

BusDC* CONVCSC_get_dc_bus(ConvCSC* conv) {
  if (conv)
    return conv->dc_bus;
  else
    return NULL;
}

REAL CONVCSC_get_angle(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->angle[t];
  else 
    return 0;
}

REAL CONVCSC_get_angle_max(ConvCSC* conv) {
  if (conv)
    return conv->angle_max;
  else
    return 0;
}

REAL CONVCSC_get_angle_min(ConvCSC* conv) {
  if (conv)
    return conv->angle_min;
  else
    return 0;
}

char CONVCSC_get_flags_vars(ConvCSC* conv) {
  if (conv)
    return conv->vars;
  else
    return 0;
}

char CONVCSC_get_flags_fixed(ConvCSC* conv) {
  if (conv)
    return conv->fixed;
  else
    return 0;
}

char CONVCSC_get_flags_bounded(ConvCSC* conv) {
  if (conv)
    return conv->bounded;
  else
    return 0;
}

char CONVCSC_get_flags_sparse(ConvCSC* conv) {
  if (conv)
    return conv->sparse;
  else
    return 0;
}

REAL CONVCSC_get_i_dc_set(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->i_dc_set[t];
  else 
    return 0;
}

int CONVCSC_get_index(ConvCSC* conv) {
  if (conv)
    return conv->index;
  else
    return -1;
}

int CONVCSC_get_index_P(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_P[t];
  else
    return -1;
}

int CONVCSC_get_index_Q(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_Q[t];
  else
    return -1;
}

int CONVCSC_get_index_P_dc(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_P_dc[t];
  else
    return -1;
}

int CONVCSC_get_index_i_dc(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_i_dc[t];
  else
    return -1;
}

int CONVCSC_get_index_ratio(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_ratio[t];
  else
    return -1;
}

int CONVCSC_get_index_angle(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_ratio[t];
  else
    return -1;
}

char* CONVCSC_get_json_string(ConvCSC* conv, char* output) {

  // Local variables
  char temp[CONVCSC_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No conv
  if (!conv)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*CONVCSC_BUFFER_SIZE*CONVCSC_NUM_JSON_FIELDS*conv->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",conv->index,FALSE);
  JSON_obj(temp,output,"ac_bus",conv->ac_bus,BUS_get_index,FALSE);
  JSON_obj(temp,output,"dc_bus",conv->dc_bus,BUSDC_get_index,FALSE);
  JSON_int(temp,output,"num_periods",conv->num_periods,FALSE);
  JSON_str(temp,output,"name",conv->name,FALSE);
  JSON_int(temp,output,"type",conv->type,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

char CONVCSC_get_mode_dc(ConvCSC* conv) {
  if (conv)
    return conv->mode_dc;
  else
    return CONVCSC_MODE_DC_NC;
}

char* CONVCSC_get_name(ConvCSC* conv) {
  if (conv)
    return conv->name;
  else
    return NULL;
}

ConvCSC* CONVCSC_get_next_ac(ConvCSC* conv) {
  if (conv)
    return conv->next_ac;
  else
    return NULL;
}

ConvCSC* CONVCSC_get_next_dc(ConvCSC* conv) {
  if (conv)
    return conv->next_dc;
  else
    return NULL;
}

int CONVCSC_get_num_bridges(ConvCSC* conv) {
  if (conv)
    return conv->num_bridges;
  else
    return 0;
}

int CONVCSC_get_num_periods(ConvCSC* conv) {
  if (conv)
    return conv->num_periods;
  else
    return 0;
}

int CONVCSC_get_num_vars(void* vconv, unsigned char var, int t_start, int t_end) {

  // Local vars
  ConvCSC* conv = (ConvCSC*)vconv;
  int num_vars = 0;
  int dt;

  // Checks
  if (!conv)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > conv->num_periods-1)
    t_end = conv->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & CONVCSC_VAR_P) && (conv->vars & CONVCSC_VAR_P)) // active power
    num_vars += dt;
  if ((var & CONVCSC_VAR_Q) && (conv->vars & CONVCSC_VAR_Q)) // reactive power
    num_vars += dt;
  if ((var & CONVCSC_VAR_PDC) && (conv->vars & CONVCSC_VAR_PDC)) // DC power
    num_vars += dt;
  if ((var & CONVCSC_VAR_RATIO) && (conv->vars & CONVCSC_VAR_RATIO)) // taps ratio
    num_vars += dt;
  if ((var & CONVCSC_VAR_ANGLE) && (conv->vars & CONVCSC_VAR_ANGLE)) // ignition or extinction angle
    num_vars += dt;
  return num_vars;
}

char CONVCSC_get_obj_type(void* conv) {
  if (conv)
    return OBJ_CONVCSC;
  else
    return OBJ_UNKNOWN;
}

REAL CONVCSC_get_P(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P[t];
  else 
    return 0;
}

REAL CONVCSC_get_P_dc(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc[t];
  else 
    return 0;
}

REAL CONVCSC_get_i_dc(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc[t]/(BUSDC_get_v(conv->dc_bus,t));
  else 
    return 0;
}

REAL CONVCSC_get_P_dc_set(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc_set[t];
  else 
    return 0;
}

REAL CONVCSC_get_Q(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->Q[t];
  else 
    return 0;
}

REAL CONVCSC_get_ratio(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->ratio[t];
  else 
    return 0;
}

REAL CONVCSC_get_ratio_max(ConvCSC* conv) {
  if (conv)
    return conv->ratio_max;
  else
    return 0;
}

REAL CONVCSC_get_ratio_min(ConvCSC* conv) {
  if (conv)
    return conv->ratio_min;
  else
    return 0;
}

REAL CONVCSC_get_v_base_p(ConvCSC* conv) {
  if (conv)
    return conv->v_base_p;
  else
    return 0;
}

REAL CONVCSC_get_v_base_s(ConvCSC* conv) {
  if (conv)
    return conv->v_base_s;
  else
    return 0;
}

REAL CONVCSC_get_v_dc_set(ConvCSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->v_dc_set[t];
  else 
    return 0;
}

Vec* CONVCSC_get_var_indices(void* vconv, unsigned char var, int t_start, int t_end) {

  // Local vars
  ConvCSC* conv = (ConvCSC*)vconv;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!conv)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > conv->num_periods-1)
    t_end = conv->num_periods-1;

  // Indices
  indices = VEC_new(CONVCSC_get_num_vars(vconv,var,t_start,t_end));
  if ((var & CONVCSC_VAR_P) && (conv->vars & CONVCSC_VAR_P)) { // active power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_P[t]);
      offset++;
    }
  }
  if ((var & CONVCSC_VAR_Q) && (conv->vars & CONVCSC_VAR_Q)) { // reactive power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_Q[t]);
      offset++;
    }
  }
  if ((var & CONVCSC_VAR_PDC) && (conv->vars & CONVCSC_VAR_PDC)) { // DC power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_P_dc[t]);
      VEC_set(indices,offset+1,conv->index_i_dc[t]);
      offset += 2;
    }
  }
  if ((var & CONVCSC_VAR_RATIO) && (conv->vars & CONVCSC_VAR_RATIO)) { // taps ratio
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_ratio[t]);
      offset++;
    }
  }
  if ((var & CONVCSC_VAR_ANGLE) && (conv->vars & CONVCSC_VAR_ANGLE)) { // ignition or extinction angle
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_angle[t]);
      offset++;
    }
  }
  
  return indices;
}

char* CONVCSC_get_var_info_string(ConvCSC* conv, int index) {

  // Local variables
  char* info;

  //Check
  if (!conv)
    return NULL;

  // Active power
  if ((conv->vars & CONVCSC_VAR_P) &&
      index >= conv->index_P[0] &&
      index <= conv->index_P[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:active power:%d",conv->index,index-conv->index_P[0]);
    return info;
  }

  // Reactive power
  if ((conv->vars & CONVCSC_VAR_Q) &&
      index >= conv->index_Q[0] &&
      index <= conv->index_Q[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:reactive power:%d",conv->index,index-conv->index_Q[0]);
    return info;
  }

  // DC power
  if ((conv->vars & CONVCSC_VAR_PDC) &&
      index >= conv->index_P_dc[0] &&
      index <= conv->index_P_dc[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:dc power:%d",conv->index,index-conv->index_P_dc[0]);
    return info;
  }

  // DC current
  if ((conv->vars & CONVCSC_VAR_PDC) &&
      index >= conv->index_i_dc[0] &&
      index <= conv->index_i_dc[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:dc current:%d",conv->index,index-conv->index_i_dc[0]);
    return info;
  }

  // Taps ratio
  if ((conv->vars & CONVCSC_VAR_RATIO) &&
      index >= conv->index_ratio[0] &&
      index <= conv->index_ratio[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:tap ratio:%d",conv->index,index-conv->index_ratio[0]);
    return info;
  }

  // Ignition or extinction angle
  if ((conv->vars & CONVCSC_VAR_ANGLE) &&
      index >= conv->index_angle[0] &&
      index <= conv->index_angle[conv->num_periods-1]) {
    info = (char*)malloc(CONVCSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVCSC_BUFFER_SIZE*sizeof(char),
             "csc_conv:%d:angle:%d",conv->index,index-conv->index_angle[0]);
    return info;
  }

  // Return
  return NULL;
}

void CONVCSC_get_var_values(ConvCSC* conv, Vec* values, int code) {
 
  // Local vars
  int t;
 
  // No laod
  if (!conv)
    return;

  // Time loop
  for (t = 0; t < conv->num_periods; t++) {

    if (conv->vars & CONVCSC_VAR_P) { // active power
      switch(code) {
      case UPPER_LIMITS:
        VEC_set(values,conv->index_P[t],CONVCSC_INF_P);
        break;
      case LOWER_LIMITS:
        VEC_set(values,conv->index_P[t],-CONVCSC_INF_P);
        break;
      default:
        VEC_set(values,conv->index_P[t],conv->P[t]);
      }
    }

    if (conv->vars & CONVCSC_VAR_Q) { // reactive power
      switch(code) {
      case UPPER_LIMITS:
        VEC_set(values,conv->index_Q[t],CONVCSC_INF_Q);
        break;
      case LOWER_LIMITS:
        VEC_set(values,conv->index_Q[t],-CONVCSC_INF_Q);
        break;
      default:
        VEC_set(values,conv->index_Q[t],conv->Q[t]);
      }
    }

    if (conv->vars & CONVCSC_VAR_PDC) { // DC power
      switch(code) {
      case UPPER_LIMITS:
        VEC_set(values,conv->index_P_dc[t],CONVCSC_INF_PDC);
        VEC_set(values,conv->index_i_dc[t],CONVCSC_INF_PDC);
        break;
      case LOWER_LIMITS:
        VEC_set(values,conv->index_P_dc[t],-CONVCSC_INF_PDC);
        VEC_set(values,conv->index_i_dc[t],-CONVCSC_INF_PDC);
        break;
      default:
        VEC_set(values,conv->index_P_dc[t],conv->P_dc[t]);
        VEC_set(values,conv->index_i_dc[t],conv->P_dc[t]/(BUSDC_get_v(conv->dc_bus,t)));
      }
    }

    if (conv->vars & CONVCSC_VAR_RATIO) { // tap ratio
      switch(code) {
      case UPPER_LIMITS:
        VEC_set(values,conv->index_ratio[t],CONVCSC_INF_RATIO);
        break;
      case LOWER_LIMITS:
        VEC_set(values,conv->index_ratio[t],-CONVCSC_INF_RATIO);
        break;
      default:
        VEC_set(values,conv->index_ratio[t],conv->ratio[t]);
      }
    }

    if (conv->vars & CONVCSC_VAR_ANGLE) { // ignition or extinction angle
      switch(code) {
      case UPPER_LIMITS:
        VEC_set(values,conv->index_angle[t],CONVCSC_INF_ANGLE);
        break;
      case LOWER_LIMITS:
        VEC_set(values,conv->index_angle[t],-CONVCSC_INF_ANGLE);
        break;
      default:
        VEC_set(values,conv->index_angle[t],conv->angle[t]);
      }
    }
  }
}

REAL CONVCSC_get_x_cap(ConvCSC* conv) {
  if (conv)
    return conv->x_cap;
  else
    return 0;
}

REAL CONVCSC_get_x(ConvCSC* conv) {
  if (conv)
    return conv->x;
  else
    return 0;
}

REAL CONVCSC_get_r(ConvCSC* conv) {
  if (conv)
    return conv->r;
  else
    return 0;
}

void CONVCSC_init(ConvCSC* conv, int num_periods) {

  // Local variables
  int T;
  
  // No conv
  if (!conv)
    return;

  conv->type = CONVCSC_TYPE_UNKNOWN;
  ARRAY_clear(conv->name,char,CONVCSC_BUFFER_SIZE);

  T = num_periods;
  conv->num_periods = num_periods;

  conv->ac_bus = NULL;
  conv->dc_bus = NULL;

  conv->fixed = 0x00;
  conv->bounded = 0x00;
  conv->sparse = 0x00;
  conv->vars = 0x00;

  conv->mode_dc = CONVCSC_MODE_DC_NC;
  ARRAY_zalloc(conv->P_dc_set,REAL,T);
  ARRAY_zalloc(conv->i_dc_set,REAL,T);
  ARRAY_zalloc(conv->v_dc_set,REAL,T);

  ARRAY_zalloc(conv->P,REAL,T);
  ARRAY_zalloc(conv->Q,REAL,T);

  ARRAY_zalloc(conv->P_dc,REAL,T);

  conv->num_bridges = 1;

  conv->x_cap = 0;

  conv->x = 0;
  conv->r = 0;

  ARRAY_zalloc(conv->ratio,REAL,T);
  conv->ratio_max = 1.;
  conv->ratio_min = 1.;
  
  ARRAY_zalloc(conv->angle,REAL,T);
  conv->angle_max = 0;
  conv->angle_min = 0;

  conv->v_base_p = 0;
  conv->v_base_s = 0;

  conv->index = -1;
  ARRAY_zalloc(conv->index_P,int,T);
  ARRAY_zalloc(conv->index_Q,int,T);
  ARRAY_zalloc(conv->index_P_dc,int,T);
  ARRAY_zalloc(conv->index_i_dc,int,T);
  ARRAY_zalloc(conv->index_ratio,int,T);
  ARRAY_zalloc(conv->index_angle,int,T);
  
  conv->next_ac = NULL;
  conv->next_dc = NULL;
}

BOOL CONVCSC_is_equal(ConvCSC* conv, ConvCSC* other) {
  return conv == other;
}

BOOL CONVCSC_is_inverter(ConvCSC* conv) {
  if (conv)
    return conv->type == CONVCSC_TYPE_INV;
  else
    return FALSE;
}

BOOL CONVCSC_is_rectifier(ConvCSC* conv) {
  if (conv)
    return conv->type == CONVCSC_TYPE_REC;
  else
    return FALSE;
}

BOOL CONVCSC_is_in_P_dc_mode(ConvCSC* conv) {
  if (conv)
    return conv->mode_dc == CONVCSC_MODE_DC_CP;
  else
    return FALSE;
}

BOOL CONVCSC_is_in_i_dc_mode(ConvCSC* conv) {
  if (conv)
    return conv->mode_dc == CONVCSC_MODE_DC_CC;
  else
    return FALSE;
}

BOOL CONVCSC_is_in_v_dc_mode(ConvCSC* conv) {
  if (conv)
    return conv->mode_dc == CONVCSC_MODE_DC_CV;
  else
    return FALSE;
}

BOOL CONVCSC_has_flags(void* vconv, char flag_type, unsigned char mask) {
  ConvCSC* conv = (ConvCSC*)vconv;
  if (conv) {
    if (flag_type == FLAG_VARS)
      return (conv->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (conv->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (conv->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (conv->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL CONVCSC_has_properties(void* vconv, char prop) {
  ConvCSC* conv = (ConvCSC*)vconv;
  if (!conv)
    return FALSE;
  return TRUE;
}

ConvCSC* CONVCSC_list_ac_add(ConvCSC* conv_list, ConvCSC* conv) {
  LIST_add(ConvCSC,conv_list,conv,next_ac);
  return conv_list;
}

ConvCSC* CONVCSC_list_ac_del(ConvCSC* conv_list, ConvCSC* conv) {
  LIST_del(ConvCSC,conv_list,conv,next_ac);
  return conv_list;
}

int CONVCSC_list_ac_len(ConvCSC* conv_list) {
  int len;
  LIST_len(ConvCSC,conv_list,next_ac,len);
  return len;
}

ConvCSC* CONVCSC_list_dc_add(ConvCSC* conv_list, ConvCSC* conv) {
  LIST_add(ConvCSC,conv_list,conv,next_dc);
  return conv_list;
}

ConvCSC* CONVCSC_list_dc_del(ConvCSC* conv_list, ConvCSC* conv) {
  LIST_del(ConvCSC,conv_list,conv,next_dc);
  return conv_list;
}

int CONVCSC_list_dc_len(ConvCSC* conv_list) {
  int len;
  LIST_len(ConvCSC,conv_list,next_dc,len);
  return len;
}

ConvCSC* CONVCSC_new(int num_periods) {
  if (num_periods > 0) {
    ConvCSC* conv = (ConvCSC*)malloc(sizeof(ConvCSC));
    CONVCSC_init(conv,num_periods);
    return conv;
  }
  else
    return NULL;
}

void CONVCSC_propagate_data_in_time(ConvCSC* conv, int start, int end) {
  int t;
  if (conv) {
    if (start < 0)
      start = 0;
    if (end > conv->num_periods)
      end = conv->num_periods;
    for (t = start+1; t < end; t++) {
      conv->P_dc_set[t] = conv->P_dc_set[start];
      conv->i_dc_set[t] = conv->i_dc_set[start];
      conv->v_dc_set[t] = conv->v_dc_set[start];
      conv->P[t] = conv->P[start];
      conv->Q[t] = conv->Q[start];
      conv->P_dc[t] = conv->P_dc[start];
      conv->ratio[t] = conv->ratio[start];
      conv->angle[t] = conv->angle[start];
    }
  }
}

void CONVCSC_set_ac_bus(ConvCSC* conv, Bus* bus) {
  Bus* old_bus;
  if (conv) {
    old_bus = conv->ac_bus;
    conv->ac_bus = NULL;
    BUS_del_csc_conv(old_bus,conv);
    conv->ac_bus = bus;
    BUS_add_csc_conv(conv->ac_bus,conv);
  }
}

void CONVCSC_set_angle(ConvCSC* conv, REAL angle, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->angle[t] = angle;
}

void CONVCSC_set_angle_max(ConvCSC* conv, REAL angle_max) {
  if (conv)
    conv->angle_max = angle_max;
}

void CONVCSC_set_angle_min(ConvCSC* conv, REAL angle_min) {
  if (conv)
    conv->angle_min = angle_min;
}

void CONVCSC_set_dc_bus(ConvCSC* conv, BusDC* bus) {
  BusDC* old_bus;
  if (conv) {
    old_bus = conv->dc_bus;
    conv->dc_bus = NULL;
    BUSDC_del_csc_conv(old_bus,conv);
    conv->dc_bus = bus;
    BUSDC_add_csc_conv(conv->dc_bus,conv);
  }
}

int CONVCSC_set_flags(void* vconv, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  ConvCSC* conv = (ConvCSC*)vconv;
  int t;

  // Check conv
  if (!conv)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(conv->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(conv->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(conv->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(conv->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & CONVCSC_VAR_P) && (mask & CONVCSC_VAR_P)) { // active power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_P[t] = index+t;
    }
    (*flags_ptr) |= CONVCSC_VAR_P;
    index += conv->num_periods;
  }
  if (!((*flags_ptr) & CONVCSC_VAR_Q) && (mask & CONVCSC_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_Q[t] = index+t;
    }
    (*flags_ptr) |= CONVCSC_VAR_Q;
    index += conv->num_periods;
  }
  if (!((*flags_ptr) & CONVCSC_VAR_PDC) && (mask & CONVCSC_VAR_PDC)) { // DC power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_P_dc[t] = index+t;
      for (t = 0; t < conv->num_periods; t++)
        conv->index_i_dc[t] = index+conv->num_periods+t;
    }
    (*flags_ptr) |= CONVCSC_VAR_PDC;
    index += 2*conv->num_periods;
  }
  if (!((*flags_ptr) & CONVCSC_VAR_RATIO) && (mask & CONVCSC_VAR_RATIO)) { // tap ratio
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_ratio[t] = index+t;
    }
    (*flags_ptr) |= CONVCSC_VAR_RATIO;
    index += conv->num_periods;
  }
  if (!((*flags_ptr) & CONVCSC_VAR_ANGLE) && (mask & CONVCSC_VAR_ANGLE)) { // angle
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_angle[t] = index+t;
    }
    (*flags_ptr) |= CONVCSC_VAR_ANGLE;
    index += conv->num_periods;
  }
  return index;
}

void CONVCSC_set_i_dc_set(ConvCSC* conv, REAL i, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->i_dc_set[t] = i;
}

void CONVCSC_set_index(ConvCSC* conv, int index) { 
  if (conv)
    conv->index = index;
}

void CONVCSC_set_mode_dc(ConvCSC* conv, char mode) {
  if (conv)
    conv->mode_dc = mode;
}

void CONVCSC_set_name(ConvCSC* conv, char* name) {
  if (conv)
    strncpy(conv->name,name,(size_t)(CONVCSC_BUFFER_SIZE-1));
}

void CONVCSC_set_num_bridges(ConvCSC* conv, int num) {
  if (conv)
    conv->num_bridges = num;
}

void CONVCSC_set_P(ConvCSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P[t] = P;
}

void CONVCSC_set_P_dc(ConvCSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P_dc[t] = P;
}

void CONVCSC_set_P_dc_set(ConvCSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P_dc_set[t] = P;
}

void CONVCSC_set_Q(ConvCSC* conv, REAL Q, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->Q[t] = Q;
}

void CONVCSC_set_ratio(ConvCSC* conv, REAL ratio, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->ratio[t] = ratio;
}

void CONVCSC_set_ratio_max(ConvCSC* conv, REAL ratio_max) {
  if (conv)
    conv->ratio_max = ratio_max;
}

void CONVCSC_set_ratio_min(ConvCSC* conv, REAL ratio_min) {
  if (conv)
    conv->ratio_min = ratio_min;
}

void CONVCSC_set_type(ConvCSC* conv, char type) {
  if (conv)
    conv->type = type;
}

void CONVCSC_set_v_base_p(ConvCSC* conv, REAL v_base_p) {
  if (conv)
    conv->v_base_p = v_base_p;
}

void CONVCSC_set_v_base_s(ConvCSC* conv, REAL v_base_s) {
  if (conv)
    conv->v_base_s = v_base_s;
}

void CONVCSC_set_v_dc_set(ConvCSC* conv, REAL v, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->v_dc_set[t] = v;
}

void CONVCSC_set_var_values(ConvCSC* conv, Vec* values) {

  // Local vars
  int t;

  // No conv
  if (!conv)
    return;

  // Time loop
  for (t = 0; t < conv->num_periods; t++) {

    if (conv->vars & CONVCSC_VAR_P) // active power (p.u.)
      conv->P[t] = VEC_get(values,conv->index_P[t]);

    if (conv->vars & CONVCSC_VAR_Q) // reactive power (p.u.)
      conv->Q[t] = VEC_get(values,conv->index_Q[t]);

    if (conv->vars & CONVCSC_VAR_PDC) // DC power (p.u.)
      conv->P_dc[t] = VEC_get(values,conv->index_P_dc[t]);

    if (conv->vars & CONVCSC_VAR_RATIO) // tap ratio (p.u.)
      conv->ratio[t] = VEC_get(values,conv->index_ratio[t]);

    if (conv->vars & CONVCSC_VAR_ANGLE) // angle (rad)
      conv->angle[t] = VEC_get(values,conv->index_angle[t]);
  }
}

void CONVCSC_set_x_cap(ConvCSC* conv, REAL x_cap) {
  if (conv)
    conv->x_cap = x_cap;
}

void CONVCSC_set_x(ConvCSC* conv, REAL x) {
  if (conv)
    conv->x = x;
}

void CONVCSC_set_r(ConvCSC* conv, REAL r) {
  if (conv)
    conv->r = r;
}
