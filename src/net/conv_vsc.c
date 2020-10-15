/** @file conv_vsc.c
 *  @brief This file defines the ConvVSC data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/conv_vsc.h>
#include <pfnet/bus.h>
#include <pfnet/bus_dc.h>
#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct ConvVSC {
  
  // Properties
  char name[CONVVSC_BUFFER_SIZE]; /**< @brief Converter name */
  
  // Times
  int num_periods;     /**< @brief Number of time periods. */
  
  // Buses
  Bus* ac_bus;         /**< @brief AC bus to which the converter is connected */
  BusDC* dc_bus;       /**< @brief DC bus to which the converter is connected */
  Bus* reg_bus;        /**< @brief AC bus regulated by this converters */

  // Flags
  short int pre_cont_status;   /**< @brief Flag for indicating whether the converter was in service before applying the contingency */
  BOOL in_service;     /**< @brief Flag for indicating whether converter is in service */
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */
  
  // Control
  char mode_ac;        /**< @brief AC control mode */
  char mode_dc;        /**< @brief DC control mode */
  REAL* P_dc_set;      /**< @brief DC power set point (p.u. system base power) */
  REAL* v_dc_set;      /**< @brief DC voltage set point (p.u.) */
  
  // AC injections
  REAL* P;             /**< @brief Active power injection into AC bus (p.u. system base power) */
  REAL* Q;             /**< @brief Reactive power injection into AC bus (p.u. system base power) */

  // Loss coefficients
  REAL loss_coeff_A;   /**< @brief Loss coefficient for constant term (p.u. system base power) */
  REAL loss_coeff_B;   /**< @brief Loss coefficient for linear term (p.u. system base power / p.u. base DC current) */

  // Limits and participation
  REAL P_max;          /**< @brief Maximum active power injection into AC bus (p.u.) */
  REAL P_min;          /**< @brief Minimum active power injection into AC bus (p.u.) */
  REAL Q_max;          /**< @brief Maximum reactive power injection into AC bus (p.u.) */
  REAL Q_min;          /**< @brief Minimum reactive power injection into AC bus (p.u.) */
  REAL Q_par;          /**< @brief Converter reactive power participation factor (unitless) */
  REAL rmpct;          /**< @brief Plant Converter reactive power participation factor (Percent) */

  // Power factor
  REAL target_power_factor; /**< @brief Target power factor */

  // DC injections
  REAL* P_dc;          /**< @brief Power injection into DC bus (p.u. system base power) */
    
  // Indices
  int index;           /**< @brief Converter index */
  int* index_P;        /**< @brief Active power injection index */
  int* index_Q;        /**< @brief Reactive power injection index */
  int* index_P_dc;     /**< @brief DC power injection index */
  int* index_i_dc;     /**< @brief DC current injection index */

  // Network
  Net* net; /**< @brief Network. */
  
  // List
  ConvVSC* next_ac;       /**< @brief List of converters connected to an AC bus */
  ConvVSC* next_dc;       /**< @brief List of converters connected to a DC bus */
  ConvVSC* reg_next;      /**< @brief List of converters regulating a bus */
};

void CONVVSC_array_del(ConvVSC* conv_array, int size) {
  int i;
  ConvVSC* conv;
  if (conv_array) {
    for (i = 0; i < size; i++) {
      conv = &(conv_array[i]);
      free(conv->P_dc_set);
      free(conv->v_dc_set);
      free(conv->P);
      free(conv->Q);
      free(conv->P_dc);
      free(conv->index_P);
      free(conv->index_Q);
      free(conv->index_P_dc);
      free(conv->index_i_dc);
      CONVVSC_set_ac_bus(conv,NULL);
      CONVVSC_set_dc_bus(conv,NULL);
      CONVVSC_set_reg_bus(conv,NULL);
    }
    free(conv_array);
  }  
}

void* CONVVSC_array_get(void* conv_array, int index) { 
  if (conv_array) 
    return (void*)&(((ConvVSC*)conv_array)[index]);
  else
    return NULL;
}

ConvVSC* CONVVSC_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    ConvVSC* conv_array = (ConvVSC*)malloc(sizeof(ConvVSC)*size);
    for (i = 0; i < size; i++) {
      CONVVSC_init(&(conv_array[i]),num_periods);
      CONVVSC_set_index(&(conv_array[i]),i);
      snprintf(conv_array[i].name,(size_t)(CONVVSC_BUFFER_SIZE-1),"%d",i);
    }
    return conv_array;
  }
  else
    return NULL;
}

void CONVVSC_clear_flags(ConvVSC* conv, char flag_mask) {
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

void CONVVSC_clear_sensitivities(ConvVSC* conv) {
  // Nothing
}

void CONVVSC_copy_from_conv(ConvVSC* conv, ConvVSC* other) {

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

  // Reg bus
  // skip reg bus

  // Times
  // skip num periods

  // Properties
  strcpy(conv->name,other->name);

  // Flags
  conv->pre_cont_status = other->pre_cont_status;
  conv->in_service = other->in_service;
  conv->fixed = other->fixed;
  conv->bounded = other->bounded;
  conv->sparse = other->sparse;
  conv->vars = other->vars;

  // Control
  conv->mode_ac = other->mode_ac;
  conv->mode_dc = other->mode_dc;
  memcpy(conv->P_dc_set,other->P_dc_set,num_periods*sizeof(REAL));
  memcpy(conv->v_dc_set,other->v_dc_set,num_periods*sizeof(REAL));
  
  // AC injections
  memcpy(conv->P,other->P,num_periods*sizeof(REAL));
  memcpy(conv->Q,other->Q,num_periods*sizeof(REAL));

  // Loss coefficents
  conv->loss_coeff_A = other->loss_coeff_A;
  conv->loss_coeff_B = other->loss_coeff_B;

  // P, Q limits and participation
  conv->P_max = other->P_max;
  conv->P_min = other->P_min;
  conv->Q_max = other->Q_max;
  conv->Q_min = other->Q_min;
  conv->Q_par = other->Q_par;
  conv->rmpct = other->rmpct;

  // Power factor
  conv->target_power_factor = other->target_power_factor;

  // DC injection
  memcpy(conv->P_dc,other->P_dc,num_periods*sizeof(REAL));

  // Indices
  // skip index
  memcpy(conv->index_P,other->index_P,num_periods*sizeof(int));
  memcpy(conv->index_Q,other->index_Q,num_periods*sizeof(int));
  memcpy(conv->index_P_dc,other->index_P_dc,num_periods*sizeof(int));
  memcpy(conv->index_i_dc,other->index_i_dc,num_periods*sizeof(int));
    
  // List
  // skip next
}

short int CONVVSC_get_pre_cont_status(void* conv) {
  if (conv)
    return ((ConvVSC*)conv)->pre_cont_status;
  else
    return 0;
}

Bus* CONVVSC_get_ac_bus(ConvVSC* conv) {
  if (conv)
    return conv->ac_bus;
  else
    return NULL;
}

BusDC* CONVVSC_get_dc_bus(ConvVSC* conv) {
  if (conv)
    return conv->dc_bus;
  else
    return NULL;
}

Bus* CONVVSC_get_reg_bus(ConvVSC* conv) {
  if (conv)
    return conv->reg_bus;
  else
    return NULL;
}

char CONVVSC_get_flags_vars(ConvVSC* conv) {
  if (conv)
    return conv->vars;
  else
    return 0;
}

char CONVVSC_get_flags_fixed(ConvVSC* conv) {
  if (conv)
    return conv->fixed;
  else
    return 0;
}

char CONVVSC_get_flags_bounded(ConvVSC* conv) {
  if (conv)
    return conv->bounded;
  else
    return 0;
}

char CONVVSC_get_flags_sparse(ConvVSC* conv) {
  if (conv)
    return conv->sparse;
  else
    return 0;
}

REAL CONVVSC_get_loss_coeff_A(ConvVSC* conv) {
  if (conv)
    return conv->loss_coeff_A;
  else
    return 0;
}

REAL CONVVSC_get_loss_coeff_B(ConvVSC* conv) {
  if (conv)
    return conv->loss_coeff_B;
  else
    return 0;
}

REAL CONVVSC_get_v_dc_set(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->v_dc_set[t];
  else 
    return 1.;
}

REAL CONVVSC_get_P_dc_set(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc_set[t];
  else 
    return 0;
}

int CONVVSC_get_index(ConvVSC* conv) {
  if (conv)
    return conv->index;
  else
    return -1;
}

int CONVVSC_get_index_P(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_P[t];
  else
    return -1;
}

int CONVVSC_get_index_Q(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_Q[t];
  else
    return -1;
}

int CONVVSC_get_index_P_dc(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_P_dc[t];
  else
    return -1;
}

int CONVVSC_get_index_i_dc(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->index_i_dc[t];
  else
    return -1;
}

int* CONVVSC_get_index_P_array(ConvVSC* conv) {
  if (conv)
    return conv->index_P;
  else
    return NULL;
}

int* CONVVSC_get_index_Q_array(ConvVSC* conv) {
  if (conv)
    return conv->index_Q;
  else
    return NULL;
}

int* CONVVSC_get_index_P_dc_array(ConvVSC* conv) {
  if (conv)
    return conv->index_P_dc;
  else
    return NULL;
}

int* CONVVSC_get_index_i_dc_array(ConvVSC* conv) {
  if (conv)
    return conv->index_i_dc;
  else
    return NULL;
}

REAL CONVVSC_get_target_power_factor(ConvVSC* conv) {
  if (conv)
    return conv->target_power_factor;
  else
    return 1.;
}

char* CONVVSC_get_json_string(ConvVSC* conv, char* output) {

  // Local variables
  char temp[CONVVSC_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No conv
  if (!conv)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*CONVVSC_BUFFER_SIZE*CONVVSC_NUM_JSON_FIELDS*conv->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",conv->index,FALSE);
  JSON_str(temp,output,"name",conv->name,FALSE);
  JSON_bool(temp,output,"in_service",conv->in_service,FALSE);
  JSON_int(temp,output,"num_periods",conv->num_periods,FALSE);
  JSON_obj(temp,output,"ac_bus",conv->ac_bus,BUS_get_index,FALSE);
  JSON_obj(temp,output,"dc_bus",conv->dc_bus,BUSDC_get_index,FALSE);
  JSON_obj(temp,output,"reg_bus",conv->reg_bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"mode_ac",conv->mode_ac,FALSE);
  JSON_int(temp,output,"mode_dc",conv->mode_dc,FALSE);
  JSON_array_float(temp,output,"P_dc_set",conv->P_dc_set,conv->num_periods,FALSE);
  JSON_array_float(temp,output,"v_dc_set",conv->v_dc_set,conv->num_periods,FALSE);
  JSON_array_float(temp,output,"P",conv->P,conv->num_periods,FALSE);
  JSON_array_float(temp,output,"Q",conv->Q,conv->num_periods,FALSE);
  JSON_array_float(temp,output,"P_dc",conv->P_dc,conv->num_periods,FALSE);
  JSON_float(temp,output,"loss_coeff_A",conv->loss_coeff_A,FALSE);
  JSON_float(temp,output,"loss_coeff_B",conv->loss_coeff_B,FALSE);
  JSON_float(temp,output,"P_max",conv->P_max,FALSE);
  JSON_float(temp,output,"P_min",conv->P_min,FALSE);
  JSON_float(temp,output,"Q_max",conv->Q_max,FALSE);
  JSON_float(temp,output,"Q_min",conv->Q_min,FALSE);
  JSON_float(temp,output,"Q_par",conv->Q_par,FALSE);
  JSON_float(temp,output,"rmpct",conv->rmpct,FALSE);
  JSON_float(temp,output,"target_power_factor",conv->target_power_factor,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

char CONVVSC_get_mode_ac(ConvVSC* conv) {
  if (conv)
    return conv->mode_ac;
  else
    return CONVVSC_MODE_AC_NC;
}

char CONVVSC_get_mode_dc(ConvVSC* conv) {
  if (conv)
    return conv->mode_dc;
  else
    return CONVVSC_MODE_DC_NC;
}

char* CONVVSC_get_name(ConvVSC* conv) {
  if (conv)
    return conv->name;
  else
    return NULL;
}

ConvVSC* CONVVSC_get_next_ac(ConvVSC* conv) {
  if (conv)
    return conv->next_ac;
  else
    return NULL;
}

ConvVSC* CONVVSC_get_next_dc(ConvVSC* conv) {
  if (conv)
    return conv->next_dc;
  else
    return NULL;
}

ConvVSC* CONVVSC_get_reg_next(ConvVSC* conv) {
  if (conv)
    return conv->reg_next;
  else
    return NULL;
}

int CONVVSC_get_num_periods(ConvVSC* conv) {
  if (conv)
    return conv->num_periods;
  else
    return 0;
}

int CONVVSC_get_num_vars(void* vconv, unsigned char var, int t_start, int t_end) {

  // Local vars
  ConvVSC* conv = (ConvVSC*)vconv;
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
  if ((var & CONVVSC_VAR_P) && (conv->vars & CONVVSC_VAR_P)) // active power
    num_vars += dt;
  if ((var & CONVVSC_VAR_Q) && (conv->vars & CONVVSC_VAR_Q)) // reactive power
    num_vars += dt;
  if ((var & CONVVSC_VAR_PDC) && (conv->vars & CONVVSC_VAR_PDC)) // DC power
    num_vars += 2*dt;
  return num_vars;
}

char CONVVSC_get_obj_type(void* conv) {
  if (conv)
    return OBJ_CONVVSC;
  else
    return OBJ_UNKNOWN;
}

REAL CONVVSC_get_P(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P[t];
  else 
    return 0;
}

REAL CONVVSC_get_Q(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->Q[t];
  else 
    return 0;
}

REAL CONVVSC_get_P_dc(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc[t];
  else 
    return 0;
}

REAL CONVVSC_get_i_dc(ConvVSC* conv, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    return conv->P_dc[t]/(BUSDC_get_v(conv->dc_bus,t));
  else 
    return 0;
}

REAL CONVVSC_get_P_max(ConvVSC* conv) {
  if (conv)
    return conv->P_max;
  else
    return 0;
}

REAL CONVVSC_get_P_min(ConvVSC* conv) {
  if (conv)
    return conv->P_min;
  else
    return 0;
}

REAL CONVVSC_get_Q_max(ConvVSC* conv) {
  if (conv)
    return conv->Q_max;
  else
    return 0;
}

REAL CONVVSC_get_Q_min(ConvVSC* conv) {
  if (conv)
    return conv->Q_min;
  else
    return 0;
}

REAL CONVVSC_get_Q_par(ConvVSC* conv) {
  if (conv)
    return conv->Q_par;
  else
    return 0;
}

REAL CONVVSC_get_rmpct(ConvVSC* conv) {
  if (conv)
    return conv->rmpct;
  else
    return 0;
}

REAL* CONVVSC_get_P_array(ConvVSC* conv) {
  if (conv)
    return conv->P;
  else
    return NULL;
}

REAL* CONVVSC_get_P_dc_array(ConvVSC* conv) {
  if (conv)
    return conv->P_dc;
  else
    return NULL;
}

REAL* CONVVSC_get_Q_array(ConvVSC* conv) {
  if (conv)
    return conv->Q;
  else
    return NULL;
}

REAL* CONVVSC_get_v_dc_set_array(ConvVSC* conv) {
  if (conv)
    return conv->v_dc_set;
  else
    return NULL;
}

REAL* CONVVSC_get_P_dc_set_array(ConvVSC* conv) {
  if (conv)
    return conv->P_dc_set;
  else
    return NULL;
}

Vec* CONVVSC_get_var_indices(void* vconv, unsigned char var, int t_start, int t_end) {

  // Local vars
  ConvVSC* conv = (ConvVSC*)vconv;
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
  indices = VEC_new(CONVVSC_get_num_vars(vconv,var,t_start,t_end));
  if ((var & CONVVSC_VAR_P) && (conv->vars & CONVVSC_VAR_P)) { // active power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_P[t]);
      offset++;
    }
  }
  if ((var & CONVVSC_VAR_Q) && (conv->vars & CONVVSC_VAR_Q)) { // reactive power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_Q[t]);
      offset++;
    }
  }
  if ((var & CONVVSC_VAR_PDC) && (conv->vars & CONVVSC_VAR_PDC)) { // DC power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,conv->index_P_dc[t]);
      VEC_set(indices,offset+1,conv->index_i_dc[t]);
      offset += 2;
    }
  }
  
  return indices;
}

char* CONVVSC_get_var_info_string(ConvVSC* conv, int index) {

  // Local variables
  char* info;
  
  //Check
  if (!conv)
    return NULL;

  // Active power
  if ((conv->vars & CONVVSC_VAR_P) &&
      index >= conv->index_P[0] &&
      index <= conv->index_P[conv->num_periods-1]) {
    info = (char*)malloc(CONVVSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVVSC_BUFFER_SIZE*sizeof(char),
             "vsc_conv:%d:active power:%d",conv->index,index-conv->index_P[0]);
    return info;
  }

  // Reactive power
  if ((conv->vars & CONVVSC_VAR_Q) &&
      index >= conv->index_Q[0] &&
      index <= conv->index_Q[conv->num_periods-1]) {
    info = (char*)malloc(CONVVSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVVSC_BUFFER_SIZE*sizeof(char),
             "vsc_conv:%d:reactive power:%d",conv->index,index-conv->index_Q[0]);
    return info;
  }

  // DC power
  if ((conv->vars & CONVVSC_VAR_PDC) &&
      index >= conv->index_P_dc[0] &&
      index <= conv->index_P_dc[conv->num_periods-1]) {
    info = (char*)malloc(CONVVSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVVSC_BUFFER_SIZE*sizeof(char),
             "vsc_conv:%d:dc power:%d",conv->index,index-conv->index_P_dc[0]);
    return info;
  }

  // DC current
  if ((conv->vars & CONVVSC_VAR_PDC) &&
      index >= conv->index_i_dc[0] &&
      index <= conv->index_i_dc[conv->num_periods-1]) {
    info = (char*)malloc(CONVVSC_BUFFER_SIZE*sizeof(char));
    snprintf(info,CONVVSC_BUFFER_SIZE*sizeof(char),
             "vsc_conv:%d:dc current:%d",conv->index,index-conv->index_i_dc[0]);
    return info;
  }

  // Return
  return NULL;
}

void CONVVSC_get_var_values(ConvVSC* conv, Vec* values, int code) {
 
  // Local vars
  int t;
 
  // No laod
  if (!conv)
    return;

  // Time loop
  for (t = 0; t < conv->num_periods; t++) {

    if (conv->vars & CONVVSC_VAR_P) { // active power
      switch(code) {
      
      case UPPER_LIMITS:
        if (conv->bounded & CONVVSC_VAR_P)
          VEC_set(values,conv->index_P[t],conv->P_max);
        else
          VEC_set(values,conv->index_P[t],CONVVSC_INF_P);
        break;
      
      case LOWER_LIMITS:
        if (conv->bounded & CONVVSC_VAR_P)
          VEC_set(values,conv->index_P[t],conv->P_min);
        else
          VEC_set(values,conv->index_P[t],-CONVVSC_INF_P);
        break;
      
      default:
        VEC_set(values,conv->index_P[t],conv->P[t]);
      }
    }

    if (conv->vars & CONVVSC_VAR_Q) { // reactive power
      switch(code) {
      
      case UPPER_LIMITS:
        if (conv->bounded & CONVVSC_VAR_Q)
          VEC_set(values,conv->index_Q[t],conv->Q_max);
        else
          VEC_set(values,conv->index_Q[t],CONVVSC_INF_Q);
        break;
      
      case LOWER_LIMITS:
        if (conv->bounded & CONVVSC_VAR_Q)
          VEC_set(values,conv->index_Q[t],conv->Q_min);
        else
          VEC_set(values,conv->index_Q[t],-CONVVSC_INF_Q);
        break;
      
      default:
        VEC_set(values,conv->index_Q[t],conv->Q[t]);
      }
    }

    if (conv->vars & CONVVSC_VAR_PDC) { // DC power
      switch(code) {
        
      case UPPER_LIMITS:
        VEC_set(values,conv->index_P_dc[t],CONVVSC_INF_PDC);
        VEC_set(values,conv->index_i_dc[t],CONVVSC_INF_PDC);
        break;
        
      case LOWER_LIMITS:
        VEC_set(values,conv->index_P_dc[t],-CONVVSC_INF_PDC);
        VEC_set(values,conv->index_i_dc[t],-CONVVSC_INF_PDC);
        break;
        
      default:
        VEC_set(values,conv->index_P_dc[t],conv->P_dc[t]);
        VEC_set(values,conv->index_i_dc[t],conv->P_dc[t]/(BUSDC_get_v(conv->dc_bus,t)));
      }
    }
  }
}


void CONVVSC_init(ConvVSC* conv, int num_periods) {

  // Local variables
  int T;
  
  // No conv
  if (!conv)
    return;

  ARRAY_clear(conv->name,char,CONVVSC_BUFFER_SIZE);
  
  T = num_periods;
  conv->num_periods = num_periods;

  conv->ac_bus = NULL;
  conv->dc_bus = NULL;
  conv->reg_bus = NULL;

  conv->pre_cont_status = PRE_CONT_UNSET;
  conv->in_service = TRUE;
  conv->fixed = 0x00;
  conv->bounded = 0x00;
  conv->sparse = 0x00;
  conv->vars = 0x00;

  conv->mode_ac = CONVVSC_MODE_AC_NC;
  conv->mode_dc = CONVVSC_MODE_DC_NC;
  ARRAY_zalloc(conv->P_dc_set,REAL,T);
  ARRAY_zalloc(conv->v_dc_set,REAL,T);

  ARRAY_zalloc(conv->P,REAL,T);
  ARRAY_zalloc(conv->Q,REAL,T);

  ARRAY_zalloc(conv->P_dc,REAL,T);

  conv->P_max = 0;
  conv->P_min = 0;
  conv->Q_max = 0;
  conv->Q_min = 0;
  conv->Q_par = 1.;
  conv->rmpct = 100.;

  conv->target_power_factor = 1.;
  
  conv->index = -1;
  ARRAY_zalloc(conv->index_P,int,T);
  ARRAY_zalloc(conv->index_Q,int,T);
  ARRAY_zalloc(conv->index_P_dc,int,T);
  ARRAY_zalloc(conv->index_i_dc,int,T);

  conv->loss_coeff_A = 0;
  conv->loss_coeff_B = 0;

  conv->net = NULL;
  
  conv->next_ac = NULL;
  conv->next_dc = NULL;
  conv->reg_next = NULL;
}

BOOL CONVVSC_is_in_service(void* conv) {
  if (conv)
    return (((ConvVSC*)conv)->in_service &&
            BUS_is_in_service(((ConvVSC*)conv)->ac_bus) &&
            BUSDC_is_in_service(((ConvVSC*)conv)->dc_bus));
  else
    return FALSE;
}

BOOL CONVVSC_is_equal(ConvVSC* conv, ConvVSC* other) {
  return conv == other;
}

BOOL CONVVSC_is_in_f_ac_mode(ConvVSC* conv) {
  if (conv)
    return conv->mode_ac == CONVVSC_MODE_AC_CF;
  else
    return FALSE;
}

BOOL CONVVSC_is_in_v_ac_mode(ConvVSC* conv) {
  if (conv)
    return conv->mode_ac == CONVVSC_MODE_AC_CV;
  else
    return FALSE;
}

BOOL CONVVSC_is_in_P_dc_mode(ConvVSC* conv) {
  if (conv)
    return conv->mode_dc == CONVVSC_MODE_DC_CP;
  else
    return FALSE;
}

BOOL CONVVSC_is_in_v_dc_mode(ConvVSC* conv) {
  if (conv)
    return conv->mode_dc == CONVVSC_MODE_DC_CV;
  else
    return FALSE;
}

BOOL CONVVSC_has_flags(void* vconv, char flag_type, unsigned char mask) {
  ConvVSC* conv = (ConvVSC*)vconv;
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

BOOL CONVVSC_has_properties(void* vconv, char prop) {
  ConvVSC* conv = (ConvVSC*)vconv;
  if (!conv)
    return FALSE;
  return TRUE;
}

ConvVSC* CONVVSC_list_ac_add(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_add(ConvVSC,conv_list,conv,next_ac);
  return conv_list;
}

ConvVSC* CONVVSC_list_ac_del(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_del(ConvVSC,conv_list,conv,next_ac);
  return conv_list;
}

int CONVVSC_list_ac_len(ConvVSC* conv_list) {
  int len;
  LIST_len(ConvVSC,conv_list,next_ac,len);
  return len;
}

ConvVSC* CONVVSC_list_dc_add(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_add(ConvVSC,conv_list,conv,next_dc);
  return conv_list;
}

ConvVSC* CONVVSC_list_dc_del(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_del(ConvVSC,conv_list,conv,next_dc);
  return conv_list;
}

int CONVVSC_list_dc_len(ConvVSC* conv_list) {
  int len;
  LIST_len(ConvVSC,conv_list,next_dc,len);
  return len;
}

ConvVSC* CONVVSC_list_reg_add(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_add(ConvVSC,conv_list,conv,reg_next);
  return conv_list;
}

ConvVSC* CONVVSC_list_reg_del(ConvVSC* conv_list, ConvVSC* conv) {
  LIST_del(ConvVSC,conv_list,conv,reg_next);
  return conv_list;
}

int CONVVSC_list_reg_len(ConvVSC* conv_list) {
  int len;
  LIST_len(ConvVSC,conv_list,reg_next,len);
  return len;
}

ConvVSC* CONVVSC_new(int num_periods) {
  if (num_periods > 0) {
    ConvVSC* conv = (ConvVSC*)malloc(sizeof(ConvVSC));
    CONVVSC_init(conv,num_periods);
    return conv;
  }
  else
    return NULL;
}

void CONVVSC_propagate_data_in_time(ConvVSC* conv, int start, int end) {
  int t;
  if (conv) {
    if (start < 0)
      start = 0;
    if (end > conv->num_periods)
      end = conv->num_periods;
    for (t = start+1; t < end; t++) {
      conv->P_dc_set[t] = conv->P_dc_set[start];
      conv->v_dc_set[t] = conv->v_dc_set[start];
      conv->P[t] = conv->P[start];
      conv->Q[t] = conv->Q[start];
      conv->P_dc[t] = conv->P_dc[start];
    }
  }
}

void CONVVSC_set_pre_cont_status(ConvVSC* conv, short int pre_cont_status) {
  if (conv && BUS_is_in_service(conv->ac_bus) && BUSDC_is_in_service(conv->dc_bus))
    conv->pre_cont_status = pre_cont_status;
}

void CONVVSC_set_in_service(ConvVSC* conv, BOOL in_service) {
  if (conv && BUS_is_in_service(conv->ac_bus) && BUSDC_is_in_service(conv->dc_bus)) {
    if (conv->in_service != in_service)
      NET_inc_state_tag(conv->net);
    conv->in_service = in_service;
  }
}

void CONVVSC_set_network(ConvVSC* conv, void* net) {
  if (conv)
    conv->net = (Net*)net;
}

void CONVVSC_set_P_max(ConvVSC* conv, REAL P_max) {
  if (conv)
    conv->P_max = P_max;
}

void CONVVSC_set_P_min(ConvVSC* conv, REAL P_min) {
  if (conv)
    conv->P_min = P_min;
}

void CONVVSC_set_Q_max(ConvVSC* conv, REAL Q_max) {
  if (conv)
    conv->Q_max = Q_max;
}

void CONVVSC_set_Q_min(ConvVSC* conv, REAL Q_min) {
  if (conv)
    conv->Q_min = Q_min;
}

void CONVVSC_set_Q_par(ConvVSC* conv, REAL Q_par) {
  if (conv)
    conv->Q_par = Q_par;
}

void CONVVSC_set_rmpct(ConvVSC* conv, REAL rmpct) {

  // Local variables
  Bus* bus;
  ConvVSC* vsc;
  if (conv) {
    bus = CONVVSC_get_ac_bus(conv);
    // Change for all converters at the same bus
    for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) { 
      vsc->rmpct = rmpct;
    }
  }
}

void CONVVSC_set_ac_bus(ConvVSC* conv, Bus* bus) {
  Bus* old_bus;
  if (conv) {
    old_bus = conv->ac_bus;
    conv->ac_bus = NULL;
    BUS_del_vsc_conv(old_bus,conv);
    conv->ac_bus = bus;
    BUS_add_vsc_conv(conv->ac_bus,conv);
  }
}

void CONVVSC_set_dc_bus(ConvVSC* conv, BusDC* bus) {
  BusDC* old_bus;
  if (conv) {
    old_bus = conv->dc_bus;
    conv->dc_bus = NULL;
    BUSDC_del_vsc_conv(old_bus,conv);
    conv->dc_bus = bus;
    BUSDC_add_vsc_conv(conv->dc_bus,conv);
  }
}

void CONVVSC_set_reg_bus(ConvVSC* conv, Bus* reg_bus) {
  Bus* old_reg_bus;
  if (conv) {
    old_reg_bus = conv->reg_bus;
    conv->reg_bus = NULL;
    BUS_del_reg_vsc_conv(old_reg_bus,conv);
    conv->reg_bus = reg_bus;
    BUS_add_reg_vsc_conv(conv->reg_bus,conv);
  }
}

int CONVVSC_set_flags(void* vconv, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  ConvVSC* conv = (ConvVSC*)vconv;
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
  if (!((*flags_ptr) & CONVVSC_VAR_P) && (mask & CONVVSC_VAR_P)) { // active power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_P[t] = index+t;
    }
    (*flags_ptr) |= CONVVSC_VAR_P;
    index += conv->num_periods;
  }
  if (!((*flags_ptr) & CONVVSC_VAR_Q) && (mask & CONVVSC_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_Q[t] = index+t;
    }
    (*flags_ptr) |= CONVVSC_VAR_Q;
    index += conv->num_periods;
  }
  if (!((*flags_ptr) & CONVVSC_VAR_PDC) && (mask & CONVVSC_VAR_PDC)) { // DC power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < conv->num_periods; t++)
        conv->index_P_dc[t] = index+t;
      for (t = 0; t < conv->num_periods; t++)
        conv->index_i_dc[t] = index+conv->num_periods+t;
    }
    (*flags_ptr) |= CONVVSC_VAR_PDC;
    index += 2*conv->num_periods;
  }
  return index;
}

void CONVVSC_set_index(ConvVSC* conv, int index) { 
  if (conv)
    conv->index = index;
}

void CONVVSC_set_mode_dc(ConvVSC* conv, char mode) {
  if (conv)
    conv->mode_dc = mode;
}

void CONVVSC_set_mode_ac(ConvVSC* conv, char mode) {
  if (conv)
    conv->mode_ac = mode;
}

void CONVVSC_set_name(ConvVSC* conv, char* name) {
  if (conv)
    strncpy(conv->name,name,(size_t)(CONVVSC_BUFFER_SIZE-1));
}

void CONVVSC_set_P(ConvVSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P[t] = P;
}

void CONVVSC_set_P_dc(ConvVSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P_dc[t] = P;
}

void CONVVSC_set_Q(ConvVSC* conv, REAL Q, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->Q[t] = Q;
}

void CONVVSC_set_var_values(ConvVSC* conv, Vec* values) {

  // Local vars
  int t;

  // No conv
  if (!conv)
    return;

  // Time loop
  for (t = 0; t < conv->num_periods; t++) {

    if (conv->vars & CONVVSC_VAR_P) // active power (p.u.)
      conv->P[t] = VEC_get(values,conv->index_P[t]);

    if (conv->vars & CONVVSC_VAR_Q) // reactive power (p.u.)
      conv->Q[t] = VEC_get(values,conv->index_Q[t]);

    if (conv->vars & CONVVSC_VAR_PDC) // DC power (p.u.)
      conv->P_dc[t] = VEC_get(values,conv->index_P_dc[t]);
  }
}

void CONVVSC_set_loss_coeff_A(ConvVSC* conv, REAL A) {
  if (conv)
    conv->loss_coeff_A = A;
}

void CONVVSC_set_loss_coeff_B(ConvVSC* conv, REAL B) {
  if (conv)
    conv->loss_coeff_B = B;
}

void CONVVSC_set_v_dc_set(ConvVSC* conv, REAL v, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->v_dc_set[t] = v;
}

void CONVVSC_set_P_dc_set(ConvVSC* conv, REAL P, int t) {
  if (conv && t >= 0 && t < conv->num_periods)
    conv->P_dc_set[t] = P;
}

void CONVVSC_set_target_power_factor(ConvVSC* conv, REAL pf) {
  if (conv) {
    if (pf > 1.)
      conv->target_power_factor = 1.;
    else if (pf < -1.)
      conv->target_power_factor = -1.;
    else if (fabs(pf) < CONVVSC_MIN_TARGET_PF)
      conv->target_power_factor = CONVVSC_MIN_TARGET_PF;
    else
      conv->target_power_factor = pf;
  }
}
