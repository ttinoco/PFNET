/** @file bus_dc.c
 *  @brief This file defines the BusDC data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/bus_dc.h>
#include <pfnet/conv_csc.h>
#include <pfnet/conv_vsc.h>
#include <pfnet/branch_dc.h>
#include <pfnet/array.h>
#include <pfnet/net.h>
#include <pfnet/json_macros.h>

struct BusDC {

  // Properties
  int number;                 /**< @brief Bus number */
  char name[BUSDC_BUFFER_SIZE]; /**< @brief Bus name */
  
  // Times
  int num_periods;  /**< @brief Number of time periods */
  
  // Voltage
  REAL v_base;    /**< @brief Base voltage (kilovolts) */
  REAL* v;        /**< @brief Voltage (p.u.) */

  // Flags
  BOOL in_service; /**< @brief Flag for indicating whether the bus is in service */
  char fixed;      /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;    /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;     /**< @brief Flags for indicating which control adjustments should be sparse */
  char vars;       /**< @brief Flags for indicating which quantities should be treated as variables */

  // Indices
  int index;      /**< @brief Bus index */
  int* index_v;   /**< @brief Voltage index */

  // Componnets
  ConvCSC* csc_conv;   /**< @brief List of CSC converters connected to bus */
  ConvVSC* vsc_conv;   /**< @brief List of VSC converters connected to bus */
  BranchDC* branch_k;  /**< @brief List of DC branches having this bus on the "k" side */
  BranchDC* branch_m;  /**< @brief List of DC branches having this bus on the "m" side */

  // Mismatches
  REAL* P_mis;    /**< @brief Active power mismatch (p.u. system base power) */

  // Network
  Net* net; /**< @brief Network. */

  // Hash
  UT_hash_handle hh_number; /**< @brief Handle for bus hash table based on numbers */
  UT_hash_handle hh_name;   /**< @brief Handle for bus hash table based on names */

  // HVDC helper
  int* di_index;
};

void BUSDC_add_branch_k(BusDC* bus, BranchDC* branch) {
  if (bus) {
    bus->branch_k = BRANCHDC_list_k_add(bus->branch_k,branch);
    if (BRANCHDC_get_bus_k(branch) != bus)
      BRANCHDC_set_bus_k(branch,bus);
  }
}

void BUSDC_add_branch_m(BusDC* bus, BranchDC* branch) {
  if (bus) {
    bus->branch_m = BRANCHDC_list_m_add(bus->branch_m,branch);
    if (BRANCHDC_get_bus_m(branch) != bus)
      BRANCHDC_set_bus_m(branch,bus);
  }
}

void BUSDC_add_csc_conv(BusDC* bus, ConvCSC* conv) {
  if (bus) {
    bus->csc_conv = CONVCSC_list_dc_add(bus->csc_conv,conv);
    if (CONVCSC_get_dc_bus(conv) != bus)
      CONVCSC_set_dc_bus(conv,bus);
  }
}

void BUSDC_add_vsc_conv(BusDC* bus, ConvVSC* conv) {
  if (bus) {
    bus->vsc_conv = CONVVSC_list_dc_add(bus->vsc_conv,conv);
    if (CONVVSC_get_dc_bus(conv) != bus)
      CONVVSC_set_dc_bus(conv,bus);
  }
}

void BUSDC_array_del(BusDC* bus_array, int size) {
  int i;
  BusDC* bus;
  if (bus_array) {
    for (i = 0; i < size; i++) {
      bus = &(bus_array[i]);
      free(bus->v);
      free(bus->index_v);
      free(bus->P_mis);
      free(bus->di_index);
      BUSDC_del_all_connections(bus);
    }
    free(bus_array);
  }
}

void* BUSDC_array_get(void* bus_array, int index) {
  if (bus_array)
    return (void*)&(((BusDC*)bus_array)[index]);
  else
    return NULL;
}

void BUSDC_array_get_max_mismatches(BusDC* bus_array, int size, REAL* P, int t) {
  int i;
  if (bus_array && t >= 0 && t < bus_array->num_periods) {
    for (i = 0; i < size; i++) {
      if (fabs(bus_array[i].P_mis[t]) > *P)
        *P = fabs(bus_array[i].P_mis[t]);
    }
  }
}

BusDC* BUSDC_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    BusDC* bus_array = (BusDC*)malloc(sizeof(BusDC)*size);
    for (i = 0; i < size; i++) {
      BUSDC_init(&(bus_array[i]),num_periods);
      BUSDC_set_index(&(bus_array[i]),i);
      snprintf(bus_array[i].name,(size_t)(BUSDC_BUFFER_SIZE-1),"%d",i);
    }
    return bus_array;
  }
  else
    return NULL;
}

void BUSDC_clear_flags(BusDC* bus, char flag_mask) {
  if (bus) {
    if (flag_mask & FLAG_VARS)
      bus->vars = 0x00;
    if (flag_mask & FLAG_BOUNDED)
      bus->bounded = 0x00;
    if (flag_mask & FLAG_FIXED)
      bus->fixed = 0x00;
    if (flag_mask & FLAG_SPARSE)
      bus->sparse = 0x00;
  }
}

void BUSDC_clear_mismatches(BusDC* bus) {
  int t;
  if (bus) {
    for (t = 0; t < bus->num_periods; t++)
      bus->P_mis[t] = 0;
  }
}

void BUSDC_clear_sensitivities(BusDC* bus) {
  // Nothing
}

void BUSDC_copy_from_dc_bus(BusDC* bus, BusDC* other) {
  /** Copies data from another bus except
   *  index, hash info, and connections.
   */

  // Local variables
  int num_periods;

  // Check
  if (!bus || !other)
    return;

  // Min num periods
  if (bus->num_periods < other->num_periods)
    num_periods = bus->num_periods;
  else
    num_periods = other->num_periods;

  // Properties
  bus->number = other->number;
  strcpy(bus->name,other->name);
  
  // Time
  // skip num periods

  // Voltage
  bus->v_base = other->v_base;
  memcpy(bus->v,other->v,num_periods*sizeof(REAL));

  // Flags
  bus->in_service = other->in_service;
  bus->fixed = other->fixed;
  bus->bounded = other->bounded;
  bus->sparse = other->sparse;
  bus->vars = other->vars;

  // Indices
  // skip index
  memcpy(bus->index_v,other->index_v,num_periods*sizeof(int));

  // Connections
  // skip connections

  // Mismatches
  memcpy(bus->P_mis,other->P_mis,num_periods*sizeof(REAL));

  // Hash
  // skip hash
}

void BUSDC_del_all_connections(BusDC* bus) {
  if (bus) {
    while (bus->csc_conv)
      BUSDC_del_csc_conv(bus,bus->csc_conv);
    while (bus->vsc_conv)
      BUSDC_del_vsc_conv(bus,bus->vsc_conv);
    while (bus->branch_k)
      BUSDC_del_branch_k(bus,bus->branch_k);
    while (bus->branch_m)
      BUSDC_del_branch_m(bus,bus->branch_m);
  }
}

void BUSDC_del_branch_k(BusDC* bus, BranchDC* branch) {
  if (bus) {
    bus->branch_k = BRANCHDC_list_k_del(bus->branch_k,branch);
    if (BRANCHDC_get_bus_k(branch) == bus)
      BRANCHDC_set_bus_k(branch,NULL);
  }
}

void BUSDC_del_branch_m(BusDC* bus, BranchDC* branch) {
  if (bus) {
    bus->branch_m = BRANCHDC_list_m_del(bus->branch_m,branch);
    if (BRANCHDC_get_bus_m(branch) == bus)
      BRANCHDC_set_bus_m(branch,NULL);
  }
}

void BUSDC_del_csc_conv(BusDC* bus, ConvCSC* conv) {
  if (bus) {
    bus->csc_conv = CONVCSC_list_dc_del(bus->csc_conv,conv);
    if (CONVCSC_get_dc_bus(conv) == bus)
      CONVCSC_set_dc_bus(conv,NULL);
  }
}

void BUSDC_del_vsc_conv(BusDC* bus, ConvVSC* conv) {
  if (bus) {
    bus->vsc_conv = CONVVSC_list_dc_del(bus->vsc_conv,conv);
    if (CONVVSC_get_dc_bus(conv) == bus)
      CONVVSC_set_dc_bus(conv,NULL);
  }
}

BranchDC* BUSDC_get_branch_k(BusDC* bus) {
  if (bus)
    return bus->branch_k;
  else
    return NULL;
}

BranchDC* BUSDC_get_branch_m(BusDC* bus) {
  if (bus)
    return bus->branch_m;
  else
    return NULL;
}

ConvCSC* BUSDC_get_csc_conv(BusDC* bus) {
  if (bus)
    return bus->csc_conv;
  else
    return NULL;
}

ConvVSC* BUSDC_get_vsc_conv(BusDC* bus) {
  if (bus)
    return bus->vsc_conv;
  else
    return NULL;
}

char BUSDC_get_flags_vars(BusDC* bus) {
  if (bus)
    return bus->vars;
  else
    return 0;
}

char BUSDC_get_flags_fixed(BusDC* bus) {
  if (bus)
    return bus->fixed;
  else
    return 0;
}

char BUSDC_get_flags_bounded(BusDC* bus) {
  if (bus)
    return bus->bounded;
  else
    return 0;
}

char BUSDC_get_flags_sparse(BusDC* bus) {
  if (bus)
    return bus->sparse;
  else
    return 0;
}

int BUSDC_get_index(BusDC* bus) {
  if (bus)
    return bus->index;
  else
    return -1;
}

int BUSDC_get_index_v(BusDC* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->index_v[t];
  else
    return -1;
}

char* BUSDC_get_json_string(BusDC* bus, char* output) {
  
  // Local variables
  char temp[BUSDC_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No bus
  if (!bus)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*BUSDC_BUFFER_SIZE*BUSDC_NUM_JSON_FIELDS*bus->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",bus->index,FALSE);
  JSON_int(temp,output,"number",bus->number,FALSE);
  JSON_str(temp,output,"name",bus->name,FALSE);
  JSON_bool(temp,output,"in_service",bus->in_service,FALSE);
  JSON_int(temp,output,"num_periods",bus->num_periods,FALSE);
  JSON_float(temp,output,"v_base",bus->v_base,FALSE);
  JSON_array_float(temp,output,"v",bus->v,bus->num_periods,FALSE);
  JSON_list_int(temp,output,"branches_k",bus,BranchDC,BUSDC_get_branch_k,BRANCHDC_get_index,BRANCHDC_get_next_k,FALSE);
  JSON_list_int(temp,output,"branches_m",bus,BranchDC,BUSDC_get_branch_m,BRANCHDC_get_index,BRANCHDC_get_next_m,FALSE);
  JSON_list_int(temp,output,"csc_converters",bus,ConvCSC,BUSDC_get_csc_conv,CONVCSC_get_index,CONVCSC_get_next_dc,FALSE);
  JSON_list_int(temp,output,"vsc_converters",bus,ConvVSC,BUSDC_get_vsc_conv,CONVVSC_get_index,CONVVSC_get_next_dc,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

int BUSDC_get_num_csc_convs(BusDC* bus, BOOL only_in_service) {
  ConvCSC* conv;
  int n = 0;
  if (bus) {
    for (conv = bus->csc_conv; conv != NULL; conv = CONVCSC_get_next_dc(conv)) {
      if (CONVCSC_is_in_service(conv) || !only_in_service)
        n++;
    }
    return n;
  }
  else
    return 0;
}

int BUSDC_get_num_vsc_convs(BusDC* bus, BOOL only_in_service) {
  ConvVSC* conv;
  int n = 0;
  if (bus) {
    for (conv = bus->vsc_conv; conv != NULL; conv = CONVVSC_get_next_dc(conv)) {
      if (CONVVSC_is_in_service(conv) || !only_in_service)
        n++;
    }
    return n;
  }
  else
    return 0;
}

int BUSDC_get_num_periods(BusDC* bus) {
  if (bus)
    return bus->num_periods;
  else
    return 0;
}

int BUSDC_get_num_vars(void* vbus, unsigned char var, int t_start, int t_end) {

  // Local vars
  BusDC* bus = (BusDC*)vbus;
  int num_vars = 0;
  int dt;

  // Checks
  if (!bus)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > bus->num_periods-1)
    t_end = bus->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & BUSDC_VAR_V) && (bus->vars & BUSDC_VAR_V))
    num_vars += dt;
  return num_vars;
}

int BUSDC_get_number(BusDC* bus) {
  if (bus)
    return bus->number;
  else
    return 0;
}

char* BUSDC_get_name(BusDC* bus) {
  if (bus)
    return bus->name;
  else
    return NULL;
}

char BUSDC_get_obj_type(void* bus) {
  if (bus)
    return OBJ_BUSDC;
  else
    return OBJ_UNKNOWN;
}

REAL BUSDC_get_P_mis(BusDC* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->P_mis[t];
  else
    return 0;
}

REAL BUSDC_get_v(BusDC* bus, int t) {
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  else
    return bus->v[t];
}

REAL BUSDC_get_v_base(BusDC* bus) {
  if (bus)
    return bus->v_base;
  else
    return 0;
}

Vec* BUSDC_get_var_indices(void* vbus, unsigned char var, int t_start, int t_end) {

  // Local vars
  BusDC* bus = (BusDC*)vbus;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!bus)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > bus->num_periods-1)
    t_end = bus->num_periods-1;

  // Indices
  indices = VEC_new(BUSDC_get_num_vars(vbus,var,t_start,t_end));
  if ((var & BUSDC_VAR_V) && (bus->vars & BUSDC_VAR_V)) { // voltage
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,bus->index_v[t]);
      offset++;
    }
  }
  return indices;
}

char* BUSDC_get_var_info_string(BusDC* bus, int index) {

  // Local variables
  char* info;

  //Check
  if (!bus)
    return NULL;

  // Voltage
  if ((bus->vars & BUSDC_VAR_V) &&
      index >= bus->index_v[0] &&
      index <= bus->index_v[bus->num_periods-1]) {
    info = (char*)malloc(BUSDC_BUFFER_SIZE*sizeof(char));
    snprintf(info,BUSDC_BUFFER_SIZE*sizeof(char),
             "dc_bus:%d:voltage:%d",bus->index,index-bus->index_v[0]);
    return info;
  }

  // Return
  return NULL;
}

void BUSDC_get_var_values(BusDC* bus, Vec* values, int code) {

  // Local vars
  int t;
  
  // No bus
  if (!bus)
    return;

  for (t = 0; t < bus->num_periods; t++) {
    
    // Voltage
    if (bus->vars & BUSDC_VAR_V) {
      switch (code) {

      case UPPER_LIMITS:
        VEC_set(values,bus->index_v[t],BUSDC_INF_V);
        break;

      case LOWER_LIMITS:
        VEC_set(values,bus->index_v[t],-BUSDC_INF_V);
        break;

      default:
        VEC_set(values,bus->index_v[t],bus->v[t]);
      }
    }
  }
}

int BUSDC_get_di_index(BusDC* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->di_index[t];
  else
    return -1;
}

BOOL BUSDC_has_flags(void* vbus, char flag_type, unsigned char mask) {
  BusDC* bus = (BusDC*)vbus;
  if (bus) {
    if (flag_type == FLAG_VARS)
      return (bus->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (bus->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (bus->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (bus->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL BUSDC_has_properties(void* vbus, char prop) {
  BusDC* bus = (BusDC*)vbus;
  if (!bus)
    return FALSE;
  return TRUE;
}

BusDC* BUSDC_hash_number_add(BusDC* bus_hash,BusDC* bus) {
  HASH_ADD(hh_number,bus_hash,number,sizeof(int),bus);
  return bus_hash;
}

void BUSDC_hash_number_del(BusDC* bus_hash) {
  while (bus_hash != NULL)
    HASH_DELETE(hh_number,bus_hash,bus_hash);
}

BusDC* BUSDC_hash_number_find(BusDC* bus_hash,int number) {
  BusDC* bus;
  HASH_FIND(hh_number,bus_hash,&number,sizeof(int),bus);
  return bus;
}

int BUSDC_hash_number_len(BusDC* bus_hash) {
  return HASH_CNT(hh_number,bus_hash);
}

BusDC* BUSDC_hash_name_add(BusDC* bus_hash, BusDC* bus) {
  HASH_ADD(hh_name,bus_hash,name[0],strlen(bus->name),bus);
  return bus_hash;
}

void BUSDC_hash_name_del(BusDC* bus_hash) {
  while (bus_hash != NULL)
    HASH_DELETE(hh_name,bus_hash,bus_hash);
}

BusDC* BUSDC_hash_name_find(BusDC* bus_hash, char* name) {
  BusDC* bus;
  HASH_FIND(hh_name,bus_hash,name,(unsigned)strlen(name),bus);
  return bus;
}

int BUSDC_hash_name_len(BusDC* bus_hash) {
  return HASH_CNT(hh_name,bus_hash);
}

void BUSDC_init(BusDC* bus, int num_periods) {

  // Local vars
  int T;
  int t;
  int i;

  // No bus
  if (!bus)
    return;

  T = num_periods;
  bus->num_periods = num_periods;

  bus->number = 0;
  for (i = 0; i < BUSDC_BUFFER_SIZE; i++)
    bus->name[i] = 0;

  bus->v_base = 0.;

  bus->index = -1;

  bus->in_service = TRUE;
  bus->fixed = 0x00;
  bus->bounded = 0x00;
  bus->sparse = 0x00;
  bus->vars = 0x00;

  bus->csc_conv = NULL;
  bus->vsc_conv = NULL;
  bus->branch_k = NULL;
  bus->branch_m = NULL;

  bus->net = NULL;

  ARRAY_zalloc(bus->v,REAL,T);

  ARRAY_zalloc(bus->index_v,int,T);

  ARRAY_zalloc(bus->P_mis,REAL,T);
  ARRAY_zalloc(bus->di_index,int,T);
  
  for (t = 0; t < bus->num_periods; t++)
    bus->v[t] = 1.;
}

BOOL BUSDC_is_in_service(void* bus) {
  if (bus)
    return ((BusDC*)bus)->in_service;
  else
    return FALSE;
}

BOOL BUSDC_is_equal(BusDC* bus, BusDC* other) {
  return bus == other;
}

BusDC* BUSDC_new(int num_periods) {
  if (num_periods > 0) {
    BusDC* bus = (BusDC*)malloc(sizeof(BusDC));
    BUSDC_init(bus,num_periods);
    return bus;
  }
  else
    return NULL;
}

void BUSDC_propagate_data_in_time(BusDC* bus, int start, int end) {
  int t;
  if (bus) {
    if (start < 0)
      start = 0;
    if (end > bus->num_periods)
      end = bus->num_periods;
    for (t = start+1; t < end; t++) {
      bus->v[t] = bus->v[start];
    }
  }
}

int BUSDC_set_flags(void* vbus, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  BusDC* bus = (BusDC*)vbus;
  int t;

  // Check bus
  if (!bus)
    return index;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(bus->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(bus->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(bus->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(bus->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & BUSDC_VAR_V) && (mask & BUSDC_VAR_V)) { // voltage
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < bus->num_periods; t++)
        bus->index_v[t] = index+t;
    }
    (*flags_ptr) |= BUSDC_VAR_V;
    index += bus->num_periods;
  }
  return index;
}

void BUSDC_set_index(BusDC* bus, int index) {
  if (bus)
    bus->index = index;
}

void BUSDC_set_number(BusDC* bus, int number) {
  if (bus)
    bus->number = number;
}

void BUSDC_set_network(BusDC* bus, void* net) {
  if (bus)
    bus->net = (Net*)net;
}

void BUSDC_set_in_service(BusDC* bus, BOOL in_service) {
  if (bus) {
    if (bus->in_service != in_service)
      NET_inc_state_tag(bus->net);
    bus->in_service = in_service;
  }
}

void BUSDC_set_name(BusDC* bus, char* name) {
  if (bus)
    strncpy(bus->name,name,(size_t)(BUSDC_BUFFER_SIZE-1));
}

void BUSDC_set_v(BusDC* bus, REAL v, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->v[t] = v;
}

void BUSDC_set_v_base(BusDC* bus, REAL v_base) {
  if (bus)
    bus->v_base = v_base;
}

void BUSDC_set_var_values(BusDC* bus, Vec* values) {

  // Local vars
  int t;

  // No bus
  if (!bus)
    return;

  // Time loop
  for (t = 0; t < bus->num_periods; t++) {

    // Voltage
    if (bus->vars & BUSDC_VAR_V)      // voltage (p.u.)
      bus->v[t] = VEC_get(values,bus->index_v[t]);
  }
}

void BUSDC_set_di_index(BusDC* bus, int idx, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->di_index[t] = idx;
}
