/** @file branch_dc.c
 *  @brief This file defines the BranchDC data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/branch_dc.h>
#include <pfnet/bus_dc.h>
#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

// BranchDC
struct BranchDC {

  // Properties
  char name[BRANCHDC_BUFFER_SIZE]; /**< @brief Branch name */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Buses
  BusDC* bus_k;      /**< @brief Bus connected to the "k" side */
  BusDC* bus_m;      /**< @brief Bus connected to the "m" side */
  
  // Resistance
  REAL r;            /**< @brief Series resistance (p.u.) */

  // Flags
  BOOL pre_cont_status; /**< @brief Flag for indicating whether the branch was in service before applying the contingency */
  BOOL in_service;   /**< @brief Flag for indicating whether the branch is in service */
  char vars;         /**< @brief Flags for indicating which quantities should be treated as variables */
  char fixed;        /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;      /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;       /**< @brief Flags for indicating which control adjustments should be sparse */

  // Indices
  int index;         /**< @brief Branch index */

  // Network
  Net* net; /**< @brief Network. */

  // List
  BranchDC* next_k;     /**< @brief List of branches connected to a bus on the "k" side */
  BranchDC* next_m;     /**< @brief List of branches connected to a bus in the "m" side */
};

void BRANCHDC_array_del(BranchDC* br_array, int size) {
  int i;
  BranchDC* br;
  if (br_array) {
    for (i = 0; i < size; i++) {
      br = &(br_array[i]);
      BRANCHDC_set_bus_k(br,NULL);
      BRANCHDC_set_bus_m(br,NULL);
    }
    free(br_array);
  }
}

void* BRANCHDC_array_get(void* br_array, int index) {
  if (br_array)
    return (void*)&(((BranchDC*)br_array)[index]);
  else
    return NULL;
}

BranchDC* BRANCHDC_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    BranchDC* br_array = (BranchDC*)malloc(sizeof(BranchDC)*size);
    for (i = 0; i < size; i++) {
      BRANCHDC_init(&(br_array[i]),num_periods);
      BRANCHDC_set_index(&(br_array[i]),i);
      snprintf(br_array[i].name,(size_t)(BRANCHDC_BUFFER_SIZE-1),"%d",i);
    }
    return br_array;
  }
  else
    return NULL;
}

void BRANCHDC_clear_flags(BranchDC* br, char flag_mask) {
  if (br) {
    if (flag_mask & FLAG_VARS)
      br->vars = 0x00;
    if (flag_mask & FLAG_BOUNDED)
      br->bounded = 0x00;
    if (flag_mask & FLAG_FIXED)
      br->fixed = 0x00;
    if (flag_mask & FLAG_SPARSE)
      br->sparse = 0x00;
  }
}

void BRANCHDC_clear_sensitivities(BranchDC* br) {
  // nothing
}

void BRANCHDC_copy_from_dc_branch(BranchDC* br, BranchDC* other) {

  // Check
  if (!br || !other)
    return;

  // Properties
  strcpy(br->name,other->name);

  // Time
  // skip time

  // Buses
  // skip buses

  // Resistance
  br->r = other->r;

  // Flags
  br->pre_cont_status = other->pre_cont_status;
  br->in_service = other->in_service;
  br->fixed = other->fixed;
  br->bounded = other->bounded;
  br->sparse = other->sparse;
  br->vars = other->vars;
  
  // Indices
  // skip index
  
  // List
  // skip next
}


BusDC* BRANCHDC_get_bus_k(BranchDC* br) {
  if (br)
    return br->bus_k;
  else
    return NULL;
}

BusDC* BRANCHDC_get_bus_m(BranchDC* br) {
  if (br)
    return br->bus_m;
  else
    return NULL;
}

char BRANCHDC_get_flags_vars(BranchDC* br) {
  if (br)
    return br->vars;
  else
    return 0;
}

char BRANCHDC_get_flags_fixed(BranchDC* br) {
  if (br)
    return br->fixed;
  else
    return 0;
}

char BRANCHDC_get_flags_bounded(BranchDC* br) {
  if (br)
    return br->bounded;
  else
    return 0;
}

char BRANCHDC_get_flags_sparse(BranchDC* br) {
  if (br)
    return br->sparse;
  else
    return 0;
}

REAL BRANCHDC_get_i_km(BranchDC* br, Vec* var_values, int t) {

  // Local variables
  BusDC* bus_k;
  BusDC* bus_m;
  REAL vk;
  REAL vm;

  if (br) {
    bus_k = BRANCHDC_get_bus_k(br);
    bus_m = BRANCHDC_get_bus_m(br);
    if (BUSDC_has_flags(bus_k,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vk = VEC_get(var_values,BUSDC_get_index_v(bus_k,t));
    else
      vk = BUSDC_get_v(bus_k,t);
    if (BUSDC_has_flags(bus_m,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vm = VEC_get(var_values,BUSDC_get_index_v(bus_m,t));
    else
      vm = BUSDC_get_v(bus_m,t);
    return (vk-vm)/br->r;
  }
  else
    return 0;  
}

REAL BRANCHDC_get_i_mk(BranchDC* br, Vec* var_values, int t) {
  return -BRANCHDC_get_i_km(br,var_values,t);
}

REAL BRANCHDC_get_P_km(BranchDC* br, Vec* var_values, int t) {

  // Local variables
  BusDC* bus_k;
  BusDC* bus_m;
  REAL vk;
  REAL vm;

  if (br) {
    bus_k = BRANCHDC_get_bus_k(br);
    bus_m = BRANCHDC_get_bus_m(br);
    if (BUSDC_has_flags(bus_k,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vk = VEC_get(var_values,BUSDC_get_index_v(bus_k,t));
    else
      vk = BUSDC_get_v(bus_k,t);
    if (BUSDC_has_flags(bus_m,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vm = VEC_get(var_values,BUSDC_get_index_v(bus_m,t));
    else
      vm = BUSDC_get_v(bus_m,t);
    return vk*(vk-vm)/br->r;
  }
  else
    return 0; 
}

REAL BRANCHDC_get_P_mk(BranchDC* br, Vec* var_values, int t) {

  // Local variables
  BusDC* bus_k;
  BusDC* bus_m;
  REAL vk;
  REAL vm;

  if (br) {
    bus_k = BRANCHDC_get_bus_k(br);
    bus_m = BRANCHDC_get_bus_m(br);
    if (BUSDC_has_flags(bus_k,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vk = VEC_get(var_values,BUSDC_get_index_v(bus_k,t));
    else
      vk = BUSDC_get_v(bus_k,t);
    if (BUSDC_has_flags(bus_m,FLAG_VARS,BUSDC_VAR_V) && var_values)
      vm = VEC_get(var_values,BUSDC_get_index_v(bus_m,t));
    else
      vm = BUSDC_get_v(bus_m,t);
    return vm*(vm-vk)/br->r;
  }
  else
    return 0;
}

int BRANCHDC_get_index(BranchDC* br) {
  if (br)
    return br->index;
  else
    return -1;
}

char* BRANCHDC_get_json_string(BranchDC* branch, char* output) {

  // Local variables
  char temp[BRANCHDC_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No branch
  if (!branch)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*BRANCHDC_BUFFER_SIZE*BRANCHDC_NUM_JSON_FIELDS*branch->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",branch->index,FALSE);
  JSON_int(temp,output,"num_periods",branch->num_periods,FALSE);
  JSON_str(temp,output,"name",branch->name,FALSE);
  JSON_bool(temp,output,"in_service",branch->in_service,FALSE);
  JSON_obj(temp,output,"bus_k",branch->bus_k,BUSDC_get_index,FALSE);
  JSON_obj(temp,output,"bus_m",branch->bus_m,BUSDC_get_index,FALSE);
  JSON_float(temp,output,"r",branch->r,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

char* BRANCHDC_get_name(BranchDC* br) {
  if (br)
    return br->name;
  else
    return NULL;
}

BranchDC* BRANCHDC_get_next_k(BranchDC* br) {
  if (br)
    return br->next_k;
  else
    return NULL;
}

BranchDC* BRANCHDC_get_next_m(BranchDC* br) {
  if (br)
    return br->next_m;
  else
    return NULL;
}

int BRANCHDC_get_num_periods(BranchDC* br) {
  if (br)
    return br->num_periods;
  else
    return 0;
}

int BRANCHDC_get_num_vars(void* vbr, unsigned char var, int t_start, int t_end) {
  // nothing
  return 0;
}

char BRANCHDC_get_obj_type(void* br) {
  if (br)
    return OBJ_BRANCHDC;
  else
    return OBJ_UNKNOWN;
}

REAL BRANCHDC_get_r(BranchDC* br) {
  if (br)
    return br->r;
  else
    return 0;
}

Vec* BRANCHDC_get_var_indices(void* vbr, unsigned char var, int t_start, int t_end) {
  // nothing
  return NULL;
}

char* BRANCHDC_get_var_info_string(BranchDC* br, int index) {
  // nothing
  return NULL;
}

void BRANCHDC_get_var_values(BranchDC* br, Vec* values, int code) {
  // nothing
  return;
}

BOOL BRANCHDC_has_flags(void* vbr, char flag_type, unsigned char mask) {
  BranchDC* br = (BranchDC*)vbr;
  if (br) {
    if (flag_type == FLAG_VARS)
      return (br->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (br->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (br->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (br->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL BRANCHDC_has_properties(void* vbr, char prop) {
  BranchDC* br = (BranchDC*)vbr;
  if (!br)
    return FALSE;
  return TRUE;
}

void BRANCHDC_init(BranchDC* br, int num_periods) {

  // No branch
  if (!br)
    return;

  br->num_periods = num_periods;

  ARRAY_clear(br->name,char,BRANCHDC_BUFFER_SIZE);

  br->bus_k = NULL;
  br->bus_m = NULL;

  br->r = 0;

  br->pre_cont_status = PRE_CONT_UNSET;
  br->in_service = TRUE;
  br->vars = 0x00;
  br->fixed = 0x00;
  br->bounded = 0x00;
  br->sparse = 0x00;

  br->index = -1;

  br->net = NULL;

  br->next_k = NULL;
  br->next_m = NULL;
};

BOOL BRANCHDC_is_in_service(void* br) {
  if (br)
    return (((BranchDC*)br)->in_service &&
            BUSDC_is_in_service(((BranchDC*)br)->bus_k) &&
            BUSDC_is_in_service(((BranchDC*)br)->bus_m));
  else
    return FALSE;
}

BOOL BRANCHDC_is_equal(BranchDC* br, BranchDC* other) {
  return br == other;
}

BranchDC* BRANCHDC_list_k_add(BranchDC* k_br_list, BranchDC* br) {
  LIST_add(BranchDC,k_br_list,br,next_k);
  return k_br_list;
}

BranchDC* BRANCHDC_list_k_del(BranchDC* k_br_list, BranchDC* br) {
  LIST_del(BranchDC,k_br_list,br,next_k);
  return k_br_list;
}

int BRANCHDC_list_k_len(BranchDC* k_br_list) {
  int len;
  LIST_len(BranchDC,k_br_list,next_k,len);
  return len;
}

BranchDC* BRANCHDC_list_m_add(BranchDC* m_br_list, BranchDC* br) {
  LIST_add(BranchDC,m_br_list,br,next_m);
  return m_br_list;
}

BranchDC* BRANCHDC_list_m_del(BranchDC* m_br_list, BranchDC* br) {
  LIST_del(BranchDC,m_br_list,br,next_m);
  return m_br_list;
}

int BRANCHDC_list_m_len(BranchDC* m_br_list) {
  int len;
  LIST_len(BranchDC,m_br_list,next_m,len);
  return len;
}

BranchDC* BRANCHDC_new(int num_periods) {
  if (num_periods > 0) {
    BranchDC* branch = (BranchDC*)malloc(sizeof(BranchDC));
    BRANCHDC_init(branch,num_periods);
    return branch;
  }
  else
    return NULL;
}

void BRANCHDC_propagate_data_in_time(BranchDC* br, int start, int end) {
  // nothing
}

void BRANCHDC_set_in_service(BranchDC* br, BOOL in_service) {
  if (br && BUSDC_is_in_service(br->bus_k) && BUSDC_is_in_service(br->bus_m)) {
    if (br->in_service != in_service)
      NET_inc_state_tag(br->net);
    br->in_service = in_service;
  }
}

void BRANCHDC_set_network(BranchDC* br, void* net) {
  if (br)
    br->net = (Net*)net;
}

void BRANCHDC_set_bus_k(BranchDC* br, BusDC* bus_k) {
  BusDC* old_bus_k;
  if (br) {
    old_bus_k = br->bus_k;
    br->bus_k = NULL;
    BUSDC_del_branch_k(old_bus_k,br);
    br->bus_k = bus_k;
    BUSDC_add_branch_k(br->bus_k,br);
  }
}

void BRANCHDC_set_bus_m(BranchDC* br, BusDC* bus_m) {
  BusDC* old_bus_m;
  if (br) {
    old_bus_m = br->bus_m;
    br->bus_m = NULL;
    BUSDC_del_branch_m(old_bus_m,br);
    br->bus_m = bus_m;
    BUSDC_add_branch_m(br->bus_m,br);
  }
}

int BRANCHDC_set_flags(void* vbr, char flag_type, unsigned char mask, int index) {
  // nothing
  return index;
}

void BRANCHDC_set_index(BranchDC* br, int index) {
  if (br)
    br->index = index;
}

void BRANCHDC_set_name(BranchDC* br, char* name) {
  if (br)
    strncpy(br->name,name,(size_t)(BRANCHDC_BUFFER_SIZE-1));
}

void BRANCHDC_set_r(BranchDC* br, REAL r) {
  if (br)
    br->r = r;
}

void BRANCHDC_set_var_values(BranchDC* br, Vec* values) {
  // nothing
}
