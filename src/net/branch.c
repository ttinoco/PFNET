/** @file branch.c
 *  @brief This file defines the Branch data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/branch.h>
#include <pfnet/bus.h>

// Branch
struct Branch {

  // Properties
  char type;         /**< @brief %Branch type */
  char obj_type;     /**< @brief %Object type */

  // Buses
  Bus* bus_from;     /**< @brief Bus connected to the "from" side */
  Bus* bus_to;       /**< @brief Bus connected to the "to" side */
  Bus* reg_bus;      /**< @brief Bus regulated by this transformer */

  // Conductance
  REAL g;            /**< @brief Series conductance (p.u.) */
  REAL g_from;       /**< @brief %Shunt conductance on "from" side (p.u.) */
  REAL g_to;         /**< @brief %Shunt conductance on "to" side (p.u.) */

  // Susceptance
  REAL b;            /**< @brief Series susceptance (p.u.) */
  REAL b_from;       /**< @brief %Shunt susceptance on "from" side (p.u.) */
  REAL b_to;         /**< @brief %Shunt shunt susceptance "to" side (p.u.) */
  
  // Tap ratio
  REAL ratio;        /**< @brief Transformer taps ratio (p.u.) */
  REAL ratio_max;    /**< @brief Maximum transformer taps ratio (p.u.) */ 
  REAL ratio_min;    /**< @brief Minimum transformer taps ratio (p.u.) */
  char num_ratios;   /**< @brief Number of tap positions */

  // Phase shift
  REAL phase;        /**< @brief Transformer phase shift (radians) */
  REAL phase_max;    /**< @brief Maximum transformer phase shift (radians) */
  REAL phase_min;    /**< @brief Minimum transformer phase shift (radians) */
  
  // Flow bounds
  REAL P_max;        /**< @brief Maximum active power flow (p.u.) */
  REAL P_min;        /**< @brief Minimum active power flow (p.u.) */
  REAL Q_max;        /**< @brief Maximum reactive power flow (p.u.) */
  REAL Q_min;        /**< @brief Minimum reactive power flow (p.u.) */

  // Power ratings
  REAL ratingA;      /**< @brief Power rating A (p.u. system base MVA) */
  REAL ratingB;      /**< @brief Power rating B (p.u. system base MVA) */
  REAL ratingC;      /**< @brief Power rating C (p.u. system base MVA) */  
   
  // Flags
  BOOL pos_ratio_v_sens; /**< @brief Flag for positive ratio-voltage sensitivity */
  char vars;             /**< @brief Flags for indicating which quantities should be treated as variables */
  char fixed;            /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;          /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;           /**< @brief Flags for indicating which control adjustments should be sparse */

  // Indices
  int index;         /**< @brief Branch index */
  int index_ratio;   /**< @brief Taps ratio index */
  int index_ratio_y; /**< @brief Taps ratio positive deviation index */
  int index_ratio_z; /**< @brief Taps ratio negative deviation index */
  int index_phase;   /**< @brief Phase shift index */
  int index_P;       /**< @brief Branch active power flow index */
  int index_Q;       /**< @brief Branch reactive power flow index */

  // List
  Branch* reg_next;  /**< @brief List of branches regulating a bus voltage magnitude */
  Branch* from_next; /**< @brief List of branches connected to a bus on the "from" side */
  Branch* to_next;   /**< @brief List of branches connected to a bus in the "to" side */
};

void* BRANCH_array_get(void* branch, int index) {
  if (branch)
    return (void*)&(((Branch*)branch)[index]);
  else
    return NULL;
}

Branch* BRANCH_array_new(int num) {
  int i;
  Branch* branch = (Branch*)malloc(sizeof(Branch)*num);
  for (i = 0; i < num; i++) {
    BRANCH_init(&(branch[i]));
    BRANCH_set_index(&(branch[i]),i);
  }
  return branch;
}

void BRANCH_array_show(Branch* branch, int num) {
  int i;
  for (i = 0; i < num; i++)
    BRANCH_show(&(branch[i]));
}

void BRANCH_clear_flags(Branch* br, char flag_type) {
  if (br) {
    if (flag_type == FLAG_VARS)
      br->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      br->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      br->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      br->sparse = 0x00;
  }
}

char BRANCH_get_obj_type(void* b) {
  if (b)
    return ((Branch*)b)->obj_type;
  else
    return OBJ_UNKNOWN;
}

int BRANCH_get_index(Branch* b) {
  if (b)
    return b->index;
  else
    return 0;
}

int BRANCH_get_index_ratio(Branch* b) {
  if (b)
    return b->index_ratio;
  else
    return 0;
}

int BRANCH_get_index_ratio_y(Branch* b) {
  if (b)
    return b->index_ratio_y;
  else
    return 0;
}

int BRANCH_get_index_ratio_z(Branch* b) {
  if (b)
    return b->index_ratio_z;
  else
    return 0;
}

int BRANCH_get_index_phase(Branch* b) {
  if (b)
    return b->index_phase;
  else
    return 0;
}

REAL BRANCH_get_ratio(Branch* b) {
  if (b)
    return b->ratio;
  else
    return 0;
}

REAL BRANCH_get_ratio_max(Branch* b) {
  if (b)
    return b->ratio_max;
  else
    return 0;
}

REAL BRANCH_get_ratio_min(Branch* b) {
  if (b)
    return b->ratio_min;
  else
    return 0;
}

REAL BRANCH_get_b(Branch* b) {
  if (b)
    return b->b;
  else
    return 0;
}

REAL BRANCH_get_b_from(Branch* b) {
  if (b)
    return b->b_from;
  else
    return 0;
}

REAL BRANCH_get_b_to(Branch* b) {
  if (b)
    return b->b_to;
  else
    return 0;
}

REAL BRANCH_get_g(Branch* b) {
  if (b)
    return b->g;
  else
    return 0;
}

REAL BRANCH_get_g_from(Branch* b) {
  if (b)
    return b->g_from;
  else
    return 0;
}

REAL BRANCH_get_g_to(Branch* b) {
  if (b)
    return b->g_to;
  else
    return 0;
}

Bus* BRANCH_get_bus_from(Branch* br) {
  if (br)
    return br->bus_from;
  else
    return NULL;
}

Bus* BRANCH_get_bus_to(Branch* br) {
  if (br)
    return br->bus_to;
  else
    return NULL;
}

Bus* BRANCH_get_reg_bus(Branch* br) {
  if (br)
    return br->reg_bus;
  else
    return NULL;
}

Branch* BRANCH_get_reg_next(Branch* br) {
  if (br)
    return br->reg_next;
  else
    return NULL;
}

Branch* BRANCH_get_from_next(Branch* br) {
  if (br)
    return br->from_next;
  else
    return NULL;
}

Branch* BRANCH_get_to_next(Branch* br) {
  if (br)
    return br->to_next;
  else
    return NULL;
}

REAL BRANCH_get_phase(Branch* b) {
  if (b)
    return b->phase;
  else
    return 0;
}

REAL BRANCH_get_phase_max(Branch* b) {
  if (b)
    return b->phase_max;
  else
    return 0;
}

REAL BRANCH_get_phase_min(Branch* b) {
  if (b)
    return b->phase_min;
  else
    return 0;
}

REAL BRANCH_get_ratingA(Branch* b) {
  if (b)
    return b->ratingA;
  else
    return 0;
}

REAL BRANCH_get_ratingB(Branch* b) {
  if (b)
    return b->ratingB;
  else
    return 0;
}

REAL BRANCH_get_ratingC(Branch* b) {
  if (b)
    return b->ratingC;
  else
    return 0;
}

void BRANCH_get_var_values(Branch* br, Vec* values, int code) {

  // No branch
  if (!br)
    return;

  if (br->vars & BRANCH_VAR_RATIO) { // taps ratio
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,br->index_ratio,br->ratio_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,br->index_ratio,br->ratio_min);
      break;
    default:
      VEC_set(values,br->index_ratio,br->ratio);
    }
  }
  if (br->vars & BRANCH_VAR_PHASE) { // phase shift
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,br->index_phase,br->phase_max);
      break;
    case LOWER_LIMITS:
      VEC_set(values,br->index_phase,br->phase_min);
      break;
    default:
      VEC_set(values,br->index_phase,br->phase);
    }
  }
  if (br->vars & BRANCH_VAR_RATIO_DEV) { // tap ratio deviations
    switch(code) {
    case UPPER_LIMITS:
      VEC_set(values,br->index_ratio_y,INF);
      VEC_set(values,br->index_ratio_z,INF);
      break;
    case LOWER_LIMITS:
      VEC_set(values,br->index_ratio_y,-INF);
      VEC_set(values,br->index_ratio_z,-INF);
      break;
    default:
      VEC_set(values,br->index_ratio_y,0.);
      VEC_set(values,br->index_ratio_z,0.);
    }
  }   
}

int BRANCH_get_var_index(void* vbr, char var) {
  Branch* br = (Branch*)vbr;
  if (!br)
    return 0;
  if (var == BRANCH_VAR_RATIO)
    return br->index_ratio;
  if (var == BRANCH_VAR_PHASE)
    return br->index_phase;
  return 0;
}

BOOL BRANCH_has_pos_ratio_v_sens(Branch* branch) {
  if (branch)
    return branch->pos_ratio_v_sens;
  else
    return FALSE;
}

BOOL BRANCH_has_flags(void* vbr, char flag_type, char mask) {
  Branch* br = (Branch*)vbr;
  if (br) {
    if (flag_type == FLAG_VARS)
      return (br->vars & mask);
    else if (flag_type == FLAG_BOUNDED)
      return (br->bounded & mask);
    else if (flag_type == FLAG_FIXED)
      return (br->fixed & mask);
    else if (flag_type == FLAG_SPARSE)
      return (br->sparse & mask);
    return FALSE;
  }
  else
    return FALSE;
}

BOOL BRANCH_has_properties(void* vbr, char prop) {
  Branch* br = (Branch*)vbr;
  if (!br)
    return FALSE;
  if ((prop & BRANCH_PROP_TAP_CHANGER) && !BRANCH_is_tap_changer(br))
    return FALSE;
  if ((prop & BRANCH_PROP_TAP_CHANGER_V) && !BRANCH_is_tap_changer_v(br))
    return FALSE;
  if ((prop & BRANCH_PROP_TAP_CHANGER_Q) && !BRANCH_is_tap_changer_Q(br))
    return FALSE;
  if ((prop & BRANCH_PROP_PHASE_SHIFTER) && !BRANCH_is_phase_shifter(br))
    return FALSE;
  return TRUE;
}

void BRANCH_init(Branch* br) {  

  br->type = BRANCH_TYPE_LINE;
  br->obj_type = OBJ_BRANCH;

  br->bus_from = NULL;
  br->bus_to = NULL;
  br->reg_bus = NULL;

  br->g = 0;
  br->g_from = 0;
  br->g_to = 0;
  br->b = 0;
  br->b_from = 0;
  br->b_to = 0;

  br->ratio = 1;
  br->ratio_max = 1;
  br->ratio_min = 1;
  br->num_ratios = 1;

  br->phase = 0;
  br->phase_max = 0;
  br->phase_min = 0;

  br->P_max = 0;
  br->P_min = 0;
  br->Q_max = 0;
  br->Q_min = 0;

  br->ratingA = 0;
  br->ratingB = 0;
  br->ratingC = 0;

  br->pos_ratio_v_sens = TRUE;
  br->vars = 0x00;
  br->fixed = 0x00;
  br->bounded = 0x00;
  br->sparse = 0x00;

  br->index = 0;
  br->index_ratio = 0;
  br->index_ratio_y = 0;
  br->index_ratio_z = 0;
  br->index_phase = 0;
  br->index_P = 0;
  br->index_Q = 0;

  br->reg_next = NULL;
  br->from_next = NULL;
  br->to_next = NULL;
};

BOOL BRANCH_is_fixed_tran(Branch* branch) {
  return (branch->type == BRANCH_TYPE_TRAN_FIXED);
}

BOOL BRANCH_is_line(Branch* branch) {
  return (branch->type == BRANCH_TYPE_LINE);
}

BOOL BRANCH_is_phase_shifter(Branch* branch) {
  return (branch->type == BRANCH_TYPE_TRAN_PHASE);
}

BOOL BRANCH_is_tap_changer(Branch* br) {
  return (BRANCH_is_tap_changer_v(br) || BRANCH_is_tap_changer_Q(br));
}

BOOL BRANCH_is_tap_changer_v(Branch* br) {
  if (br)
    return (br->type == BRANCH_TYPE_TRAN_TAP_V);
  else
    return FALSE;
}

BOOL BRANCH_is_tap_changer_Q(Branch* br) {
  if (br)
    return (br->type == BRANCH_TYPE_TRAN_TAP_Q);
  else
    return FALSE;
}

Branch* BRANCH_list_reg_add(Branch* reg_br_list, Branch* reg_br) {
  LIST_add(reg_br_list,reg_br,reg_next);
  return reg_br_list;
}

int BRANCH_list_reg_len(Branch* reg_br_list) {
  int len;
  LIST_len(Branch,reg_br_list,reg_next,len);
  return len;
}

Branch* BRANCH_list_from_add(Branch* from_br_list, Branch* br) {
  LIST_add(from_br_list,br,from_next);
  return from_br_list;
}

int BRANCH_list_from_len(Branch* from_br_list) {
  int len;
  LIST_len(Branch,from_br_list,from_next,len);
  return len;
}

Branch* BRANCH_list_to_add(Branch* to_br_list, Branch* br) {
  LIST_add(to_br_list,br,to_next);
  return to_br_list;
}

int BRANCH_list_to_len(Branch* to_br_list) {
  int len;
  LIST_len(Branch,to_br_list,to_next,len);
  return len;
}

Branch* BRANCH_new(void) {
  Branch* branch = (Branch*)malloc(sizeof(Branch));
  BRANCH_init(branch);
  return branch;
}

void BRANCH_set_index(Branch* br, int index) {
  if (br) 
    br->index = index;
}

void BRANCH_set_type(Branch* br, int type) {
  if (br)
    br->type = type;
}

void BRANCH_set_bus_from(Branch* branch, Bus* bus_from) {
  branch->bus_from = (Bus*)bus_from;
}

void BRANCH_set_bus_to(Branch* branch, Bus* bus_to) {
  branch->bus_to = (Bus*)bus_to;
}

void BRANCH_set_reg_bus(Branch* branch, Bus* reg_bus) {
  branch->reg_bus = (Bus*)reg_bus;
}

void BRANCH_set_g(Branch* branch, REAL g) {
  branch->g = g;
}

void BRANCH_set_g_from(Branch* branch, REAL g_from) {
  branch->g_from = g_from;
}

void BRANCH_set_g_to(Branch* branch, REAL g_to) {
  branch->g_to = g_to;
}

void BRANCH_set_b(Branch* branch, REAL b) {
  branch->b = b;
}

void BRANCH_set_b_from(Branch* branch, REAL b_from) {
  branch->b_from = b_from;
}

void BRANCH_set_b_to(Branch* branch, REAL b_to) {
  branch->b_to = b_to;
}

void BRANCH_set_ratio(Branch* branch, REAL ratio) {
  branch->ratio = ratio;
}

void BRANCH_set_ratio_max(Branch* branch, REAL ratio) {
  branch->ratio_max = ratio;
}

void BRANCH_set_ratio_min(Branch* branch, REAL ratio) {
  branch->ratio_min = ratio;
}

void BRANCH_set_pos_ratio_v_sens(Branch* branch, BOOL flag) {
  if (branch)
    branch->pos_ratio_v_sens = flag;
}

void BRANCH_set_phase(Branch* br, REAL phase) {
  if (br)
    br->phase = phase;
}

void BRANCH_set_phase_max(Branch* br, REAL phase) {
  if (br)
    br->phase_max = phase;
}

void BRANCH_set_phase_min(Branch* br, REAL phase) {
  if (br)
    br->phase_min = phase;
}

void BRANCH_set_P_max(Branch* br, REAL P_max) {
  if (br)
    br->P_max = P_max;
}

void BRANCH_set_P_min(Branch* br, REAL P_min) {
  if (br)
    br->P_min = P_min;
}

void BRANCH_set_Q_max(Branch* br, REAL Q_max) {
  if (br)
    br->Q_max = Q_max;
}

void BRANCH_set_Q_min(Branch* br, REAL Q_min) {
  if (br)
    br->Q_min = Q_min;
}

void BRANCH_set_ratingA(Branch* br, REAL r) {
  if (br)
    br->ratingA = r;
}

void BRANCH_set_ratingB(Branch* br, REAL r) {
  if (br)
    br->ratingB = r;
}

void BRANCH_set_ratingC(Branch* br, REAL r) {
  if (br)
    br->ratingC = r;
}

void BRANCH_set_var_values(Branch* br, Vec* values) {

  // No branch
  if (!br)
    return;

  // Set variable values
  if (br->vars & BRANCH_VAR_RATIO)    // taps ratio
    br->ratio = VEC_get(values,br->index_ratio); 
  if (br->vars & BRANCH_VAR_PHASE)    // phase shift
    br->phase = VEC_get(values,br->index_phase);
}

int BRANCH_set_flags(void* vbr, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Branch* br = (Branch*)vbr;

  // Check branch
  if (!br)
    return index;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(br->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(br->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(br->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(br->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & BRANCH_VAR_RATIO) && (mask & BRANCH_VAR_RATIO)) { // taps ratio
    if (flag_type == FLAG_VARS)
      br->index_ratio = index;
    (*flags_ptr) |= BRANCH_VAR_RATIO;
    index++;
  }
  if (!((*flags_ptr) & BRANCH_VAR_PHASE) && (mask & BRANCH_VAR_PHASE)) { // phase shift
    if (flag_type == FLAG_VARS)
      br->index_phase = index;
    (*flags_ptr) |= BRANCH_VAR_PHASE;
    index++;
  }
  if (!((*flags_ptr) & BRANCH_VAR_RATIO_DEV) && (mask & BRANCH_VAR_RATIO_DEV)) { // taps ratio deviations
    if (flag_type == FLAG_VARS) {
      br->index_ratio_y = index;
      br->index_ratio_z = index+1;
    }
    (*flags_ptr) |= BRANCH_VAR_RATIO_DEV;
    index += 2;
  }
  return index;
}

void BRANCH_show(Branch* br) {
  printf("branch %d\t%d\t%d\n",
	 BUS_get_number(br->bus_from),
	 BUS_get_number(br->bus_to),
	 br->type);
}

