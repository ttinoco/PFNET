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
#include <pfnet/array.h>

// Branch
struct Branch {

  // Properties
  char type;         /**< @brief %Branch type */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

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
  REAL* ratio;       /**< @brief Transformer taps ratio (p.u.) */
  REAL ratio_max;    /**< @brief Maximum transformer taps ratio (p.u.) */ 
  REAL ratio_min;    /**< @brief Minimum transformer taps ratio (p.u.) */
  char num_ratios;   /**< @brief Number of tap positions */

  // Phase shift
  REAL* phase;       /**< @brief Transformer phase shift (radians) */
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
  BOOL outage;           /**< @brief Flag for indicating that branch in on outage */
  BOOL pos_ratio_v_sens; /**< @brief Flag for positive ratio-voltage sensitivity */
  char vars;             /**< @brief Flags for indicating which quantities should be treated as variables */
  char fixed;            /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;          /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;           /**< @brief Flags for indicating which control adjustments should be sparse */

  // Indices
  int index;          /**< @brief Branch index */
  int* index_ratio;   /**< @brief Taps ratio index */
  int* index_ratio_y; /**< @brief Taps ratio positive deviation index */
  int* index_ratio_z; /**< @brief Taps ratio negative deviation index */
  int* index_phase;   /**< @brief Phase shift index */

  // Sensitivities
  REAL* sens_P_u_bound;  /**< @brief Sensitivity of active power flow upper bound */
  REAL* sens_P_l_bound;  /**< @brief Sensitivity of active power flow lower bound */

  // List
  Branch* reg_next;  /**< @brief List of branches regulating a bus voltage magnitude */
  Branch* from_next; /**< @brief List of branches connected to a bus on the "from" side */
  Branch* to_next;   /**< @brief List of branches connected to a bus in the "to" side */
};

void* BRANCH_array_get(void* br_array, int index) {
  if (br_array)
    return (void*)&(((Branch*)br_array)[index]);
  else
    return NULL;
}

void BRANCH_array_del(Branch* br_array, int size) {
  int i;
  Branch* br;
  if (br_array) {
    for (i = 0; i < size; i++) {
      br = &(br_array[i]);
      free(br->ratio);
      free(br->phase);
      free(br->index_ratio);
      free(br->index_ratio_y);
      free(br->index_ratio_z);
      free(br->index_phase);
      free(br->sens_P_u_bound);
      free(br->sens_P_l_bound);
    }
    free(br_array);
  }
}

Branch* BRANCH_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    Branch* br_array = (Branch*)malloc(sizeof(Branch)*size);
    for (i = 0; i < size; i++) {
      BRANCH_init(&(br_array[i]),num_periods);
      BRANCH_set_index(&(br_array[i]),i);
    }
    return br_array;
  }
  else
    return NULL;
}

void BRANCH_array_show(Branch* br_array, int size, int t) {
  int i;
  if (br_array) {
    for (i = 0; i < size; i++)
      BRANCH_show(&(br_array[i]),t);
  }
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

void BRANCH_clear_sensitivities(Branch* br) {
  int t;
  if (br) {
    for (t = 0; t < br->num_periods; t++) {
      br->sens_P_u_bound[t] = 0;
      br->sens_P_l_bound[t] = 0;
    }
  }
}

int BRANCH_get_num_periods(Branch* br) {
  if (br)
    return br->num_periods;
  else
    return 0;
}

char BRANCH_get_type(Branch* br) {
  if (br)
    return br->type;
  else
    return BRANCH_TYPE_LINE;
}

char BRANCH_get_obj_type(void* br) {
  if (br)
    return OBJ_BRANCH;
  else
    return OBJ_UNKNOWN;
}

REAL BRANCH_get_sens_P_u_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_P_u_bound[t];
  else
    return 0;
}

REAL BRANCH_get_sens_P_l_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_P_l_bound[t];
  else
    return 0;
}

int BRANCH_get_index(Branch* br) {
  if (br)
    return br->index;
  else
    return 0;
}

int BRANCH_get_index_ratio(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_ratio[t];
  else
    return 0;
}

int BRANCH_get_index_ratio_y(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_ratio_y[t];
  else
    return 0;
}

int BRANCH_get_index_ratio_z(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_ratio_z[t];
  else
    return 0;
}

int BRANCH_get_index_phase(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_phase[t];
  else
    return 0;
}

REAL BRANCH_get_ratio(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->ratio[t];
  else
    return 0;
}

REAL BRANCH_get_ratio_max(Branch* br) {
  if (br)
    return br->ratio_max;
  else
    return 0;
}

REAL BRANCH_get_ratio_min(Branch* br) {
  if (br)
    return br->ratio_min;
  else
    return 0;
}

REAL BRANCH_get_b(Branch* br) {
  if (br)
    return br->b;
  else
    return 0;
}

REAL BRANCH_get_b_from(Branch* br) {
  if (br)
    return br->b_from;
  else
    return 0;
}

REAL BRANCH_get_b_to(Branch* br) {
  if (br)
    return br->b_to;
  else
    return 0;
}

REAL BRANCH_get_g(Branch* br) {
  if (br)
    return br->g;
  else
    return 0;
}

REAL BRANCH_get_g_from(Branch* br) {
  if (br)
    return br->g_from;
  else
    return 0;
}

REAL BRANCH_get_g_to(Branch* br) {
  if (br)
    return br->g_to;
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

REAL BRANCH_get_phase(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->phase[t];
  else
    return 0;
}

REAL BRANCH_get_phase_max(Branch* br) {
  if (br)
    return br->phase_max;
  else
    return 0;
}

REAL BRANCH_get_phase_min(Branch* br) {
  if (br)
    return br->phase_min;
  else
    return 0;
}

REAL BRANCH_get_ratingA(Branch* br) {
  if (br)
    return br->ratingA;
  else
    return 0;
}

REAL BRANCH_get_ratingB(Branch* br) {
  if (br)
    return br->ratingB;
  else
    return 0;
}

REAL BRANCH_get_ratingC(Branch* br) {
  if (br)
    return br->ratingC;
  else
    return 0;
}

REAL BRANCH_get_P_flow_DC(Branch* br, int t) {
  /* Active power flow (DC approx) from bus
     "from" to bus "to". */

  if (br && t >= 0 && t < br->num_periods) {
    return -(br->b)*(BUS_get_v_ang(br->bus_from,t)-
		     BUS_get_v_ang(br->bus_to,t)-
		     br->phase[t]);
  }
  else
    return 0;
}

void BRANCH_get_var_values(Branch* br, Vec* values, int code) {

  // Local vars
  int t;

  // No branch
  if (!br)
    return;

  // Time loop
  for (t = 0; t < br->num_periods; t++) {

    if (br->vars & BRANCH_VAR_RATIO) { // taps ratio
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,br->index_ratio[t],br->ratio_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,br->index_ratio[t],br->ratio_min);
	break;
      default:
	VEC_set(values,br->index_ratio[t],br->ratio[t]);
      }
    }
    if (br->vars & BRANCH_VAR_PHASE) { // phase shift
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,br->index_phase[t],br->phase_max);
	break;
      case LOWER_LIMITS:
	VEC_set(values,br->index_phase[t],br->phase_min);
	break;
      default:
	VEC_set(values,br->index_phase[t],br->phase[t]);
      }
    }
    if (br->vars & BRANCH_VAR_RATIO_DEV) { // tap ratio deviations
      switch(code) {
      case UPPER_LIMITS:
	VEC_set(values,br->index_ratio_y[t],BRANCH_INF_RATIO);
	VEC_set(values,br->index_ratio_z[t],BRANCH_INF_RATIO);
	break;
      case LOWER_LIMITS:
	VEC_set(values,br->index_ratio_y[t],0.);
	VEC_set(values,br->index_ratio_z[t],0.);
	break;
      default:
	VEC_set(values,br->index_ratio_y[t],0.);
	VEC_set(values,br->index_ratio_z[t],0.);
      }
    }
  }   
}

int BRANCH_get_num_vars(void* vbr, char var) {
  Branch* br = (Branch*)vbr;
  int num_vars = 0;
  if (!br)
    return 0;
  if ((var & BRANCH_VAR_RATIO) && (br->vars & BRANCH_VAR_RATIO))
    num_vars += br->num_periods;
  if ((var & BRANCH_VAR_PHASE) && (br->vars & BRANCH_VAR_PHASE))
    num_vars += br->num_periods;
  if ((var & BRANCH_VAR_RATIO_DEV) && (br->vars & BRANCH_VAR_RATIO_DEV))
    num_vars += 2*br->num_periods;
  return num_vars;
}

Vec* BRANCH_get_var_indices(void* vbr, char var) {
  Branch* br = (Branch*)vbr;
  Vec* indices;
  int num_vars;
  int offset = 0;
  int t;
  if (!br)
    return NULL;
  indices = VEC_new(BRANCH_get_num_vars(vbr,var));
  if ((var & BRANCH_VAR_RATIO) && (br->vars & BRANCH_VAR_RATIO)) {
    for (t = 0; t < br->num_periods; t++)
      VEC_set(indices,offset+t,br->index_ratio[t]);
    offset += br->num_periods;
  }
  if ((var & BRANCH_VAR_PHASE) && (br->vars & BRANCH_VAR_PHASE)) {
    for (t = 0; t < br->num_periods; t++)
      VEC_set(indices,offset+t,br->index_phase[t]);
    offset += br->num_periods;
  }
  if ((var & BRANCH_VAR_RATIO_DEV) && (br->vars & BRANCH_VAR_RATIO_DEV)) {
    for (t = 0; t < br->num_periods; t++) {
      VEC_set(indices,offset+2*t,br->index_ratio_y[t]);
      VEC_set(indices,offset+2*t+1,br->index_ratio_z[t]);
    }
    offset += 2*br->num_periods;
  }
  return indices;
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
  if ((prop & BRANCH_PROP_NOT_OUT) && BRANCH_is_on_outage(br))
    return FALSE;
  return TRUE;
}

void BRANCH_init(Branch* br, int num_periods) {  

  // Local vars
  int T;
  int t;

  // No branch
  if (!br)
    return;

  T = num_periods;
  br->num_periods = num_periods;

  br->type = BRANCH_TYPE_LINE;

  br->bus_from = NULL;
  br->bus_to = NULL;
  br->reg_bus = NULL;

  br->g = 0;
  br->g_from = 0;
  br->g_to = 0;
  br->b = 0;
  br->b_from = 0;
  br->b_to = 0;

  br->ratio_max = 1;
  br->ratio_min = 1;
  br->num_ratios = 1;

  br->phase_max = 0;
  br->phase_min = 0;

  br->P_max = 0;
  br->P_min = 0;
  br->Q_max = 0;
  br->Q_min = 0;

  br->ratingA = 0;
  br->ratingB = 0;
  br->ratingC = 0;

  br->outage = FALSE;
  br->pos_ratio_v_sens = TRUE;
  br->vars = 0x00;
  br->fixed = 0x00;
  br->bounded = 0x00;
  br->sparse = 0x00;

  br->index = 0;

  ARRAY_zalloc(br->ratio,REAL,T); 
  ARRAY_zalloc(br->phase,REAL,T);

  ARRAY_zalloc(br->index_ratio,int,T);
  ARRAY_zalloc(br->index_ratio_y,int,T);
  ARRAY_zalloc(br->index_ratio_z,int,T);
  ARRAY_zalloc(br->index_phase,int,T);

  ARRAY_zalloc(br->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(br->sens_P_l_bound,REAL,T);

  for (t = 0; t < br->num_periods; t++)
    br->ratio[t] = 1.;

  br->reg_next = NULL;
  br->from_next = NULL;
  br->to_next = NULL;
}

BOOL BRANCH_is_equal(Branch* br, Branch* other) {
  return br == other;
}

BOOL BRANCH_is_on_outage(Branch* br) {
  if (br)
    return br->outage;
  else
    return FALSE;
}

BOOL BRANCH_is_fixed_tran(Branch* br) {
  if (br)
    return br->type == BRANCH_TYPE_TRAN_FIXED;
  else
    return FALSE;
}

BOOL BRANCH_is_line(Branch* br) {
  if (br)
    return br->type == BRANCH_TYPE_LINE;
  else
    return FALSE;
}

BOOL BRANCH_is_phase_shifter(Branch* br) {
  if (br)
    return br->type == BRANCH_TYPE_TRAN_PHASE;
  else
    return FALSE;
}

BOOL BRANCH_is_tap_changer(Branch* br) {
  if (br)
    return BRANCH_is_tap_changer_v(br) || BRANCH_is_tap_changer_Q(br);
  else
    return FALSE;
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
  LIST_add(Branch,reg_br_list,reg_br,reg_next);
  return reg_br_list;
}

Branch* BRANCH_list_reg_del(Branch* reg_br_list, Branch* reg_br) {
  LIST_del(Branch,reg_br_list,reg_br,reg_next);
  return reg_br_list;
}

int BRANCH_list_reg_len(Branch* reg_br_list) {
  int len;
  LIST_len(Branch,reg_br_list,reg_next,len);
  return len;
}

Branch* BRANCH_list_from_add(Branch* from_br_list, Branch* br) {
  LIST_add(Branch,from_br_list,br,from_next);
  return from_br_list;
}

Branch* BRANCH_list_from_del(Branch* from_br_list, Branch* br) {
  LIST_del(Branch,from_br_list,br,from_next);
  return from_br_list;
}

int BRANCH_list_from_len(Branch* from_br_list) {
  int len;
  LIST_len(Branch,from_br_list,from_next,len);
  return len;
}

Branch* BRANCH_list_to_add(Branch* to_br_list, Branch* br) {
  LIST_add(Branch,to_br_list,br,to_next);
  return to_br_list;
}

Branch* BRANCH_list_to_del(Branch* to_br_list, Branch* br) {
  LIST_del(Branch,to_br_list,br,to_next);
  return to_br_list;
}

int BRANCH_list_to_len(Branch* to_br_list) {
  int len;
  LIST_len(Branch,to_br_list,to_next,len);
  return len;
}

Branch* BRANCH_new(int num_periods) {
  if (num_periods > 0) {
    Branch* branch = (Branch*)malloc(sizeof(Branch));
    BRANCH_init(branch,num_periods);
    return branch;
  }
  else
    return NULL;
}

void BRANCH_set_sens_P_u_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_P_u_bound[t] = value;
}

void BRANCH_set_sens_P_l_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_P_l_bound[t] = value;
}

void BRANCH_set_index(Branch* br, int index) {
  if (br) 
    br->index = index;
}

void BRANCH_set_type(Branch* br, int type) {
  if (br)
    br->type = type;
}

void BRANCH_set_bus_from(Branch* br, Bus* bus_from) {
  if (br)
    br->bus_from = bus_from;
}

void BRANCH_set_bus_to(Branch* br, Bus* bus_to) {
  if (br)
    br->bus_to = bus_to;
}

void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus) {
  if (br)  
    br->reg_bus = reg_bus;
}

void BRANCH_set_g(Branch* br, REAL g) {
  if (br)  
    br->g = g;
}

void BRANCH_set_g_from(Branch* br, REAL g_from) {
  if (br)
    br->g_from = g_from;
}

void BRANCH_set_g_to(Branch* br, REAL g_to) {
  if (br)
    br->g_to = g_to;
}

void BRANCH_set_b(Branch* br, REAL b) {
  if (br)
    br->b = b;
}

void BRANCH_set_b_from(Branch* br, REAL b_from) {
  if (br)  
    br->b_from = b_from;
}

void BRANCH_set_b_to(Branch* br, REAL b_to) {
  if (br)
    br->b_to = b_to;
}

void BRANCH_set_ratio(Branch* br, REAL ratio, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->ratio[t] = ratio;
}

void BRANCH_set_ratio_max(Branch* br, REAL ratio) {
  if (br)
    br->ratio_max = ratio;
}

void BRANCH_set_ratio_min(Branch* br, REAL ratio) {
  if (br)
    br->ratio_min = ratio;
}

void BRANCH_set_pos_ratio_v_sens(Branch* br, BOOL flag) {
  if (br)
    br->pos_ratio_v_sens = flag;
}

void BRANCH_set_outage(Branch* br, BOOL outage) {
  if (br)
    br->outage = outage;
}

void BRANCH_set_phase(Branch* br, REAL phase, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->phase[t] = phase;
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

  // Local vars
  int t;

  // No branch
  if (!br)
    return;

  // Time loop
  for (t = 0; t < br->num_periods; t++) {
    
    // Ratio and phase
    if (br->vars & BRANCH_VAR_RATIO)    // taps ratio
      br->ratio[t] = VEC_get(values,br->index_ratio[t]); 
    if (br->vars & BRANCH_VAR_PHASE)    // phase shift
      br->phase[t] = VEC_get(values,br->index_phase[t]);
  }
}

int BRANCH_set_flags(void* vbr, char flag_type, char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Branch* br = (Branch*)vbr;
  int t;

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
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < br->num_periods; t++)
	br->index_ratio[t] = index+t;
    }
    (*flags_ptr) |= BRANCH_VAR_RATIO;
    index += br->num_periods;
  }
  if (!((*flags_ptr) & BRANCH_VAR_PHASE) && (mask & BRANCH_VAR_PHASE)) { // phase shift
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < br->num_periods; t++)
	br->index_phase[t] = index+t;
    }
    (*flags_ptr) |= BRANCH_VAR_PHASE;
    index += br->num_periods;
  }
  if (!((*flags_ptr) & BRANCH_VAR_RATIO_DEV) && (mask & BRANCH_VAR_RATIO_DEV)) { // taps ratio deviations
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < br->num_periods; t++) {
	br->index_ratio_y[t] = index+2*t;
	br->index_ratio_z[t] = index+2*t+1;
      }
    }
    (*flags_ptr) |= BRANCH_VAR_RATIO_DEV;
    index += 2*br->num_periods;
  }
  return index;
}

void BRANCH_show(Branch* br, int t) {
  printf("branch %d\t%d\t%d\n",
	 BUS_get_number(br->bus_from),
	 BUS_get_number(br->bus_to),
	 br->type);
}

void BRANCH_propagate_data_in_time(Branch* br) {
  int t;
  if (br) {
    for (t = 1; t < br->num_periods; t++) {
      br->ratio[t] = br->ratio[0];
      br->phase[t] = br->phase[0];
    }
  }
}
