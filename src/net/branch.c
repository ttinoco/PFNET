/** @file branch.c
 *  @brief This file defines the Branch data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/branch.h>
#include <pfnet/bus.h>
#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

// Branch
struct Branch {

  // Properties
  char type;         /**< @brief %Branch type */
  char name[BRANCH_BUFFER_SIZE]; /**< @brief Branch name */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Buses
  Bus* bus_k;        /**< @brief Bus connected to the "k" side */
  Bus* bus_m;        /**< @brief Bus connected to the "m" side */
  Bus* reg_bus;      /**< @brief Bus regulated by this transformer */

  // Conductance
  REAL g;            /**< @brief Series conductance (p.u.) */
  REAL g_k;          /**< @brief %Shunt conductance on "k" side (p.u.) */
  REAL g_m;          /**< @brief %Shunt conductance on "m" side (p.u.) */

  // Susceptance
  REAL b;            /**< @brief Series susceptance (p.u.) */
  REAL b_k;          /**< @brief %Shunt susceptance on "k" side (p.u.) */
  REAL b_m;          /**< @brief %Shunt shunt susceptance "m" side (p.u.) */

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

  // Thermal ratings
  REAL ratingA;      /**< @brief Thermal rating A (p.u. system base MVA) */
  REAL ratingB;      /**< @brief Thermal rating B (p.u. system base MVA) */
  REAL ratingC;      /**< @brief Thermal rating C (p.u. system base MVA) */

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
  int* index_phase;   /**< @brief Phase shift index */

  // Sensitivities
  REAL* sens_P_u_bound;     /**< @brief Sensitivity of active power flow upper bound */
  REAL* sens_P_l_bound;     /**< @brief Sensitivity of active power flow lower bound */
  REAL* sens_ratio_u_bound; /**< @brief Sensitivity of tap ratio upper bound */
  REAL* sens_ratio_l_bound; /**< @brief Sensitivity of tap ratio lower bound */
  REAL* sens_phase_u_bound; /**< @brief Sensitivity of phase shift upper bound */
  REAL* sens_phase_l_bound; /**< @brief Sensitivity of phase shift lower bound */
  REAL* sens_i_mag_u_bound; /**< @brief Sensitivity of current magnitude upper bound */

  // Network
  Net* net; /**< @brief Network. */

  // List
  Branch* reg_next;   /**< @brief List of branches regulating a bus voltage magnitude */
  Branch* next_k;     /**< @brief List of branches connected to a bus on the "k" side */
  Branch* next_m;     /**< @brief List of branches connected to a bus in the "m" side */
};

void BRANCH_set_network(Branch* br, void* net) {
  if (br)
    br->net = (Net*)net;
}

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
      free(br->index_phase);
      free(br->sens_P_u_bound);
      free(br->sens_P_l_bound);
      free(br->sens_ratio_u_bound);
      free(br->sens_ratio_l_bound);
      free(br->sens_phase_u_bound);
      free(br->sens_phase_l_bound);
      free(br->sens_i_mag_u_bound);
      BRANCH_set_bus_k(br,NULL);
      BRANCH_set_bus_m(br,NULL);
      BRANCH_set_reg_bus(br,NULL);
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
      snprintf(br_array[i].name,(size_t)(BRANCH_BUFFER_SIZE-1),"%d",i);
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
      br->sens_ratio_u_bound[t] = 0;
      br->sens_ratio_l_bound[t] = 0;
      br->sens_phase_u_bound[t] = 0;
      br->sens_phase_l_bound[t] = 0;
      br->sens_i_mag_u_bound[t] = 0;
    }
  }
}

void BRANCH_copy_from_branch(Branch* br, Branch* other) {

  // Local variables
  int num_periods;

  // Check
  if (!br || !other)
    return;

  // Min num periods
  if (br->num_periods < other->num_periods)
    num_periods = br->num_periods;
  else
    num_periods = other->num_periods;

  // Properties
  br->type = other->type;
  strcpy(br->name,other->name);

  // Time
  // skip time

  // Buses
  // skip buses

  // Conductance
  br->g = other->g;
  br->g_k = other->g_k;
  br->g_m = other->g_m;

  // Susceptance
  br->b = other->b;
  br->b_k = other->b_k;
  br->b_m = other->b_m;

  // Tap ratio
  memcpy(br->ratio,other->ratio,num_periods*sizeof(REAL));
  br->ratio_max = other->ratio_max;
  br->ratio_min = other->ratio_min;
  br->num_ratios = other->num_ratios;
  
  // Phase shift
  memcpy(br->phase,other->phase,num_periods*sizeof(REAL));
  br->phase_max = other->phase_max;
  br->phase_min = other->phase_min;

  // Flow bounds
  br->P_max = other->P_max;
  br->P_min = other->P_min;
  br->Q_max = other->Q_max;
  br->Q_min = other->Q_min;

  // Thermal ratings
  br->ratingA = other->ratingA;
  br->ratingB = other->ratingB;
  br->ratingC = other->ratingC;
  
  // Flags
  br->outage = other->outage;
  br->pos_ratio_v_sens = other->pos_ratio_v_sens;
  br->fixed = other->fixed;
  br->bounded = other->bounded;
  br->sparse = other->sparse;
  br->vars = other->vars;
  
  // Indices
  // skip index
  memcpy(br->index_ratio,other->index_ratio,num_periods*sizeof(int));
  memcpy(br->index_phase,other->index_phase,num_periods*sizeof(int));
  
  // Sensitivities
  memcpy(br->sens_P_u_bound,other->sens_P_u_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_P_l_bound,other->sens_P_l_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_ratio_u_bound,other->sens_ratio_u_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_ratio_l_bound,other->sens_ratio_l_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_phase_u_bound,other->sens_phase_u_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_phase_l_bound,other->sens_phase_l_bound,num_periods*sizeof(REAL));
  memcpy(br->sens_i_mag_u_bound,other->sens_i_mag_u_bound,num_periods*sizeof(REAL));
  
  // List
  // skip next
}

char BRANCH_get_flags_vars(Branch* br) {
  if (br)
    return br->vars;
  else
    return 0;
}

char BRANCH_get_flags_fixed(Branch* br) {
  if (br)
    return br->fixed;
  else
    return 0;
}

char BRANCH_get_flags_bounded(Branch* br) {
  if (br)
    return br->bounded;
  else
    return 0;
}

char BRANCH_get_flags_sparse(Branch* br) {
  if (br)
    return br->sparse;
  else
    return 0;
}

char* BRANCH_get_name(Branch* br) {
  if (br)
    return br->name;
  else
    return NULL;
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

REAL* BRANCH_get_sens_P_u_bound_array(Branch* br) {
  if (br)
    return br->sens_P_u_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_P_l_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_P_l_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_P_l_bound_array(Branch* br) {
  if (br)
    return br->sens_P_l_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_ratio_u_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_ratio_u_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_ratio_u_bound_array(Branch* br) {
  if (br)
    return br->sens_ratio_u_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_ratio_l_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_ratio_l_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_ratio_l_bound_array(Branch* br) {
  if (br)
    return br->sens_ratio_l_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_phase_u_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_phase_u_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_phase_u_bound_array(Branch* br) {
  if (br)
    return br->sens_phase_u_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_phase_l_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_phase_l_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_phase_l_bound_array(Branch* br) {
  if (br)
    return br->sens_phase_l_bound;
  else
    return NULL;
}

REAL BRANCH_get_sens_i_mag_u_bound(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->sens_i_mag_u_bound[t];
  else
    return 0;
}

REAL* BRANCH_get_sens_i_mag_u_bound_array(Branch* br) {
  if (br)
    return br->sens_i_mag_u_bound;
  else
    return NULL;
}

int BRANCH_get_index(Branch* br) {
  if (br)
    return br->index;
  else
    return -1;
}

int BRANCH_get_index_ratio(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_ratio[t];
  else
    return -1;
}

int BRANCH_get_index_phase(Branch* br, int t) {
  if (br && t >= 0 && t < br->num_periods)
    return br->index_phase[t];
  else
    return -1;
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

REAL BRANCH_get_b_k(Branch* br) {
  if (br)
    return br->b_k;
  else
    return 0;
}

REAL BRANCH_get_b_m(Branch* br) {
  if (br)
    return br->b_m;
  else
    return 0;
}

REAL BRANCH_get_g(Branch* br) {
  if (br)
    return br->g;
  else
    return 0;
}

REAL BRANCH_get_g_k(Branch* br) {
  if (br)
    return br->g_k;
  else
    return 0;
}

REAL BRANCH_get_g_m(Branch* br) {
  if (br)
    return br->g_m;
  else
    return 0;
}

Bus* BRANCH_get_bus_k(Branch* br) {
  if (br)
    return br->bus_k;
  else
    return NULL;
}

Bus* BRANCH_get_bus_m(Branch* br) {
  if (br)
    return br->bus_m;
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

Branch* BRANCH_get_next_k(Branch* br) {
  if (br)
    return br->next_k;
  else
    return NULL;
}

Branch* BRANCH_get_next_m(Branch* br) {
  if (br)
    return br->next_m;
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

REAL BRANCH_get_P_max(Branch* br) {
  if (br)
    return br->P_max;
  else
    return 0;
}

REAL BRANCH_get_P_min(Branch* br) {
  if (br)
    return br->P_min;
  else
    return 0;
}
REAL BRANCH_get_Q_max(Branch* br) {
  if (br)
    return br->Q_max;
  else
    return 0;
}

REAL BRANCH_get_Q_min(Branch* br) {
  if (br)
    return br->Q_min;
  else
    return 0;
}

void BRANCH_compute_flows(Branch* br, Vec* var_values, int t, REAL* flows) {
  /** Compute the flows in this branch's pi model equivalent
   *  including the flow from the bus, the flow in the shunt elements,
   *  and the flow in the series element. These values are returned in 
   *  the 'flows' argument.
   */
    
  // Buses
  Bus* bus_k;
  Bus* bus_m;
  
  // Voltage magnitudes
  REAL v_k;
  REAL v_m;
  
  // Voltage angles
  REAL w_k;
  REAL w_m;
  
  // Phase shift
  REAL phi;
  
  // Tap ratios
  REAL a_km;
  REAL a_mk = 1.;
  
  // Series conductance and susceptance
  REAL g_km;
  REAL b_km;
  REAL g_mk;
  REAL b_mk;
  
  // Shunt conductance and susceptance
  REAL g_k_sh;
  REAL b_k_sh;
  REAL g_m_sh;
  REAL b_m_sh;

  // Intermediate values
  REAL v_k_tap_squared;
  REAL v_m_tap_squared;
  REAL v_k_v_m_tap;
  REAL theta_km;
  REAL theta_mk;

  // Check inputs
  if (!flows || !br)
    return;

  // Populate the local variables
  bus_k = BRANCH_get_bus_k(br);
  bus_m = BRANCH_get_bus_m(br);

  // Get voltage angles
  if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG) && var_values)
    w_k = VEC_get(var_values,BUS_get_index_v_ang(bus_k,t));
  else
    w_k = BUS_get_v_ang(bus_k,t);
  if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG) && var_values)
    w_m = VEC_get(var_values,BUS_get_index_v_ang(bus_m,t));
  else
    w_m = BUS_get_v_ang(bus_m,t);

  // Get voltage magnitudes
  if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VMAG) && var_values)
    v_k = VEC_get(var_values,BUS_get_index_v_mag(bus_k,t));
  else
    v_k = BUS_get_v_mag(bus_k,t);
  if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VMAG) && var_values)
    v_m = VEC_get(var_values,BUS_get_index_v_mag(bus_m,t));
  else
    v_m = BUS_get_v_mag(bus_m,t);

  // Get series conductance and susceptance
  g_km = BRANCH_get_g(br);
  b_km = BRANCH_get_b(br);
  g_mk = g_km;
  b_mk = b_km;
  g_k_sh = BRANCH_get_g_k(br);
  b_k_sh = BRANCH_get_b_k(br);
  g_m_sh = BRANCH_get_g_m(br);
  b_m_sh = BRANCH_get_b_m(br);

  // Get tap ratio from k, a_mk = 1 always
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) && var_values)
    a_km = VEC_get(var_values,BRANCH_get_index_ratio(br,t));
  else
    a_km = BRANCH_get_ratio(br,t);

  // Get phase shift angle, this should be subtracted on side "m" (to)
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) && var_values)
    phi = VEC_get(var_values,BRANCH_get_index_phase(br,t));
  else
    phi = BRANCH_get_phase(br,t);

  // Repeated calculations
  v_k_tap_squared = a_km*a_km*v_k*v_k;
  v_m_tap_squared = a_mk*a_mk*v_m*v_m;
  v_k_v_m_tap = a_km*a_mk*v_k*v_m;
  theta_km = w_k-w_m-phi;
  theta_mk = w_m-w_k+phi;

  // Calculate series elements
  // P_km_series = a_km^2*v_k^2*g_km - a_km*a_mk*v_k*v_m*( g_km*cos(w_k-w_m-phi) + b_km*sin(w_k-w_m-phi))
  flows[BRANCH_P_KM_SERIES] = (v_k_tap_squared*g_km -
			       v_k_v_m_tap*( g_km*cos(theta_km) + b_km*sin(theta_km)));
  
  // Q_km_series = -a_km^2*v_k^2*b_km - a_km*a_mk*v_k*v_m*( g_km*sin(w_k-w_m-phi) - b_km*cos(w_k-w_m-phi))
  flows[BRANCH_Q_KM_SERIES] = (-v_k_tap_squared*b_km -
			       v_k_v_m_tap*( g_km*sin(theta_km) - b_km*cos(theta_km)));
  
  // P_mk_series = a_mk^2*v_m^2*g_mk - a_mk*a_km*v_k*v_m*( g_mk*cos(w_k-w_m+phi) + b_mk*sin(w_k-w_m+phi))
  flows[BRANCH_P_MK_SERIES] = (v_m_tap_squared*g_mk -
			       v_k_v_m_tap*( g_mk*cos(theta_mk) + b_mk*sin(theta_mk)));
  
  // Q_mk_series = -a_mk^2*v_m^2*b_mk - a_mk*a_km*v_k*v_m*( g_mk*sin(w_k-w_m+phi) - b_mk*cos(w_k-w_m+phi))
  flows[BRANCH_Q_MK_SERIES] = (-v_m_tap_squared*b_mk -
			       v_k_v_m_tap*( g_mk*sin(theta_mk) - b_mk*cos(theta_mk)));

  // Calculate shunt elements
  // P_k_shunt = v_k^2*a_km^2*g_k_sh
  flows[BRANCH_P_K_SHUNT] = v_k_tap_squared*g_k_sh;

  // Q_k_shunt = -v_k^2*a_km^2*b_k_sh
  flows[BRANCH_Q_K_SHUNT] = -v_k_tap_squared*b_k_sh;

  // P_m_shunt = v_m^2*a_mk^2*g_m_sh
  flows[BRANCH_P_M_SHUNT] = v_m_tap_squared*g_m_sh;

  // Q_m_shunt = -v_m^2*a_mk^2*b_m_sh
  flows[BRANCH_Q_M_SHUNT] = -v_m_tap_squared*b_m_sh;

  // Total flows from the given bus
  // P_km = P_km_series + P_k_shunt
  flows[BRANCH_P_KM] = flows[BRANCH_P_KM_SERIES] + flows[BRANCH_P_K_SHUNT];
  
  // Q_km = Q_km_series + Q_k_shunt
  flows[BRANCH_Q_KM] = flows[BRANCH_Q_KM_SERIES] + flows[BRANCH_Q_K_SHUNT];
  
  // P_mk = P_mk_series + P_m_shunt
  flows[BRANCH_P_MK] = flows[BRANCH_P_MK_SERIES] + flows[BRANCH_P_M_SHUNT];
  
  // Q_mk = Q_mk_series + Q_m_shunt
  flows[BRANCH_Q_MK] = flows[BRANCH_Q_MK_SERIES] + flows[BRANCH_Q_M_SHUNT];
}

REAL BRANCH_get_i_km_mag(Branch* br, Vec* var_values, int t, REAL eps) {

  REAL flows[BRANCH_FLOW_SIZE];
  REAL R_km;
  REAL I_km;
  Bus* bus_k;
  REAL vk;

  if (br) {
    bus_k = BRANCH_get_bus_k(br);
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VMAG) && var_values)
      vk = VEC_get(var_values,BUS_get_index_v_mag(bus_k,t));
    else
      vk = BUS_get_v_mag(bus_k,t);
    BRANCH_compute_flows(br,var_values,t,flows);
    R_km = flows[BRANCH_P_KM]/vk;
    I_km = flows[BRANCH_Q_KM]/vk;
    return sqrt(R_km*R_km+I_km*I_km+eps);
  }
  else
    return 0;
}

sxREAL BRANCH_get_i_mk_mag(Branch* br, Vec* var_values, int t, REAL eps) {
  
  REAL flows[BRANCH_FLOW_SIZE];
  REAL R_mk;
  REAL I_mk;
  Bus* bus_m;
  REAL vm;

  if (br) {
    bus_m = BRANCH_get_bus_m(br);
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VMAG) && var_values)
      vm = VEC_get(var_values,BUS_get_index_v_mag(bus_m,t));
    else
      vm = BUS_get_v_mag(bus_m,t);
    BRANCH_compute_flows(br,var_values,t,flows);
    R_mk = flows[BRANCH_P_MK]/vm;
    I_mk = flows[BRANCH_Q_MK]/vm;
    return sqrt(R_mk*R_mk+I_mk*I_mk+eps);
  }
  else
    return 0;
}

REAL BRANCH_get_S_km_mag(Branch* br, Vec* var_values, int t) {

  REAL flows[BRANCH_FLOW_SIZE];
  REAL P_km;
  REAL Q_km;

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    P_km = flows[BRANCH_P_KM];
    Q_km = flows[BRANCH_Q_KM];
    return sqrt(P_km*P_km+Q_km*Q_km);
  }
  else
    return 0;
}

REAL BRANCH_get_S_mk_mag(Branch* br, Vec* var_values, int t) {

  REAL flows[BRANCH_FLOW_SIZE];
  REAL P_mk;
  REAL Q_mk;

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    P_mk = flows[BRANCH_P_MK];
    Q_mk = flows[BRANCH_Q_MK];
    return sqrt(P_mk*P_mk+Q_mk*Q_mk);
  }
  else
    return 0;
}

REAL BRANCH_get_P_km(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow measured at bus "k" towards bus "m".
   *  P_km = P_km_series + P_k_shunt
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_KM];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_km(Branch* br, Vec* var_values, int t) {
  /** Gets the reactive power flow measured at bus "k" towards bus "m".
   *  Q_km = Q_km_series + Q_k_shunt
   */
  
  REAL flows[BRANCH_FLOW_SIZE];
  
  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_KM];
  }
  else
    return 0;
}

REAL BRANCH_get_P_mk(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow measured at bus "m" towards bus "k".
   *  P_mk = P_mk_series + P_m_shunt
   */

  REAL flows[BRANCH_FLOW_SIZE];
  
  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_MK];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_mk(Branch* br, Vec* var_values, int t) {
  /** Gets the reactive power flow measured at bus "m" towards bus "k".
   *  Q_mk = Q_mk_series + Q_m_shunt
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_MK];
  }
  else
    return 0;
}

REAL BRANCH_get_P_km_series(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow across the series element from bus "k" to bus "m".
   *  P_km_series = a_km^2*v_k^2*g_km - a_km*a_mk*v_k*v_m*( g_km*cos(w_k-w_m-phi) + b_km*sin(w_k-w_m-phi))
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_KM_SERIES];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_km_series(Branch* br, Vec* var_values, int t) {
  /** Gets the reactive power flow across the series element from bus "k" to bus "m".
   *  Q_km_series = -a_km^2*v_k^2*b_km - a_km*a_mk*v_k*v_m*( g_km*sin(w_k-w_m-phi) - b_km*cos(w_k-w_m-phi))
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_KM_SERIES];
  }
  else
    return 0;
}

REAL BRANCH_get_P_mk_series(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow across the series element from bus "m" to bus "k".
   *  P_mk_series = -a_mk^2*v_m^2*g_mk - a_mk*a_km*v_k*v_m*( g_mk*cos(w_k-w_m+phi) + b_mk*sin(w_k-w_m+phi))
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_MK_SERIES];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_mk_series(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow across the series element from bus "m" to bus "k".
   *  Q_mk_series = -a_mk^2*v_m^2*b_mk - a_mk*a_km*v_k*v_m*( g_mk*sin(w_k-w_m+phi) - b_mk*cos(w_k-w_m+phi))
   */
  
  REAL flows[BRANCH_FLOW_SIZE];
  
  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_MK_SERIES];
  }
  else
    return 0;
}

REAL BRANCH_get_P_k_shunt(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow to the shunt element from bus "k".
   *  P_k_shunt = v_k^2*a_km^2*g_k_sh
   */
  
  REAL flows[BRANCH_FLOW_SIZE];
  
  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_K_SHUNT];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_k_shunt(Branch* br, Vec* var_values, int t) {
  /** Gets the reactive power flow to the shunt element from bus "k".
   *  Q_k_shunt = v_k^2*a_km^2*b_k_sh
   */
  
  REAL flows[BRANCH_FLOW_SIZE];
  
  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_K_SHUNT];
  }
  else
    return 0;
}

REAL BRANCH_get_P_m_shunt(Branch* br, Vec* var_values, int t) {
  /** Gets the real power flow to the shunt element from bus "m".
   *  P_m_shunt = v_m^2*a_mk^2*g_m_sh
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_P_M_SHUNT];
  }
  else
    return 0;
}

REAL BRANCH_get_Q_m_shunt(Branch* br, Vec* var_values, int t) {
  /** Gets the reactive power flow to the shunt element from bus "m".
   *  Q_m_shunt = v_m^2*a_mk^2*b_m_sh
   */
  
  REAL flows[BRANCH_FLOW_SIZE];

  if (br) {
    BRANCH_compute_flows(br,var_values,t,flows);
    return flows[BRANCH_Q_M_SHUNT];
  }
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

REAL BRANCH_get_P_km_DC(Branch* br, int t) {
  /** Active power flow (DC approx) from bus "k" to bus "m". 
   *  P_km_DC = -b (w_k - w_m - Phi_km) 
   */

  if (br && t >= 0 && t < br->num_periods) {
    return -(br->b)*(BUS_get_v_ang(br->bus_k,t)-
		     BUS_get_v_ang(br->bus_m,t)-
		     br->phase[t]);
  }
  else
    return 0;
}

// 
REAL BRANCH_get_P_mk_DC(Branch* br, int t) {
  /** Active power flow (DC approx) from bus "m" to bus "k". 
   *  P_mk_DC = -b (w_m - w_m - Phi_km)
   */

  if (br && t >= 0 && t < br->num_periods) {
    return -(br->b)*(BUS_get_v_ang(br->bus_m,t)-
		     BUS_get_v_ang(br->bus_k,t)+
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
	if (br->bounded & BRANCH_VAR_RATIO)
	  VEC_set(values,br->index_ratio[t],br->ratio_max);
	else
	  VEC_set(values,br->index_ratio[t],BRANCH_INF_RATIO);
	break;

      case LOWER_LIMITS:
	if (br->bounded & BRANCH_VAR_RATIO)
	  VEC_set(values,br->index_ratio[t],br->ratio_min);
	else
	  VEC_set(values,br->index_ratio[t],-BRANCH_INF_RATIO);
	break;

      default:
	VEC_set(values,br->index_ratio[t],br->ratio[t]);
      }
    }
    if (br->vars & BRANCH_VAR_PHASE) { // phase shift
      switch(code) {

      case UPPER_LIMITS:
	if (br->bounded & BRANCH_VAR_PHASE)
	  VEC_set(values,br->index_phase[t],br->phase_max);
	else
	  VEC_set(values,br->index_phase[t],BRANCH_INF_PHASE);
	break;

      case LOWER_LIMITS:
	if (br->bounded & BRANCH_VAR_PHASE)
	  VEC_set(values,br->index_phase[t],br->phase_min);
	else
	  VEC_set(values,br->index_phase[t],-BRANCH_INF_PHASE);
	break;

      default:
	VEC_set(values,br->index_phase[t],br->phase[t]);
      }
    }
  }
}

char* BRANCH_get_var_info_string(Branch* br, int index) {

  // Local variables
  char* info;

  //Check
  if (!br)
    return NULL;

  // Taps ratio
  if ((br->vars & BRANCH_VAR_RATIO) &&
      index >= br->index_ratio[0] &&
      index <= br->index_ratio[br->num_periods-1]) {
    info = (char*)malloc(BRANCH_BUFFER_SIZE*sizeof(char));
    snprintf(info,BRANCH_BUFFER_SIZE*sizeof(char),
	     "branch:%d:tap ratio:%d",br->index,index-br->index_ratio[0]);
    return info;
  }

  // Phase shift
  if ((br->vars & BRANCH_VAR_PHASE) &&
      index >= br->index_phase[0] &&
      index <= br->index_phase[br->num_periods-1]) {
    info = (char*)malloc(BRANCH_BUFFER_SIZE*sizeof(char));
    snprintf(info,BRANCH_BUFFER_SIZE*sizeof(char),
	     "branch:%d:phase shift:%d",br->index,index-br->index_phase[0]);
    return info;
  }

  // Return
  return NULL;
}

int BRANCH_get_num_vars(void* vbr, unsigned char var, int t_start, int t_end) {

  // Local vars
  Branch* br = (Branch*)vbr;
  int num_vars = 0;
  int dt;

  // Checks
  if (!br)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > br->num_periods-1)
    t_end = br->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & BRANCH_VAR_RATIO) && (br->vars & BRANCH_VAR_RATIO)) // taps ratio
    num_vars += dt;
  if ((var & BRANCH_VAR_PHASE) && (br->vars & BRANCH_VAR_PHASE)) // phase shifts
    num_vars += dt;
  return num_vars;
}

Vec* BRANCH_get_var_indices(void* vbr, unsigned char var, int t_start, int t_end) {

  // Local vars
  Branch* br = (Branch*)vbr;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!br)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > br->num_periods-1)
    t_end = br->num_periods-1;

  // Indices
  indices = VEC_new(BRANCH_get_num_vars(vbr,var,t_start,t_end));
  if ((var & BRANCH_VAR_RATIO) && (br->vars & BRANCH_VAR_RATIO)) { // taps ratio
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,br->index_ratio[t]);
      offset++;
    }
  }
  if ((var & BRANCH_VAR_PHASE) && (br->vars & BRANCH_VAR_PHASE)) { // phase shift
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,br->index_phase[t]);
      offset++;
    }
  }
  return indices;
}

char* BRANCH_get_json_string(Branch* branch, char* output) {

  // Local variables
  char temp[BRANCH_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No branch
  if (!branch)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*BRANCH_BUFFER_SIZE*BRANCH_NUM_JSON_FIELDS*branch->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",branch->index,FALSE);
  JSON_int(temp,output,"type",branch->type,FALSE);
  JSON_int(temp,output,"num_periods",branch->num_periods,FALSE);
  JSON_str(temp,output,"name",branch->name,FALSE);
  JSON_obj(temp,output,"bus_k",branch->bus_k,BUS_get_index,FALSE);
  JSON_obj(temp,output,"bus_m",branch->bus_m,BUS_get_index,FALSE);
  JSON_obj(temp,output,"reg_bus",branch->reg_bus,BUS_get_index,FALSE);
  JSON_float(temp,output,"g",branch->g,FALSE);
  JSON_float(temp,output,"g_k",branch->g_k,FALSE);
  JSON_float(temp,output,"g_m",branch->g_m,FALSE);
  JSON_float(temp,output,"b",branch->b,FALSE);
  JSON_float(temp,output,"b_k",branch->b_k,FALSE);
  JSON_float(temp,output,"b_m",branch->b_m,FALSE);
  JSON_array_float(temp,output,"ratio",branch->ratio,branch->num_periods,FALSE);
  JSON_float(temp,output,"ratio_max",branch->ratio_max,FALSE);
  JSON_float(temp,output,"ratio_min",branch->ratio_min,FALSE);
  JSON_array_float(temp,output,"phase",branch->phase,branch->num_periods,FALSE);
  JSON_float(temp,output,"phase_max",branch->phase_max,FALSE);
  JSON_float(temp,output,"phase_min",branch->phase_min,FALSE);
  JSON_float(temp,output,"ratingA",branch->ratingA,FALSE);
  JSON_float(temp,output,"ratingB",branch->ratingB,FALSE);
  JSON_float(temp,output,"ratingC",branch->ratingC,FALSE);
  JSON_bool(temp,output,"outage",branch->outage,FALSE);
  JSON_bool(temp,output,"pos_ratio_v_sens",branch->pos_ratio_v_sens,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL BRANCH_has_pos_ratio_v_sens(Branch* branch) {
  if (branch)
    return branch->pos_ratio_v_sens;
  else
    return FALSE;
}

BOOL BRANCH_has_flags(void* vbr, char flag_type, unsigned char mask) {
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
  ARRAY_clear(br->name,char,BRANCH_BUFFER_SIZE);

  br->bus_k = NULL;
  br->bus_m = NULL;
  br->reg_bus = NULL;

  br->g = 0;
  br->g_k = 0;
  br->g_m = 0;
  br->b = 0;
  br->b_k = 0;
  br->b_m = 0;

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

  br->index = -1;

  ARRAY_zalloc(br->ratio,REAL,T);
  ARRAY_zalloc(br->phase,REAL,T);

  ARRAY_zalloc(br->index_ratio,int,T);
  ARRAY_zalloc(br->index_phase,int,T);

  ARRAY_zalloc(br->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(br->sens_P_l_bound,REAL,T);
  ARRAY_zalloc(br->sens_ratio_u_bound,REAL,T);
  ARRAY_zalloc(br->sens_ratio_l_bound,REAL,T);
  ARRAY_zalloc(br->sens_phase_u_bound,REAL,T);
  ARRAY_zalloc(br->sens_phase_l_bound,REAL,T);
  ARRAY_zalloc(br->sens_i_mag_u_bound,REAL,T);

  for (t = 0; t < br->num_periods; t++)
    br->ratio[t] = 1.;

  br->net = NULL;

  br->reg_next = NULL;
  br->next_k = NULL;
  br->next_m = NULL;
};

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
    return br->type == BRANCH_TYPE_TRAN_TAP_V;
  else
    return FALSE;
}

BOOL BRANCH_is_tap_changer_Q(Branch* br) {
  if (br)
    return br->type == BRANCH_TYPE_TRAN_TAP_Q;
  else
    return FALSE;
}

BOOL BRANCH_is_part_of_3_winding_transformer(Branch* br) {
  if (br)
    return BUS_is_star(br->bus_k) || BUS_is_star(br->bus_m);
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

Branch* BRANCH_list_k_add(Branch* k_br_list, Branch* br) {
  LIST_add(Branch,k_br_list,br,next_k);
  return k_br_list;
}

Branch* BRANCH_list_k_del(Branch* k_br_list, Branch* br) {
  LIST_del(Branch,k_br_list,br,next_k);
  return k_br_list;
}

int BRANCH_list_k_len(Branch* k_br_list) {
  int len;
  LIST_len(Branch,k_br_list,next_k,len);
  return len;
}

Branch* BRANCH_list_m_add(Branch* m_br_list, Branch* br) {
  LIST_add(Branch,m_br_list,br,next_m);
  return m_br_list;
}

Branch* BRANCH_list_m_del(Branch* m_br_list, Branch* br) {
  LIST_del(Branch,m_br_list,br,next_m);
  return m_br_list;
}

int BRANCH_list_m_len(Branch* m_br_list) {
  int len;
  LIST_len(Branch,m_br_list,next_m,len);
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

void BRANCH_set_name(Branch* br, char* name) {
  if (br)
    strncpy(br->name,name,(size_t)(BRANCH_BUFFER_SIZE-1));
}

void BRANCH_set_sens_P_u_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_P_u_bound[t] = value;
}

void BRANCH_set_sens_P_l_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_P_l_bound[t] = value;
}

void BRANCH_set_sens_ratio_u_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_ratio_u_bound[t] = value;
}

void BRANCH_set_sens_ratio_l_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_ratio_l_bound[t] = value;
}

void BRANCH_set_sens_phase_u_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_phase_u_bound[t] = value;
}

void BRANCH_set_sens_phase_l_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_phase_l_bound[t] = value;
}

void BRANCH_set_sens_i_mag_u_bound(Branch* br, REAL value, int t) {
  if (br && t >= 0 && t < br->num_periods)
    br->sens_i_mag_u_bound[t] = value;
}

void BRANCH_set_index(Branch* br, int index) {
  if (br)
    br->index = index;
}

void BRANCH_set_type(Branch* br, int type) {
  if (br)
    br->type = type;
}

void BRANCH_set_bus_k(Branch* br, Bus* bus_k) {
  Bus* old_bus_k;
  if (br) {
    old_bus_k = br->bus_k;
    br->bus_k = NULL;
    BUS_del_branch_k(old_bus_k,br);
    br->bus_k = bus_k;
    BUS_add_branch_k(br->bus_k,br);
  }
}

void BRANCH_set_bus_m(Branch* br, Bus* bus_m) {
  Bus* old_bus_m;
  if (br) {
    old_bus_m = br->bus_m;
    br->bus_m = NULL;
    BUS_del_branch_m(old_bus_m,br);
    br->bus_m = bus_m;
    BUS_add_branch_m(br->bus_m,br);
  }
}

void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus) {
  Bus* old_reg_bus;
  if (br) {
    old_reg_bus = br->reg_bus;
    br->reg_bus = NULL;
    BUS_del_reg_tran(old_reg_bus,br);
    br->reg_bus = reg_bus;
    BUS_add_reg_tran(br->reg_bus,br);
    if (reg_bus)
      br->type = BRANCH_TYPE_TRAN_TAP_V;
    else if (br->type == BRANCH_TYPE_TRAN_TAP_V)
      br->type = BRANCH_TYPE_TRAN_FIXED;
  }
}

void BRANCH_set_g(Branch* br, REAL g) {
  if (br)
    br->g = g;
}

void BRANCH_set_g_k(Branch* br, REAL g_k) {
  if (br)
    br->g_k = g_k;
}

void BRANCH_set_g_m(Branch* br, REAL g_m) {
  if (br)
    br->g_m = g_m;
}

void BRANCH_set_b(Branch* br, REAL b) {
  if (br)
    br->b = b;
}

void BRANCH_set_b_k(Branch* br, REAL b_k) {
  if (br)
    br->b_k = b_k;
}

void BRANCH_set_b_m(Branch* br, REAL b_m) {
  if (br)
    br->b_m = b_m;
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
    if (br->outage != outage)
      NET_inc_state_tag(br->net);
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

int BRANCH_set_flags(void* vbr, char flag_type, unsigned char mask, int index) {

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
  return index;
}

void BRANCH_show(Branch* br, int t) {
  printf("branch %d\t%d\t%d\n",
	 BUS_get_number(br->bus_k),
	 BUS_get_number(br->bus_m),
	 br->type);
}

void BRANCH_propagate_data_in_time(Branch* br, int start, int end) {
  int t;
  if (br) {
    if (start < 0)
      start = 0;
    if (end > br->num_periods)
      end = br->num_periods;
    for (t = start+1; t < end; t++) {
      br->ratio[t] = br->ratio[start];
      br->phase[t] = br->phase[start];
    }
  }
}

void BRANCH_power_flow_count(Branch* br, int* J_nnz, int* H_nnz, int t, BOOL km) {

  // Local variables
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int k;
  int m;

  // Check
  if (!br || !J_nnz || !H_nnz)
    return;

  // Outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
  }

  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);

  // Direction
  if (km) {
    k = 0;
    m = 1;
  }
  else {
    k = 1;
    m = 0;
  }

  //***********
  if (var_w[k]) { // wk var

    // J
    (*J_nnz)++; // Pkm/dwk
    (*J_nnz)++; // Qkm/dwk
    
    // H
    (*H_nnz)++;   // wk and wk
    if (var_v[k]) 
      (*H_nnz)++; // wk and vk
    if (var_w[m]) 
      (*H_nnz)++; // wk and wm
    if (var_v[m]) 
      (*H_nnz)++; // wk and vm
    if (var_a)    
      (*H_nnz)++; // wk and a
    if (var_phi)  
      (*H_nnz)++; // wk and phi
  }
  
  //**********
  if (var_v[k]) { // vk var
    
    // J
    (*J_nnz)++; // Pkm/dvk
    (*J_nnz)++; // Qkm/dvk
    
    // H
    (*H_nnz)++;   // vk and vk
    if (var_w[m]) 
      (*H_nnz)++; // vk and wm
    if (var_v[m]) 
      (*H_nnz)++; // vk and vm
    if (var_a)    
      (*H_nnz)++; // vk and a
    if (var_phi)  
      (*H_nnz)++; // vk and phi
  }
  
  //***********
  if (var_w[m]) { // wm var
    
    // J
    (*J_nnz)++; // Pkm/dwm
    (*J_nnz)++; // Qkm/dwm
    
    // H
    (*H_nnz)++;   // wm and wm
    if (var_v[m])
      (*H_nnz)++; // wm and vm
    if (var_a)     
      (*H_nnz)++; // wm and a
    if (var_phi)    
      (*H_nnz)++; // wm and phi
  }
  
  //***********
  if (var_v[m]) { // vm var
    
    // J
    (*J_nnz)++; // Pkm/dvm
    (*J_nnz)++; // Qkm/dvm
    
    // H
    (*H_nnz)++;   // vm and vm
    if (var_a)   
      (*H_nnz)++; // vm and a
    if (var_phi) 
      (*H_nnz)++; // vm and phi
  }
  
  //********
  if (var_a) { // a var
    
    // J
    (*J_nnz)++; // Pkm/da
    (*J_nnz)++; // Qkm/da
    
    // H
    if (k == 0)  
      (*H_nnz)++; // a and a (important check k==0) 
    if (var_phi) 
      (*H_nnz)++; // a and phi
  }
  
  //**********
  if (var_phi) { // phi var
    
    // J
    (*J_nnz)++; // Pkm/dphi
    (*J_nnz)++; // Qkm/dphi
    
    // H
    (*H_nnz)++; // phi and phi
  }  
}

void BRANCH_power_flow_analyze(Branch* br, int* J_nnz, Mat* J, int index_P, int index_Q, int* H_nnz, Mat* H, int t, BOOL km) {

  // Local variables
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int w_index[2];
  int v_index[2];
  int a_index;
  int phi_index;
  int k;
  int m;

  // Check
  if (!br || !J_nnz || !J || !H_nnz || !H) 
    return;

  // Outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    w_index[k] = BUS_get_index_v_ang(bus[k],t);
    v_index[k] = BUS_get_index_v_mag(bus[k],t);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
  }

  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);
  a_index = BRANCH_get_index_ratio(br,t);
  phi_index = BRANCH_get_index_phase(br,t);

  // Direction
  if (km) {
    k = 0;
    m = 1;
  }
  else {
    k = 1;
    m = 0;
  }
  
  //***********
  if (var_w[k]) { // wk var
    
    // J
    MAT_set_i(J,*J_nnz,index_P);
    MAT_set_j(J,*J_nnz,w_index[k]);
    (*J_nnz)++; // Pkm/dwk
    
    MAT_set_i(J,*J_nnz,index_Q);
    MAT_set_j(J,*J_nnz,w_index[k]);
    (*J_nnz)++; // Qkm/dwk
    
    // H
    MAT_set_i(J,*H_nnz,w_index[k]);
    MAT_set_j(J,*H_nnz,w_index[k]);
    (*H_nnz)++;    // wk and wk
    if (var_v[k]) { 
      MAT_set_i(H[k],*H_nnz,w_index[k]);
      MAT_set_j(H[k],*H_nnz,v_index[k]);
      (*H_nnz)++; // wk and vk
    }
    if (var_w[m]) { 
      MAT_set_i(H[k],H_nnz_val,w_index[k]);
      MAT_set_j(H[k],H_nnz_val,w_index[m]);
      (*H_nnz)++; // wk and wm
    }
    if (var_v[m]) { 
      MAT_set_i(H[k],H_nnz_val,w_index[k]);
      MAT_set_j(H[k],H_nnz_val,v_index[m]);
      (*H_nnz)++; // wk and vm
    }
    if (var_a) {  
      MAT_set_i(H[k],H_nnz_val,w_index[k]);
      MAT_set_j(H[k],H_nnz_val,a_index);
      (*H_nnz)++; // wk and a
    }
    if (var_phi) { 
      MAT_set_i(H[k],H_nnz_val,w_index[k]);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      (*H_nnz)++; // wk and phi
    }
  }
  
  //***********
  if (var_v[k]) { // vk var
    
    // J
    MAT_set_i(J,*J_nnz,P_index[m]); // dPm/dvk
    MAT_set_j(J,*J_nnz,v_index[k]);
    (*J_nnz)++;
    
    MAT_set_i(J,*J_nnz,Q_index[m]); // dQm/dvk
    MAT_set_j(J,*J_nnz,v_index[k]);
    (*J_nnz)++;
    
    // H
    if (var_w[m]) { // vk and wm
      MAT_set_i(H[k],H_nnz_val,v_index[k]);
      MAT_set_j(H[k],H_nnz_val,w_index[m]);
      (*H_nnz)++;
    }
    if (var_v[m]) { // vk and vm
      MAT_set_i(H[k],H_nnz_val,v_index[k]);
      MAT_set_j(H[k],H_nnz_val,v_index[m]);
      (*H_nnz)++;
    }
    if (var_a) {   // vk and a
      MAT_set_i(H[k],H_nnz_val,v_index[k]);
      MAT_set_j(H[k],H_nnz_val,a_index);
      (*H_nnz)++;
    }
    if (var_phi) { // vk and phi
      MAT_set_i(H[k],H_nnz_val,v_index[k]);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      (*H_nnz)++;
    }
  }

  //***********
  if (var_w[m]) { // wm var
    
    // J
    // Nothing
    
    // H
    MAT_set_i(H[k],H_nnz_val,w_index[m]); // wm and wm
    MAT_set_j(H[k],H_nnz_val,w_index[m]);
    (*H_nnz)++;
    if (var_v[m]) {   // wm and vm
      MAT_set_i(H[k],H_nnz_val,w_index[m]);
      MAT_set_j(H[k],H_nnz_val,v_index[m]);
      (*H_nnz)++;
    }
    if (var_a) {      // wm and a
      MAT_set_i(H[k],H_nnz_val,w_index[m]);
      MAT_set_j(H[k],H_nnz_val,a_index);
      (*H_nnz)++;
    }
    if (var_phi) {    // wm and phi
      MAT_set_i(H[k],H_nnz_val,w_index[m]);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      (*H_nnz)++;
    }
  }
  
  //***********
  if (var_v[m]) { // vm var
    
    // J
    // Nothing
    
    // H
    if (var_a) {   // vm and a
      MAT_set_i(H[k],H_nnz_val,v_index[m]);
      MAT_set_j(H[k],H_nnz_val,a_index);
      (*H_nnz)++;
    }
    if (var_phi) { // vm and phi
      MAT_set_i(H[k],H_nnz_val,v_index[m]);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      (*H_nnz)++;
    }
  }
  
  //********
  if (var_a) { // a var
    
    // J
    MAT_set_i(J,*J_nnz,P_index[k]); // dPk/da
    MAT_set_j(J,*J_nnz,a_index);
    (*J_nnz)++;
    
    MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/da
    MAT_set_j(J,*J_nnz,a_index);
    (*J_nnz)++;
    
    // H
    if (k == 0) { // a and a (important check k==0)
      MAT_set_i(H[k],H_nnz_val,a_index);
      MAT_set_j(H[k],H_nnz_val,a_index);
      (*H_nnz)++;
    }
    if (var_phi) { // a and phi
      MAT_set_i(H[k],H_nnz_val,a_index);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      (*H_nnz)++;
    }
  }
  
  //**********
  if (var_phi) { // phi var
    
    // J
    MAT_set_i(J,*J_nnz,P_index[k]); // dPk/dphi
    MAT_set_j(J,*J_nnz,phi_index);
    (*J_nnz)++;
    
    MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/dphi
    MAT_set_j(J,*J_nnz,phi_index);
    (*J_nnz)++;
    
    // H
    MAT_set_i(H[k],H_nnz_val,phi_index);
    MAT_set_j(H[k],H_nnz_val,phi_index);
    (*H_nnz)++; // phi and phi
  }
  
}

void BRANCH_power_flow_eval(Branch* br, REAL* P, REAL* Q, int* J_nnz, REAL* JP, REAL* JQ, int* H_nnz, REAL* HP, REAL* HQ, Vec* x, int t, BOOL km) {

}
