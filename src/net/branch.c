/** @file branch.c
 *  @brief This file defines the Branch data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
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
  int* index_phase;   /**< @brief Phase shift index */

  // Sensitivities
  REAL* sens_P_u_bound;  /**< @brief Sensitivity of active power flow upper bound */
  REAL* sens_P_l_bound;  /**< @brief Sensitivity of active power flow lower bound */

  // List
  Branch* reg_next;   /**< @brief List of branches regulating a bus voltage magnitude */
  Branch* next_k;     /**< @brief List of branches connected to a bus on the "k" side */
  Branch* next_m;     /**< @brief List of branches connected to a bus in the "m" side */
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

void BRANCH_compute_flows(Branch* br, Vec* var_values, int t, REAL* flows) {
  /** Compute the flows in this branch's pi model equivalent
   *  including the flow from the bus, the flow in the shunt elements,
   *  and the flow in the series element. These values are returned in 
   *  the 'flows' argument.
   */

  int i;
  
  // Buses
  Bus* bus_k;
  Bus* bus_m;
  
  // Voltages
  REAL v_k;
  REAL v_m;
  
  // Angles
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
  if (!flows)
    return;
  else if (!br || t < 0) {
    for (i=0; i < BRANCH_FLOW_SIZE; i++)
      flows[i] = 0;
    return;
  }

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

REAL BRANCH_get_P_km(Branch* br, Vec* var_values, int t) {
  /** Get the real power flow measured at bus "k" towards bus "m".
   *  P_km = P_km_series + P_k_shunt
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_P_KM];
}

REAL BRANCH_get_Q_km(Branch* br, Vec* var_values, int t) {
  /** Get the reactive power flow measured at bus "k" towards bus "m".
   *  Q_km = Q_km_series + Q_k_shunt
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_Q_KM];
}

REAL BRANCH_get_P_mk(Branch* br, Vec* var_values, int t) {
  /** Get the real power flow measured at bus "m" towards bus "k".
   *  P_mk = P_mk_series + P_m_shunt
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_P_MK];
}

REAL BRANCH_get_Q_mk(Branch* br, Vec* var_values, int t) {
  /** Get the reactive power flow measured at bus "m" towards bus "k".
   *  Q_mk = Q_mk_series + Q_m_shunt
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_Q_MK];
}

REAL BRANCH_get_P_km_series(Branch* br, Vec* var_values, int t) {
  /** Get the real power flow across the series element from bus "k" to bus "m".
   *  P_km_series = a_km^2*v_k^2*g_km - a_km*a_mk*v_k*v_m*( g_km*cos(w_k-w_m-phi) + b_km*sin(w_k-w_m-phi))
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_P_KM_SERIES];
}

REAL BRANCH_get_Q_km_series(Branch* br, Vec* var_values, int t) {
  /** Get the reactive power flow across the series element from bus "k" to bus "m".
   *  Q_km_series = -a_km^2*v_k^2*b_km - a_km*a_mk*v_k*v_m*( g_km*sin(w_k-w_m-phi) - b_km*cos(w_k-w_m-phi))
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_Q_KM_SERIES];
}

REAL BRANCH_get_P_mk_series(Branch* br, Vec* var_values, int t) {
  /** Get the real power flow across the series element from bus "m" to bus "k".
   *  P_mk_series = -a_mk^2*v_m^2*g_mk - a_mk*a_km*v_k*v_m*( g_mk*cos(w_k-w_m+phi) + b_mk*sin(w_k-w_m+phi))
   */
  REAL flows[BRANCH_FLOW_SIZE];
  BRANCH_compute_flows(br,var_values,t,flows);
  return flows[BRANCH_P_MK_SERIES];
}

REAL BRANCH_get_Q_mk_series(Branch* br, Vec* var_values, int t) {
   /** Get the real power flow across the series element from bus "m" to bus "k".
    *  Q_mk_series = -a_mk^2*v_m^2*b_mk - a_mk*a_km*v_k*v_m*( g_mk*sin(w_k-w_m+phi) - b_mk*cos(w_k-w_m+phi))
    */
   REAL flows[BRANCH_FLOW_SIZE];
   BRANCH_compute_flows(br,var_values,t,flows);
   return flows[BRANCH_Q_MK_SERIES];
}

REAL BRANCH_get_P_k_shunt(Branch* br, Vec* var_values, int t) {
   /** Get the real power flow to the shunt element from bus "k".
    *  P_k_shunt = v_k^2*a_km^2*g_k_sh
    */
   REAL flows[BRANCH_FLOW_SIZE];
   BRANCH_compute_flows(br,var_values,t,flows);
   return flows[BRANCH_P_K_SHUNT];
}

REAL BRANCH_get_Q_k_shunt(Branch* br, Vec* var_values, int t) {
   /** Get the reactive power flow to the shunt element from bus "k".
    *  Q_k_shunt = v_k^2*a_km^2*b_k_sh
    */
   REAL flows[BRANCH_FLOW_SIZE];
   BRANCH_compute_flows(br,var_values,t,flows);
   return flows[BRANCH_Q_K_SHUNT];
}

REAL BRANCH_get_P_m_shunt(Branch* br, Vec* var_values, int t) {
   /** Get the real power flow to the shunt element from bus "m".
    *  P_m_shunt = v_m^2*a_mk^2*g_m_sh
    */
   REAL flows[BRANCH_FLOW_SIZE];
   BRANCH_compute_flows(br,var_values,t,flows);
   return flows[BRANCH_P_M_SHUNT];
}

REAL BRANCH_get_Q_m_shunt(Branch* br, Vec* var_values, int t) {
   /** Get the reactive power flow to the shunt element from bus "m".
    *  Q_m_shunt = v_m^2*a_mk^2*b_m_sh
    */
   REAL flows[BRANCH_FLOW_SIZE];
   BRANCH_compute_flows(br,var_values,t,flows);
   return flows[BRANCH_Q_M_SHUNT];
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
  BOOL resize;
  int i;

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
  
  // Start
  strcpy(output,"{ ");

  // Type
  sprintf(temp,"\"type\" : %d", branch->type);
  strcat(output,temp);
  strcat(output,", ");

  // Num periods
  sprintf(temp,"\"num_periods\" : %d", branch->num_periods);
  strcat(output,temp);
  strcat(output,", ");
  
  // Bus k
  if (branch->bus_k)
    sprintf(temp,"\"bus_k\" : %d", BUS_get_index(branch->bus_k));
  else
    sprintf(temp,"\"bus_k\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // Bus m
  if (branch->bus_m)
    sprintf(temp,"\"bus_m\" : %d", BUS_get_index(branch->bus_m));
  else
    sprintf(temp,"\"bus_m\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // Reg bus
  if (branch->reg_bus)
    sprintf(temp,"\"reg_bus\" : %d", BUS_get_index(branch->reg_bus));
  else
    sprintf(temp,"\"reg_bus\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // g
  sprintf(temp,"\"g\" : %.10e", branch->g);
  strcat(output,temp);
  strcat(output,", ");

  // g k
  sprintf(temp,"\"g_k\" : %.10e", branch->g_k);
  strcat(output,temp);
  strcat(output,", ");

  // g m
  sprintf(temp,"\"g_m\" : %.10e", branch->g_m);
  strcat(output,temp);
  strcat(output,", ");

  // b
  sprintf(temp,"\"b\" : %.10e", branch->b);
  strcat(output,temp);
  strcat(output,", ");

  // b k
  sprintf(temp,"\"b_k\" : %.10e", branch->b_k);
  strcat(output,temp);
  strcat(output,", ");

  // b m
  sprintf(temp,"\"b_m\" : %.10e", branch->b_m);
  strcat(output,temp);
  strcat(output,", ");

  // Ratio
  strcat(output,"\"ratio\" : [ ");
  for (i = 0; i < branch->num_periods; i++) {
    sprintf(temp,"%.10e", branch->ratio[i]);
    strcat(output,temp);
    if (i < branch->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  
  // Ratio max
  sprintf(temp,"\"ratio_max\" : %.10e", branch->ratio_max);
  strcat(output,temp);
  strcat(output,", ");

  // Ratio min
  sprintf(temp,"\"ratio_min\" : %.10e", branch->ratio_min);
  strcat(output,temp);
  strcat(output,", ");

  // Num ratios
  sprintf(temp,"\"num_ratios\" : %d", branch->num_ratios);
  strcat(output,temp);
  strcat(output,", ");

  // Phase
  strcat(output,"\"phase\" : [ ");
  for (i = 0; i < branch->num_periods; i++) {
    sprintf(temp,"%.10e", branch->phase[i]);
    strcat(output,temp);
    if (i < branch->num_periods-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
 
  // Phase max
  sprintf(temp,"\"phase_max\" : %.10e", branch->phase_max);
  strcat(output,temp);
  strcat(output,", ");

  // Phase min
  sprintf(temp,"\"phase_min\" : %.10e", branch->phase_min);
  strcat(output,temp);
  strcat(output,", ");

  // Rating A
  sprintf(temp,"\"ratingA\" : %.10e", branch->ratingA);
  strcat(output,temp);
  strcat(output,", ");
  
  // Rating B
  sprintf(temp,"\"ratingB\" : %.10e", branch->ratingB);
  strcat(output,temp);
  strcat(output,", ");
  
  // Rating C
  sprintf(temp,"\"ratingC\" : %.10e", branch->ratingC);
  strcat(output,temp);
  strcat(output,", ");

  // Outage
  sprintf(temp,"\"outage\" : %s", branch->outage ? "true" : "false");
  strcat(output,temp);
  strcat(output,", ");

  // Pos ratio-voltage sensitivity
  sprintf(temp,"\"pos_ratio_v_sens\" : %s", branch->pos_ratio_v_sens ? "true" : "false");
  strcat(output,temp);
  strcat(output,", ");

  // Index
  sprintf(temp,"\"index\" : %d", branch->index);
  strcat(output,temp);
  strcat(output,", ");
  
  // Reg next
  if (branch->reg_next)
    sprintf(temp,"\"reg_next\" : %d", BRANCH_get_index(branch->reg_next));
  else
    sprintf(temp,"\"reg_next\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // Next k
  if (branch->next_k)
    sprintf(temp,"\"next_k\" : %d", BRANCH_get_index(branch->next_k));
  else
    sprintf(temp,"\"next_k\" : %s", "null");
  strcat(output,temp);
  strcat(output,", ");

  // Next m
  if (branch->next_m)
    sprintf(temp,"\"next_m\" : %d", BRANCH_get_index(branch->next_m));
  else
    sprintf(temp,"\"next_m\" : %s", "null");
  strcat(output,temp);
  strcat(output,"");

  // End
  strcat(output," }");
  
  // Resiye
  if (resize)
    output = (char*)realloc(output,sizeof(char)*(strlen(output)+1)); // +1 important!

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

  br->index = 0;

  ARRAY_zalloc(br->ratio,REAL,T);
  ARRAY_zalloc(br->phase,REAL,T);

  ARRAY_zalloc(br->index_ratio,int,T);
  ARRAY_zalloc(br->index_phase,int,T);

  ARRAY_zalloc(br->sens_P_u_bound,REAL,T);
  ARRAY_zalloc(br->sens_P_l_bound,REAL,T);

  for (t = 0; t < br->num_periods; t++)
    br->ratio[t] = 1.;

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

void BRANCH_set_bus_k(Branch* br, Bus* bus_k) {
  if (br)
    br->bus_k = bus_k;
}

void BRANCH_set_bus_m(Branch* br, Bus* bus_m) {
  if (br)
    br->bus_m = bus_m;
}

void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus) {
  if (br)
    br->reg_bus = reg_bus;
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

void BRANCH_propagate_data_in_time(Branch* br) {
  int t;
  if (br) {
    for (t = 1; t < br->num_periods; t++) {
      br->ratio[t] = br->ratio[0];
      br->phase[t] = br->phase[0];
    }
  }
}
