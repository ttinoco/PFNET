/** @file facts.c
 *  @brief This file defines the Facts data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/facts.h>
#include <pfnet/bus.h>
#include <pfnet/array.h>
#include <pfnet/net.h>
#include <pfnet/json_macros.h>

struct Facts {

  // Bus
  Bus* bus_k;          /**< @brief Bus connected to "k" side */
  Bus* bus_m;          /**< @brief Bus connected to "m" side */
  Bus* reg_bus;        /**< @brief Bus regulated by this FACTS device */
  
  // Times
  int num_periods;     /**< @brief Number of time periods. */

  // Properties
  char name[FACTS_BUFFER_SIZE]; /**< @brief Facts name */
  
  // Flags
  short int pre_cont_status;   /**< @brief Flag for indicating whether the FACTS was in service before applying the contingency */
  BOOL in_service;     /**< @brief Flag for indicating whether FACTS is in service */
  char fixed;          /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;        /**< @brief Flags for indicating which quantities should be bounded */
  char vars;           /**< @brief Flags for indicating which quantities should be treated as variables */
  char sparse;         /**< @brief Flags for indicating which control adjustments should be sparse */

  // Series mode
  char mode_s;         /**< @brief Series link mode (bypass, constant impedance, constant voltage) */
  
  // Interface
  REAL* P_k;           /**< @brief Active power injected into the "k" bus (p.u.) */
  REAL* Q_k;           /**< @brief Reactive power injected into the "k" bus (p.u.) */ 
  REAL* P_m;           /**< @brief Active power injected into the "m" bus (p.u.) */
  REAL* Q_m;           /**< @brief Reactive power injected into the "m" bus (p.u.) */
  
  // Converters
  REAL* Q_sh;          /**< @brief Reactive power provided by shunt converter (p.u.) */
  REAL* Q_s;           /**< @brief Reactive power provided by series converter (p.u.) */
  REAL* P_dc;          /**< @brief DC power exchanged from shunt to series converter (p.u.) */
  REAL Q_par;          /**< @brief Reactive power participation factor of shunt converter for regulating bus voltage magnitude (unitless) */
  REAL rmpct;          /**< @brief Plant Reactive power participation factor of shunt converter for regulating bus voltage magnitude (percent) */

  // Set points 
  REAL* P_set;          /**< @brief Active power set-point at the "m" bus (p.u.) */
  REAL* Q_set;          /**< @brief Reactive power set-point at the "m" bus (p.u.) */

  // Limits
  REAL Q_max_s;        /**< @brief Maximum series reactive power limit (p.u.) */
  REAL Q_min_s;        /**< @brief Minimum series reactive power limit (p.u.) */
  REAL Q_max_sh;       /**< @brief Maximum shunt reactive power limit (p.u.) */
  REAL Q_min_sh;       /**< @brief Minimum shunt reactive power limit (p.u.) */
  REAL i_max_s;        /**< @brief Maximum series converter current (p.u.) */
  REAL i_max_sh;       /**< @brief Maximum shunt converter current (p.u.) */
  REAL P_max_dc;       /**< @brief Maximum DC power transfer (p.u.) */
  REAL v_min_m;        /**< @brief Minimum voltage magnitude for bus "m" (p.u.) */
  REAL v_max_m;        /**< @brief Maximum voltage magnitude for bus "m" (p.u.) */
  
  // Series voltage
  REAL* v_mag_s;       /**< @brief Series voltage magnitude (p.u.) */
  REAL* v_ang_s;       /**< @brief Series voltage angel (radians wrt bus_k voltage angle) */
  REAL v_max_s;        /**< @brief Maximum series voltage magnitude (p.u.) */

  // Series admittance set-point
  REAL g;              /**< @brief Conductance set-point for constant impedance mode (p.u.) */
  REAL b;              /**< @brief Susceptance set-point for constant impedance mode (p.u.) */
  
  // Indices
  int index;           /**< @brief Facts index */
  int* index_v_mag_s;  /**< @brief Series voltage magnitude index */
  int* index_v_ang_s;  /**< @brief Series voltage angle index */

  int* index_P_k;      /**< @brief Index of active power injected into "k" bus. */
  int* index_Q_k;      /**< @brief Index of reactive power injected into "k" bus. */
  int* index_P_m;      /**< @brief Index of active power injected into "m" bus. */
  int* index_Q_m;      /**< @brief Index of reactive power injected into "m" bus. */
  int* index_Q_sh;     /**< @brief Index of reactive power injected by shunt converter. */
  int* index_Q_s;      /**< @brief Index of active power injected by series converter. */
  int* index_P_dc;     /**< @brief Index of DC power transferred from shunt to series converter. */

  // Network
  Net* net; /**< @brief Network. */

  // List
  Facts* reg_next;       /**< @brief List of facts devices regulating a bus voltage magnitude */
  Facts* next_k;         /**< @brief List of facts devices connected to a bus on the "k" side */
  Facts* next_m;         /**< @brief List of facts devices connected to a bus on the "m" side */
};

void* FACTS_array_get(void* facts_array, int index) { 
  if (facts_array) 
    return (void*)&(((Facts*)facts_array)[index]);
  else
    return NULL;
}

void FACTS_array_del(Facts* facts_array, int size) {
  int i;
  Facts* facts;
  if (facts_array) {
    for (i = 0; i < size; i++) {
      facts = &(facts_array[i]);
      free(facts->v_mag_s);
      free(facts->v_ang_s);
      free(facts->P_k);
      free(facts->Q_k);
      free(facts->P_m);
      free(facts->Q_m);
      free(facts->Q_s);
      free(facts->Q_sh);
      free(facts->P_dc);
      free(facts->Q_set);
      free(facts->P_set);
      free(facts->index_v_mag_s);
      free(facts->index_v_ang_s);
      free(facts->index_P_k);
      free(facts->index_Q_k);
      free(facts->index_P_m);
      free(facts->index_Q_m);
      free(facts->index_Q_s);
      free(facts->index_Q_sh);
      free(facts->index_P_dc);
      FACTS_set_bus_k(facts,NULL);
      FACTS_set_bus_m(facts,NULL);
      FACTS_set_reg_bus(facts,NULL);
    }
    free(facts_array);
  }  
}

Facts* FACTS_array_new(int size, int num_periods) { 
  int i;
  if (num_periods > 0) {
    Facts* facts_array = (Facts*)malloc(sizeof(Facts)*size);
    for (i = 0; i < size; i++) {
      FACTS_init(&(facts_array[i]),num_periods);
      FACTS_set_index(&(facts_array[i]),i);
      snprintf(facts_array[i].name,(size_t)(FACTS_BUFFER_SIZE-1),"%d",i);
    }
    return facts_array;
  }
  else
    return NULL;
}

void FACTS_array_show(Facts* facts_array, int size, int t) { 
  int i;
  if (facts_array) {
    for (i = 0; i < size; i++) 
      FACTS_show(&(facts_array[i]),t);
  }
}

void FACTS_clear_sensitivities(Facts* facts) {
  // Nothing
}

void FACTS_clear_flags(Facts* facts, char flag_type) {
  if (facts) {
    if (flag_type & FLAG_VARS)
      facts->vars = 0x00;
    if (flag_type & FLAG_BOUNDED)
      facts->bounded = 0x00;
    if (flag_type & FLAG_FIXED)
      facts->fixed = 0x00;
    if (flag_type & FLAG_SPARSE)
      facts->sparse = 0x00;
  }
}

void FACTS_copy_from_facts(Facts* facts, Facts* other) {

  // Local variables
  int num_periods;

  // Check
  if (!facts || !other)
    return;

  // Min num periods
  if (facts->num_periods < other->num_periods)
    num_periods = facts->num_periods;
  else
    num_periods = other->num_periods;

  // Buses
  // skip bus k
  // skip bus m
  // skip reg bus

  // Times
  // skip num periods

  // Properties
  strcpy(facts->name,other->name);

  // Flags
  facts->pre_cont_status = other->pre_cont_status;
  facts->in_service = other->in_service;
  facts->fixed = other->fixed;
  facts->bounded = other->bounded;
  facts->sparse = other->sparse;
  facts->vars = other->vars;

  // Series voltage
  memcpy(facts->v_mag_s,other->v_mag_s,num_periods*sizeof(REAL));
  memcpy(facts->v_ang_s,other->v_ang_s,num_periods*sizeof(REAL));
  facts->v_max_s = other->v_max_s;

  // Series mode
  facts->mode_s = other->mode_s;

  // Interface
  memcpy(facts->P_k,other->P_k,num_periods*sizeof(REAL));
  memcpy(facts->Q_k,other->Q_k,num_periods*sizeof(REAL));
  memcpy(facts->P_m,other->P_m,num_periods*sizeof(REAL));
  memcpy(facts->Q_m,other->Q_m,num_periods*sizeof(REAL));

  // Converters
  memcpy(facts->Q_sh,other->Q_sh,num_periods*sizeof(REAL));
  memcpy(facts->Q_s,other->Q_s,num_periods*sizeof(REAL));
  memcpy(facts->P_dc,other->P_dc,num_periods*sizeof(REAL));
  facts->Q_par = other->Q_par;
  facts->rmpct = other->rmpct;

  // Set points
  memcpy(facts->P_set,other->P_set,num_periods*sizeof(REAL));
  memcpy(facts->Q_set,other->Q_set,num_periods*sizeof(REAL));

  // Admittance
  facts->g = other->g;
  facts->b = other->b;

  // Limits
  facts->Q_max_s = other->Q_max_s;
  facts->Q_max_sh = other->Q_max_sh;
  facts->Q_min_s = other->Q_min_s;
  facts->Q_min_sh = other->Q_min_sh;
  facts->i_max_s = other->i_max_s;
  facts->i_max_sh = other->i_max_sh;
  facts->P_max_dc = other->P_max_dc;
  facts->v_min_m = other->v_min_m;
  facts->v_max_m = other->v_max_m;

  // Indices
  // skip index
  memcpy(facts->index_v_mag_s,other->index_v_mag_s,num_periods*sizeof(int));
  memcpy(facts->index_v_ang_s,other->index_v_ang_s,num_periods*sizeof(int));
  memcpy(facts->index_P_k,other->index_P_k,num_periods*sizeof(int));
  memcpy(facts->index_Q_k,other->index_Q_k,num_periods*sizeof(int));
  memcpy(facts->index_P_m,other->index_P_m,num_periods*sizeof(int));
  memcpy(facts->index_Q_m,other->index_Q_m,num_periods*sizeof(int));
  memcpy(facts->index_Q_s,other->index_Q_s,num_periods*sizeof(int));
  memcpy(facts->index_Q_sh,other->index_Q_sh,num_periods*sizeof(int));
  memcpy(facts->index_P_dc,other->index_P_dc,num_periods*sizeof(int));
  
  // List
  // skip next
}

short int FACTS_get_pre_cont_status(void* facts) {
  if (facts)
    return ((Facts*)facts)->pre_cont_status;
  else
    return 0;
}

char FACTS_get_flags_vars(Facts* facts) {
  if (facts)
    return facts->vars;
  else
    return 0;
}

char FACTS_get_flags_fixed(Facts* facts) {
  if (facts)
    return facts->fixed;
  else
    return 0;
}

char FACTS_get_flags_bounded(Facts* facts) {
  if (facts)
    return facts->bounded;
  else
    return 0;
}

char FACTS_get_flags_sparse(Facts* facts) {
  if (facts)
    return facts->sparse;
  else
    return 0;
}

char* FACTS_get_name(Facts* facts) {
  if (facts)
    return facts->name;
  else
    return NULL;
}

int FACTS_get_num_periods(Facts* facts) {
  if (facts)
    return facts->num_periods;
  else
    return 0;
}

char FACTS_get_obj_type(void* facts) {
  if (facts)
    return OBJ_FACTS;
  else
    return OBJ_UNKNOWN;
}

Bus* FACTS_get_bus_k(Facts* facts) {
  if (facts)
    return facts->bus_k;
  else
    return NULL;
}

Bus* FACTS_get_bus_m(Facts* facts) {
  if (facts)
    return facts->bus_m;
  else
    return NULL;
}

Bus* FACTS_get_reg_bus(Facts* facts) {
  if (facts)
    return facts->reg_bus;
  else
    return NULL;
}

int FACTS_get_index(Facts* facts) {
  if (facts)
    return facts->index;
  else
    return -1;
}

int FACTS_get_index_v_mag_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_v_mag_s[t];
  else
    return -1;
}

int FACTS_get_index_v_ang_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_v_ang_s[t];
  else
    return -1;
}

int FACTS_get_index_P_k(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_P_k[t];
  else
    return -1;
}

int FACTS_get_index_Q_k(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_Q_k[t];
  else
    return -1;
}

int FACTS_get_index_P_m(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_P_m[t];
  else
    return -1;
}

int FACTS_get_index_Q_m(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_Q_m[t];
  else
    return -1;
}

int FACTS_get_index_Q_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_Q_s[t];
  else
    return -1;
}

int FACTS_get_index_Q_sh(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_Q_sh[t];
  else
    return -1;
}

int FACTS_get_index_P_dc(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->index_P_dc[t];
  else
    return -1;
}


int* FACTS_get_index_v_mag_s_array(Facts* facts) {
  if (facts)
    return facts->index_v_mag_s;
  else
    return NULL;
}

int* FACTS_get_index_v_ang_s_array(Facts* facts) {
  if (facts)
    return facts->index_v_ang_s;
  else
    return NULL;
}

int* FACTS_get_index_P_k_array(Facts* facts) {
  if (facts)
    return facts->index_P_k;
  else
    return NULL;
}

int* FACTS_get_index_P_m_array(Facts* facts) {
  if (facts)
    return facts->index_P_m;
  else
    return NULL;
}

int* FACTS_get_index_Q_k_array(Facts* facts) {
  if (facts)
    return facts->index_Q_k;
  else
    return NULL;
}

int* FACTS_get_index_Q_m_array(Facts* facts) {
  if (facts)
    return facts->index_Q_m;
  else
    return NULL;
}

int* FACTS_get_index_P_dc_array(Facts* facts) {
  if (facts)
    return facts->index_P_dc;
  else
    return NULL;
}

int* FACTS_get_index_Q_s_array(Facts* facts) {
  if (facts)
    return facts->index_Q_s;
  else
    return NULL;
}

int* FACTS_get_index_Q_sh_array(Facts* facts) {
  if (facts)
    return facts->index_Q_sh;
  else
    return NULL;
}

Facts* FACTS_get_next_k(Facts* facts) {
  if (facts)
    return facts->next_k;
  else
    return NULL;
}

Facts* FACTS_get_next_m(Facts* facts) {
  if (facts)
    return facts->next_m;
  else
    return NULL;
}

Facts* FACTS_get_reg_next(Facts* facts) {
  if (facts)
    return facts->reg_next;
  else
    return NULL;
}

REAL FACTS_get_v_max_s(Facts* facts) {
  if (facts)
    return facts->v_max_s;
  else
    return 0;
}

REAL FACTS_get_g(Facts* facts) {
  if (facts)
    return facts->g;
  else
    return 0;
}

REAL FACTS_get_b(Facts* facts) {
  if (facts)
    return facts->b;
  else
    return 0;
}

char FACTS_get_mode_s(Facts* facts) {
  if (facts)
    return facts->mode_s;
  else
    return FACTS_SERIES_MODE_DISABLED;
}

REAL FACTS_get_v_mag_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->v_mag_s[t];
  else
    return 0;
}

REAL FACTS_get_v_ang_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->v_ang_s[t];
  else
    return 0;
}

REAL FACTS_get_P_k(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->P_k[t];
  else
    return 0;
}

REAL FACTS_get_P_m(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->P_m[t];
  else
    return 0;
}

REAL FACTS_get_Q_k(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->Q_k[t];
  else
    return 0;
}

REAL FACTS_get_Q_m(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->Q_m[t];
  else
    return 0;
}

REAL FACTS_get_Q_sh(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->Q_sh[t];
  else
    return 0;
}

REAL FACTS_get_Q_s(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->Q_s[t];
  else
    return 0;
}

REAL FACTS_get_P_dc(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->P_dc[t];
  else
    return 0;
}

REAL FACTS_get_P_set(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->P_set[t];
  else
    return 0;
}

REAL FACTS_get_Q_set(Facts* facts, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    return facts->Q_set[t];
  else
    return 0;
}

REAL* FACTS_get_v_mag_s_array(Facts* facts) {
  if (facts)
    return facts->v_mag_s;
  else
    return NULL;
}

REAL* FACTS_get_v_ang_s_array(Facts* facts) {
  if (facts)
    return facts->v_ang_s;
  else
    return NULL;
}

REAL* FACTS_get_P_k_array(Facts* facts) {
  if (facts)
    return facts->P_k;
  else
    return NULL;
}

REAL* FACTS_get_P_m_array(Facts* facts) {
  if (facts)
    return facts->P_m;
  else
    return NULL;
}

REAL* FACTS_get_Q_k_array(Facts* facts) {
  if (facts)
    return facts->Q_k;
  else
    return NULL;
}

REAL* FACTS_get_Q_m_array(Facts* facts) {
  if (facts)
    return facts->Q_m;
  else
    return NULL;
}

REAL* FACTS_get_Q_sh_array(Facts* facts) {
  if (facts)
    return facts->Q_sh;
  else
    return NULL;
}

REAL* FACTS_get_Q_s_array(Facts* facts) {
  if (facts)
    return facts->Q_s;
  else
    return NULL;
}

REAL* FACTS_get_P_dc_array(Facts* facts) {
  if (facts)
    return facts->P_dc;
  else
    return NULL;
}

REAL* FACTS_get_P_set_array(Facts* facts) {
  if (facts)
    return facts->P_set;
  else
    return NULL;
}

REAL* FACTS_get_Q_set_array(Facts* facts) {
  if (facts)
    return facts->Q_set;
  else
    return NULL;
}

REAL FACTS_get_Q_par(Facts* facts) {
  if (facts)
    return facts->Q_par;
  else
    return 0;
}

REAL FACTS_get_rmpct(Facts* facts) {
  if (facts)
    return facts->rmpct;
  else
    return 0;
}

REAL FACTS_get_Q_max_s(Facts* facts) {
  if (facts)
    return facts->Q_max_s;
  else
    return 0;
}

REAL FACTS_get_Q_max_sh(Facts* facts) {
  if (facts)
    return facts->Q_max_sh;
  else
    return 0;
}

REAL FACTS_get_Q_min_s(Facts* facts) {
  if (facts)
    return facts->Q_min_s;
  else
    return 0;
}

REAL FACTS_get_Q_min_sh(Facts* facts) {
  if (facts)
    return facts->Q_min_sh;
  else
    return 0;
}

REAL FACTS_get_i_max_s(Facts* facts) {
  if (facts)
    return facts->i_max_s;
  else
    return 0;
}

REAL FACTS_get_i_max_sh(Facts* facts) {
  if (facts)
    return facts->i_max_sh;
  else
    return 0;
}

REAL FACTS_get_P_max_dc(Facts* facts) {
  if (facts)
    return facts->P_max_dc;
  else
    return 0;
}

REAL FACTS_get_v_min_m(Facts* facts) {
  if (facts)
    return facts->v_min_m;
  else
    return 0;
}

REAL FACTS_get_v_max_m(Facts* facts) {
  if (facts)
    return facts->v_max_m;
  else
    return 0;
}

void FACTS_get_var_values(Facts* facts, Vec* values, int code) {
 
  // Local vars
  int t;
 
  // No laod
  if (!facts)
    return;

  // Time loop
  for (t = 0; t < facts->num_periods; t++) {

    if (facts->vars & FACTS_VAR_VMAG_S) { // series voltage magnitude
      switch(code) {
 
      case UPPER_LIMITS:
        if (facts->bounded & FACTS_VAR_VMAG_S)
          VEC_set(values,facts->index_v_mag_s[t],facts->v_max_s);
        else
          VEC_set(values,facts->index_v_mag_s[t],FACTS_INF_VMAG_S);
        break;

      case LOWER_LIMITS:
        if (facts->bounded & FACTS_VAR_VMAG_S)
          VEC_set(values,facts->index_v_mag_s[t],0.);
        else
          VEC_set(values,facts->index_v_mag_s[t],-FACTS_INF_VMAG_S);
        break;

      default:
        VEC_set(values,facts->index_v_mag_s[t],facts->v_mag_s[t]);
      }
    }
    
    if (facts->vars & FACTS_VAR_VANG_S) { // series voltage angle
      switch(code) {

      case UPPER_LIMITS:
        VEC_set(values,facts->index_v_ang_s[t],FACTS_INF_VANG_S);
        break;
        
      case LOWER_LIMITS:
        VEC_set(values,facts->index_v_ang_s[t],-FACTS_INF_VANG_S);
        break;

      default:
        VEC_set(values,facts->index_v_ang_s[t],facts->v_ang_s[t]);
      }
    }

    if (facts->vars & FACTS_VAR_P) { // active power
      switch(code) {

      case UPPER_LIMITS:
        VEC_set(values,facts->index_P_k[t],FACTS_INF_P);
        VEC_set(values,facts->index_P_m[t],FACTS_INF_P);
        if (facts->bounded & FACTS_VAR_P)
          VEC_set(values,facts->index_P_dc[t],facts->P_max_dc);
        else
          VEC_set(values,facts->index_P_dc[t],FACTS_INF_P);
        break;

      case LOWER_LIMITS:
        VEC_set(values,facts->index_P_k[t],-FACTS_INF_P);
        VEC_set(values,facts->index_P_m[t],-FACTS_INF_P);
        if (facts->bounded & FACTS_VAR_P)
          VEC_set(values,facts->index_P_dc[t],-facts->P_max_dc);
        else
          VEC_set(values,facts->index_P_dc[t],-FACTS_INF_P);
        break;

      default:
        VEC_set(values,facts->index_P_k[t],facts->P_k[t]);
        VEC_set(values,facts->index_P_m[t],facts->P_m[t]);
        VEC_set(values,facts->index_P_dc[t],facts->P_dc[t]);
      }
    }

    if (facts->vars & FACTS_VAR_Q) { // reactive power
      switch(code) {

      case UPPER_LIMITS:
        VEC_set(values,facts->index_Q_k[t],FACTS_INF_Q);
        VEC_set(values,facts->index_Q_m[t],FACTS_INF_Q);
        if (facts->bounded & FACTS_VAR_Q) {
          VEC_set(values,facts->index_Q_s[t],facts->Q_max_s);
          VEC_set(values,facts->index_Q_sh[t],facts->Q_max_sh);
        }
        else {
          VEC_set(values,facts->index_Q_s[t],FACTS_INF_Q);
          VEC_set(values,facts->index_Q_sh[t],FACTS_INF_Q);
        }
        break;

      case LOWER_LIMITS:
        VEC_set(values,facts->index_Q_k[t],-FACTS_INF_Q);
        VEC_set(values,facts->index_Q_m[t],-FACTS_INF_Q);
        if (facts->bounded & FACTS_VAR_Q) {
          VEC_set(values,facts->index_Q_s[t],facts->Q_min_s);
          VEC_set(values,facts->index_Q_sh[t],facts->Q_min_sh);
        }
        else {
          VEC_set(values,facts->index_Q_s[t],-FACTS_INF_Q);
          VEC_set(values,facts->index_Q_sh[t],-FACTS_INF_Q);
        }
        break;

      default:
        VEC_set(values,facts->index_Q_k[t],facts->Q_k[t]);
        VEC_set(values,facts->index_Q_m[t],facts->Q_m[t]);
        VEC_set(values,facts->index_Q_s[t],facts->Q_s[t]);
        VEC_set(values,facts->index_Q_sh[t],facts->Q_sh[t]);
      }
    }
  }
}

char* FACTS_get_var_info_string(Facts* facts, int index) {

  // Local variables
  char* info;

  //Check
  if (!facts)
    return NULL;

  // Series voltage magnitude
  if ((facts->vars & FACTS_VAR_VMAG_S) &&
      index >= facts->index_v_mag_s[0] &&
      index <= facts->index_v_mag_s[facts->num_periods-1]) {
    info = (char*)malloc(FACTS_BUFFER_SIZE*sizeof(char));
    snprintf(info,FACTS_BUFFER_SIZE*sizeof(char),
             "facts:%d:series voltage magnitude:%d",facts->index,index-facts->index_v_mag_s[0]);
    return info;
  }

  // Series voltage angle
  if ((facts->vars & FACTS_VAR_VANG_S) &&
      index >= facts->index_v_ang_s[0] &&
      index <= facts->index_v_ang_s[facts->num_periods-1]) {
    info = (char*)malloc(FACTS_BUFFER_SIZE*sizeof(char));
    snprintf(info,FACTS_BUFFER_SIZE*sizeof(char),
             "facts:%d:series voltage angle:%d",facts->index,index-facts->index_v_ang_s[0]);
    return info;
  }

  // Active power

  // Reactive power

  // Return
  return NULL;
}

int FACTS_get_num_vars(void* vfacts, unsigned char var, int t_start, int t_end) {

  // Local vars
  Facts* facts = (Facts*)vfacts;
  int num_vars = 0;
  int dt;

  // Checks
  if (!facts)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > facts->num_periods-1)
    t_end = facts->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & FACTS_VAR_VMAG_S) && (facts->vars & FACTS_VAR_VMAG_S)) // series voltage magnitude
    num_vars += dt;
  if ((var & FACTS_VAR_VANG_S) && (facts->vars & FACTS_VAR_VANG_S)) // series voltage angle
    num_vars += dt;
  if ((var & FACTS_VAR_P) && (facts->vars & FACTS_VAR_P)) // active power (P_k, P_m, P_dc)
    num_vars += 3*dt;
  if ((var & FACTS_VAR_Q) && (facts->vars & FACTS_VAR_Q)) // active power (Q_k, Q_m, Q_s, Q_sh)
    num_vars += 4*dt;
  return num_vars;
}

Vec* FACTS_get_var_indices(void* vfacts, unsigned char var, int t_start, int t_end) {

  // Local vars
  Facts* facts = (Facts*)vfacts;
  Vec* indices;
  int offset = 0;
  int t;

  // Checks
  if (!facts)
    return NULL;
  if (t_start < 0)
    t_start = 0;
  if (t_end > facts->num_periods-1)
    t_end = facts->num_periods-1;

  // Indices
  indices = VEC_new(FACTS_get_num_vars(vfacts,var,t_start,t_end));
  if ((var & FACTS_VAR_VMAG_S) && (facts->vars & FACTS_VAR_VMAG_S)) { // series voltage magnitude
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,facts->index_v_mag_s[t]);
      offset++;
    }
  }
  if ((var & FACTS_VAR_VANG_S) && (facts->vars & FACTS_VAR_VANG_S)) { // series voltage angle
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,facts->index_v_ang_s[t]);
      offset++;
    }
  }
  if ((var & FACTS_VAR_P) && (facts->vars & FACTS_VAR_P)) { // active power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,facts->index_P_k[t]);
      VEC_set(indices,offset+1,facts->index_P_m[t]);
      VEC_set(indices,offset+2,facts->index_P_dc[t]);
      offset += 3;
    }
  }
  if ((var & FACTS_VAR_Q) && (facts->vars & FACTS_VAR_Q)) { // reactive power
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,facts->index_Q_k[t]);
      VEC_set(indices,offset+1,facts->index_Q_m[t]);
      VEC_set(indices,offset+2,facts->index_Q_s[t]);
      VEC_set(indices,offset+3,facts->index_Q_sh[t]);
      offset += 4;
    }
  }
  return indices;
}

char* FACTS_get_json_string(Facts* facts, char* output) {

  // Local variables
  char temp[FACTS_JSON_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No facts
  if (!facts)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*FACTS_BUFFER_SIZE*FACTS_NUM_JSON_FIELDS*facts->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"index",facts->index,FALSE);
  JSON_int(temp,output,"num_periods",facts->num_periods,FALSE);
  JSON_str(temp,output,"name",facts->name,FALSE);
  JSON_bool(temp,output,"in_service",facts->in_service,FALSE);
  JSON_obj(temp,output,"bus_k",facts->bus_k,BUS_get_index,FALSE);
  JSON_obj(temp,output,"bus_m",facts->bus_m,BUS_get_index,FALSE);
  JSON_obj(temp,output,"reg_bus",facts->reg_bus,BUS_get_index,FALSE);
  JSON_int(temp,output,"mode_s",facts->mode_s,FALSE);

  JSON_array_float(temp,output,"P_k",facts->P_k,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_k",facts->Q_k,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"P_m",facts->P_m,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_m",facts->Q_m,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_sh",facts->Q_sh,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_s",facts->Q_s,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"P_dc",facts->P_dc,facts->num_periods,FALSE);
  JSON_float(temp,output,"Q_par",facts->Q_par,FALSE);
  JSON_float(temp,output,"rmpct",facts->rmpct,FALSE);
  JSON_array_float(temp,output,"P_set",facts->P_set,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"Q_set",facts->Q_set,facts->num_periods,FALSE);

  JSON_float(temp,output,"Q_max_s",facts->Q_max_s,FALSE);
  JSON_float(temp,output,"Q_min_s",facts->Q_min_s,FALSE);
  JSON_float(temp,output,"Q_max_sh",facts->Q_max_sh,FALSE);
  JSON_float(temp,output,"Q_min_sh",facts->Q_min_sh,FALSE);
  JSON_float(temp,output,"i_max_s",facts->i_max_s,FALSE);
  JSON_float(temp,output,"i_max_sh",facts->i_max_sh,FALSE);
  JSON_float(temp,output,"P_max_dc",facts->P_max_dc,FALSE);
  JSON_float(temp,output,"v_max_m",facts->v_max_m,FALSE);
  JSON_float(temp,output,"v_min_m",facts->v_min_m,FALSE);

  JSON_array_float(temp,output,"v_mag_s",facts->v_mag_s,facts->num_periods,FALSE);
  JSON_array_float(temp,output,"v_ang_s",facts->v_ang_s,facts->num_periods,FALSE);
  JSON_float(temp,output,"v_max_s",facts->v_max_s,FALSE);

  JSON_float(temp,output,"g",facts->g,FALSE);
  JSON_float(temp,output,"b",facts->b,TRUE);
  
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL FACTS_has_flags(void* vfacts, char flag_type, unsigned char mask) {
  Facts* facts = (Facts*)vfacts;
  if (facts) {
    if (flag_type == FLAG_VARS)
      return (facts->vars & mask) == mask;
    else if (flag_type == FLAG_BOUNDED)
      return (facts->bounded & mask) == mask;
    else if (flag_type == FLAG_FIXED)
      return (facts->fixed & mask) == mask;
    else if (flag_type == FLAG_SPARSE)
      return (facts->sparse & mask) == mask;
    return FALSE;
  }
  else
    return FALSE;
}

BOOL FACTS_has_properties(void* vfacts, char prop) {
  Facts* facts = (Facts*)vfacts;
  if (!facts)
    return FALSE;
  return TRUE;
}

void FACTS_init(Facts* facts, int num_periods) {

  // Local vars
  int T;

  // No facts
  if (!facts)
    return;

  T = num_periods;
  facts->num_periods = num_periods;
  ARRAY_clear(facts->name,char,FACTS_BUFFER_SIZE);
  
  facts->bus_k = NULL;
  facts->bus_m = NULL;
  facts->reg_bus = NULL;

  facts->pre_cont_status = PRE_CONT_UNSET;
  facts->in_service = TRUE;
  facts->fixed = 0x00;
  facts->bounded = 0x00;
  facts->sparse = 0x00;
  facts->vars = 0x00;
      
  facts->index = -1;

  ARRAY_zalloc(facts->v_mag_s,REAL,T);
  ARRAY_zalloc(facts->v_ang_s,REAL,T);

  facts->mode_s = FACTS_SERIES_MODE_DISABLED;

  ARRAY_zalloc(facts->P_k,REAL,T);
  ARRAY_zalloc(facts->Q_k,REAL,T);
  ARRAY_zalloc(facts->P_m,REAL,T);
  ARRAY_zalloc(facts->Q_m,REAL,T);

  ARRAY_zalloc(facts->Q_sh,REAL,T);
  ARRAY_zalloc(facts->Q_s,REAL,T);
  ARRAY_zalloc(facts->P_dc,REAL,T);
  facts->Q_par = 1.;
  facts->rmpct = 100.;

  ARRAY_zalloc(facts->P_set,REAL,T);
  ARRAY_zalloc(facts->Q_set,REAL,T);

  facts->Q_max_s = 0.;
  facts->Q_max_sh = 0.;
  facts->Q_min_s = 0.;
  facts->Q_min_sh = 0.;
  facts->i_max_s = 0.;
  facts->i_max_sh = 0.;
  facts->P_max_dc = 0.;
  facts->v_min_m = 1.;
  facts->v_max_m = 1.;  
  
  ARRAY_zalloc(facts->index_v_mag_s,int,T);
  ARRAY_zalloc(facts->index_v_ang_s,int,T);
  ARRAY_zalloc(facts->index_P_k,int,T);
  ARRAY_zalloc(facts->index_P_m,int,T);
  ARRAY_zalloc(facts->index_P_dc,int,T);
  ARRAY_zalloc(facts->index_Q_k,int,T);
  ARRAY_zalloc(facts->index_Q_m,int,T);
  ARRAY_zalloc(facts->index_Q_s,int,T);
  ARRAY_zalloc(facts->index_Q_sh,int,T);

  facts->v_max_s = 0.;

  facts->g = 0;
  facts->b = 0;

  facts->net = NULL;

  facts->reg_next = NULL;
  facts->next_k = NULL;
  facts->next_m = NULL;
}

BOOL FACTS_is_regulator(Facts* facts) {
  if (facts)
    return facts->reg_bus != NULL;
  else
    return FALSE;
}

BOOL FACTS_is_STATCOM(Facts* facts) {
  return FACTS_is_series_link_disabled(facts);
}

BOOL FACTS_is_SSSC(Facts* facts) {
  return ((!FACTS_is_series_link_disabled(facts)) &&
          (FACTS_get_P_max_dc(facts) == 0.) &&
          (FACTS_get_i_max_sh(facts) == 0.));
}

BOOL FACTS_is_UPFC(Facts* facts) {
  return  ((!FACTS_is_series_link_disabled(facts)) &&
           (FACTS_get_P_max_dc(facts) > 0.) &&
           (FACTS_get_i_max_sh(facts) > 0.));
}

BOOL FACTS_is_series_link_disabled(Facts* facts) {
  if (facts)
    return facts->mode_s == FACTS_SERIES_MODE_DISABLED;
  else
    return FALSE;
}

BOOL FACTS_is_series_link_bypassed(Facts* facts) {
  if (facts)
    return facts->mode_s == FACTS_SERIES_MODE_BYPASS;
  else
    return FALSE;
}

BOOL FACTS_is_in_normal_series_mode(Facts* facts) {
  if (facts)
    return facts->mode_s == FACTS_SERIES_MODE_NORMAL;
  else
    return FALSE;
}

BOOL FACTS_is_in_constant_series_z_mode(Facts* facts) {
  if (facts)
    return facts->mode_s == FACTS_SERIES_MODE_CZ;
  else
    return FALSE;
}

BOOL FACTS_is_in_constant_series_v_mode(Facts* facts) {
  if (facts)
    return facts->mode_s == FACTS_SERIES_MODE_CV;
  else
    return FALSE;
}

BOOL FACTS_is_in_service(void* facts) {
  if (facts)
    return (((Facts*)facts)->in_service &&
            BUS_is_in_service(((Facts*)facts)->bus_k) &&
            (BUS_is_in_service(((Facts*)facts)->bus_m) ||
             FACTS_is_series_link_disabled((Facts*)facts)));
  else
    return FALSE;
}

BOOL FACTS_is_equal(Facts* facts, Facts* other) {
  return facts == other;
}

Facts* FACTS_list_k_add(Facts* facts_list, Facts* facts) {
  LIST_add(Facts,facts_list,facts,next_k);
  return facts_list;
}

Facts* FACTS_list_k_del(Facts* facts_list, Facts* facts) {
  LIST_del(Facts,facts_list,facts,next_k);
  return facts_list;
}

int FACTS_list_k_len(Facts* facts_list) {
  int len;
  LIST_len(Facts,facts_list,next_k,len);
  return len;
}

Facts* FACTS_list_m_add(Facts* facts_list, Facts* facts) {
  LIST_add(Facts,facts_list,facts,next_m);
  return facts_list;
}

Facts* FACTS_list_m_del(Facts* facts_list, Facts* facts) {
  LIST_del(Facts,facts_list,facts,next_m);
  return facts_list;
}

int FACTS_list_m_len(Facts* facts_list) {
  int len;
  LIST_len(Facts,facts_list,next_m,len);
  return len;
}

Facts* FACTS_list_reg_add(Facts* facts_list, Facts* facts) {
  LIST_add(Facts,facts_list,facts,reg_next);
  return facts_list;
}

Facts* FACTS_list_reg_del(Facts* facts_list, Facts* facts) {
  LIST_del(Facts,facts_list,facts,reg_next);
  return facts_list;
}

int FACTS_list_reg_len(Facts* facts_list) {
  int len;
  LIST_len(Facts,facts_list,reg_next,len);
  return len;
}

Facts* FACTS_new(int num_periods) {
  if (num_periods > 0) {
    Facts* facts = (Facts*)malloc(sizeof(Facts));
    FACTS_init(facts,num_periods);
    return facts;
  }
  else
    return NULL;
}

void FACTS_set_pre_cont_status(Facts* facts, short int pre_cont_status) {
  if (facts && BUS_is_in_service(facts->bus_k) &&
      (BUS_is_in_service(facts->bus_m) || FACTS_is_series_link_disabled(facts)))
    facts->pre_cont_status = pre_cont_status;
}

void FACTS_set_in_service(Facts* facts, BOOL in_service) {
  if (facts && BUS_is_in_service(facts->bus_k) &&
      (BUS_is_in_service(facts->bus_m) || FACTS_is_series_link_disabled(facts))) {
    if (facts->in_service != in_service)
      NET_inc_state_tag(facts->net);
    facts->in_service = in_service;
  }
}

void FACTS_set_network(Facts* facts, void* net) {
  if (facts)
    facts->net = (Net*)net;
}

void FACTS_set_name(Facts* facts, char* name) {
  if (facts)
    strncpy(facts->name,name,(size_t)(FACTS_BUFFER_SIZE-1));
}

void FACTS_set_bus_k(Facts* facts, Bus* bus) {
  Bus* old_bus;
  if (facts) {
    old_bus = facts->bus_k;
    facts->bus_k = NULL;
    BUS_del_facts_k(old_bus,facts);
    facts->bus_k = bus;
    BUS_add_facts_k(facts->bus_k,facts);
  }
}

void FACTS_set_bus_m(Facts* facts, Bus* bus) {
  Bus* old_bus;
  if (facts) {
    old_bus = facts->bus_m;
    facts->bus_m = NULL;
    BUS_del_facts_m(old_bus,facts);
    facts->bus_m = bus;
    BUS_add_facts_m(facts->bus_m,facts);
  }
}

void FACTS_set_reg_bus(Facts* facts, Bus* bus) {
  Bus* old_bus;
  if (facts) {
    old_bus = facts->reg_bus;
    facts->reg_bus = NULL;
    BUS_del_reg_facts(old_bus,facts);
    facts->reg_bus = bus;
    BUS_add_reg_facts(facts->reg_bus,facts);
  }
}

void FACTS_set_index(Facts* facts, int index) { 
  if (facts)
    facts->index = index;
}

void FACTS_set_v_mag_s(Facts* facts, REAL v_mag, int t) { 
  if (facts && t >= 0 && t < facts->num_periods)
    facts->v_mag_s[t] = v_mag;
}

void FACTS_set_v_ang_s(Facts* facts, REAL v_ang, int t) { 
  if (facts && t >= 0 && t < facts->num_periods)
    facts->v_ang_s[t] = v_ang;
}

void FACTS_set_v_max_s(Facts* facts, REAL v_max) {
  if (facts)
    facts->v_max_s = v_max;
}

void FACTS_set_mode_s(Facts* facts, char mode) {
  if (facts)
    facts->mode_s = mode;
}

void FACTS_set_P_k(Facts* facts, REAL P, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->P_k[t] = P;
}

void FACTS_set_P_m(Facts* facts, REAL P, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->P_m[t] = P;
}

void FACTS_set_Q_k(Facts* facts, REAL Q, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->Q_k[t] = Q;
}

void FACTS_set_Q_m(Facts* facts, REAL Q, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->Q_m[t] = Q;
}

void FACTS_set_Q_s(Facts* facts, REAL Q, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->Q_s[t] = Q;
}

void FACTS_set_Q_sh(Facts* facts, REAL Q, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->Q_sh[t] = Q;
}

void FACTS_set_P_dc(Facts* facts, REAL P, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->P_dc[t] = P;
}

void FACTS_set_Q_par(Facts* facts, REAL Q_par) {
  if (facts)
    facts->Q_par = Q_par;
}

void FACTS_set_rmpct(Facts* facts, REAL rmpct) {

  // Local variables
  Facts* f;
  Bus* bus;

  if (facts) {
    bus = FACTS_get_bus_k(facts);
    // Change rmpct for all parallel FACTS
    for (f=BUS_get_facts_k(bus); f !=NULL; f = FACTS_get_next_k(f)) {
      f->rmpct = rmpct;
    }
    
    bus = FACTS_get_bus_m(facts);
    for (f=BUS_get_facts_m(bus); f !=NULL; f = FACTS_get_next_m(f)) {
      f->rmpct = rmpct;
    }
  }
}

void FACTS_set_P_set(Facts* facts, REAL P, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->P_set[t] = P;
}

void FACTS_set_Q_set(Facts* facts, REAL Q, int t) {
  if (facts && t >= 0 && t < facts->num_periods)
    facts->Q_set[t] = Q;
}

void FACTS_set_Q_max_s(Facts* facts, REAL Q_max) {
  if (facts)
    facts->Q_max_s = Q_max;
}

void FACTS_set_Q_max_sh(Facts* facts, REAL Q_max) {
  if (facts)
    facts->Q_max_sh = Q_max;
}

void FACTS_set_Q_min_s(Facts* facts, REAL Q_max) {
  if (facts)
    facts->Q_min_s = Q_max;
}

void FACTS_set_Q_min_sh(Facts* facts, REAL Q_max) {
  if (facts)
    facts->Q_min_sh = Q_max;
}

void FACTS_set_i_max_s(Facts* facts, REAL i_max) {
  if (facts)
    facts->i_max_s = i_max;
}

void FACTS_set_i_max_sh(Facts* facts, REAL i_max) {
  if (facts)
    facts->i_max_sh = i_max;
}

void FACTS_set_P_max_dc(Facts* facts, REAL P_max) {
  if (facts)
    facts->P_max_dc = P_max;
}

void FACTS_set_v_min_m(Facts* facts, REAL v_min) {
  if (facts)
    facts->v_min_m = v_min;
}

void FACTS_set_v_max_m(Facts* facts, REAL v_max) {
  if (facts)
    facts->v_max_m = v_max;
}

void FACTS_set_g(Facts* facts, REAL g) {
  if (facts)
    facts->g = g;
}

void FACTS_set_b(Facts* facts, REAL b) {
  if (facts)
    facts->b = b;
}

int FACTS_set_flags(void* vfacts, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Facts* facts = (Facts*)vfacts;
  int t;

  // Check facts
  if (!facts)
    return 0;

  // Set flag pointer
  if (flag_type == FLAG_VARS)
    flags_ptr = &(facts->vars);
  else if (flag_type == FLAG_FIXED)
    flags_ptr = &(facts->fixed);
  else if (flag_type == FLAG_BOUNDED)
    flags_ptr = &(facts->bounded);
  else if (flag_type == FLAG_SPARSE)
    flags_ptr = &(facts->sparse);
  else
    return index;

  // Set flags
  if (!((*flags_ptr) & FACTS_VAR_VMAG_S) && (mask & FACTS_VAR_VMAG_S)) { // series voltage magnitude
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < facts->num_periods; t++)
        facts->index_v_mag_s[t] = index+t;
    }
    (*flags_ptr) |= FACTS_VAR_VMAG_S;
    index += facts->num_periods;
  }
  if (!((*flags_ptr) & FACTS_VAR_VANG_S) && (mask & FACTS_VAR_VANG_S)) { // series voltage angle
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < facts->num_periods; t++)
        facts->index_v_ang_s[t] = index+t;
    }
    (*flags_ptr) |= FACTS_VAR_VANG_S;
    index += facts->num_periods;
  }
  if (!((*flags_ptr) & FACTS_VAR_P) && (mask & FACTS_VAR_P)) { // active power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < facts->num_periods; t++)
        facts->index_P_k[t] = index+t;
      for (t = 0; t < facts->num_periods; t++)
        facts->index_P_m[t] = index+facts->num_periods+t;
      for (t = 0; t < facts->num_periods; t++)
        facts->index_P_dc[t] = index+2*facts->num_periods+t;
    }
    (*flags_ptr) |= FACTS_VAR_P;
    index += 3*facts->num_periods;
  }
  if (!((*flags_ptr) & FACTS_VAR_Q) && (mask & FACTS_VAR_Q)) { // reactive power
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < facts->num_periods; t++)
        facts->index_Q_k[t] = index+t;
      for (t = 0; t < facts->num_periods; t++)
        facts->index_Q_m[t] = index+facts->num_periods+t;
      for (t = 0; t < facts->num_periods; t++)
        facts->index_Q_s[t] = index+2*facts->num_periods+t;
      for (t = 0; t < facts->num_periods; t++)
        facts->index_Q_sh[t] = index+3*facts->num_periods+t;
    }
    (*flags_ptr) |= FACTS_VAR_Q;
    index += 4*facts->num_periods;
  }
  return index;
}

void FACTS_set_var_values(Facts* facts, Vec* values) {

  // Local vars
  int t;

  // No facts
  if (!facts)
    return;

  // Time loop
  for (t = 0; t < facts->num_periods; t++) {

    if (facts->vars & FACTS_VAR_VMAG_S) // series voltage magnitude (p.u.)
      facts->v_mag_s[t] = VEC_get(values,facts->index_v_mag_s[t]);

    if (facts->vars & FACTS_VAR_VANG_S) // series voltage angle (p.u.)
      facts->v_ang_s[t] = VEC_get(values,facts->index_v_ang_s[t]);

    if (facts->vars & FACTS_VAR_P) { // active power (p.u.)
      facts->P_k[t] = VEC_get(values,facts->index_P_k[t]);
      facts->P_m[t] = VEC_get(values,facts->index_P_m[t]);
      facts->P_dc[t] = VEC_get(values,facts->index_P_dc[t]);      
    }

    if (facts->vars & FACTS_VAR_P) { // reactive power (p.u.)
      facts->Q_k[t] = VEC_get(values,facts->index_Q_k[t]);
      facts->Q_m[t] = VEC_get(values,facts->index_Q_m[t]);
      facts->Q_s[t] = VEC_get(values,facts->index_Q_s[t]);
      facts->Q_sh[t] = VEC_get(values,facts->index_Q_sh[t]);
    }
  }
}

void FACTS_show(Facts* facts, int t) { 
  if (facts)
    printf("facts %d\t%d\t%d\n",
           BUS_get_number(facts->bus_k),
           BUS_get_number(facts->bus_m),
           facts->index);
}

void FACTS_propagate_data_in_time(Facts* facts, int start, int end) {
  int t;
  if (facts) {
    if (start < 0)
      start = 0;
    if (end > facts->num_periods)
      end = facts->num_periods;
    for (t = start+1; t < end; t++) {
      facts->v_mag_s[t] = facts->v_mag_s[start];
      facts->v_ang_s[t] = facts->v_ang_s[start];
      facts->P_k[t] = facts->P_k[start];
      facts->P_m[t] = facts->P_m[start];
      facts->Q_k[t] = facts->Q_k[start];
      facts->Q_m[t] = facts->Q_m[start];
      facts->Q_set[t] = facts->Q_set[start];
      facts->P_set[t] = facts->P_set[start];
      facts->Q_sh[t] = facts->Q_sh[start];
      facts->Q_s[t] = facts->Q_s[start];
      facts->P_dc[t] = facts->P_dc[start];
    }
  }
}
