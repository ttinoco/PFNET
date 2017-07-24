/** @file bus.c
 *  @brief This file defines the Bus data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/bus.h>
#include <pfnet/gen.h>
#include <pfnet/branch.h>
#include <pfnet/load.h>
#include <pfnet/shunt.h>
#include <pfnet/vargen.h>
#include <pfnet/bat.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>

struct Bus {

  // Properties
  int number;                 /**< @brief Bus number */
  char name[BUS_BUFFER_SIZE]; /**< @brief Bus name */

  // Times
  int num_periods;   /**< @brief Number of time periods. */

  // Voltage
  REAL* v_mag;        /**< @brief Voltage magnitude (p.u.) */
  REAL* v_ang;        /**< @brief Voltage angle (radians) */
  REAL* v_set;        /**< @brief Voltage magnitude set point (p.u.) */
  REAL v_max_reg;     /**< @brief Regulation maximum voltage magnitude (p.u.) */
  REAL v_min_reg;     /**< @brief Regulation minimum voltage magnitude (p.u.) */
  REAL v_max_norm;    /**< @brief Normal maximum voltage magnitude (p.u.) */
  REAL v_min_norm;    /**< @brief Normal minimum voltage magnitude (p.u.) */
  REAL v_max_emer;    /**< @brief Emergency maximum voltage magnitude (p.u.) */
  REAL v_min_emer;    /**< @brief Emergency minimum voltage magnitude (p.u.) */

  // Flags
  BOOL slack;        /**< @brief Flag for indicating the the bus is a slack bus */
  char fixed;        /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;      /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;       /**< @brief Flags for indicating which control adjustments should be sparse */
  char vars;         /**< @brief Flags for indicating which quantities should be treated as variables */

  // Price
  REAL* price;        /**< @brief Energy price at bus ($/ (hr p.u.)) */

  // Indices
  int index;          /**< @brief Bus index */
  int* index_v_mag;   /**< @brief Voltage magnitude index */
  int* index_v_ang;   /**< @brief Voltage angle index */

  // Components
  Gen* gen;            /**< @brief List of generators connected to bus */
  Gen* reg_gen;        /**< @brief List of generators regulating the voltage magnitude of bus */
  Load* load;          /**< @brief List of loads connected to bus */
  Shunt* shunt;        /**< @brief List of shunt devices connected to bus */
  Shunt* reg_shunt;    /**< @brief List of shunt devices regulating the voltage magnitude of bus */
  Branch* branch_k;    /**< @brief List of branches having this bus on the "k" side */
  Branch* branch_m;    /**< @brief List of branches having this bus on the "m" side */
  Branch* reg_tran;    /**< @brief List of transformers regulating the voltage magnitude of bus */
  Vargen* vargen;      /**< @brief List of variable generators connected to bus */
  Bat* bat;            /**< @brief List of batteries connected to bus */

  // Sensitivities
  REAL* sens_P_balance;      /**< @brief Sensitivity of active power balance */
  REAL* sens_Q_balance;      /**< @brief Sensitivity of reactive power balance */
  REAL* sens_v_mag_u_bound;  /**< @brief Sensitivity of voltage magnitude upper bound */
  REAL* sens_v_mag_l_bound;  /**< @brief Sensitivity of voltage magnitude lower bound */
  REAL* sens_v_ang_u_bound;  /**< @brief Sensitivity of voltage angle upper bound */
  REAL* sens_v_ang_l_bound;  /**< @brief Sensitivity of voltage angle lower bound */
  REAL* sens_v_reg_by_gen;   /**< @brief Sensitivity of voltage regulation by generator */
  REAL* sens_v_reg_by_tran;  /**< @brief Sensitivity of voltage regulation by transformer */
  REAL* sens_v_reg_by_shunt; /**< @brief Sensitivity of voltage regulation by shunt device */

  // Mismatches
  REAL* P_mis; /**< @brief Active power mismatch (p.u. system base power) */
  REAL* Q_mis; /**< @brief Reactive power mismatch (p.u. system base power) */

  // Hash
  UT_hash_handle hh_number; /**< @brief Handle for bus hash table based on numbers */
  UT_hash_handle hh_name;   /**< @brief Handle for bus hash table based on names */

  // List
  Bus* next; /**< @brief List of buses */
};

void BUS_add_gen(Bus* bus, Gen* gen) {
  if (bus)
    bus->gen = GEN_list_add(bus->gen,gen);
}

void BUS_del_gen(Bus* bus, Gen* gen) {
  if (bus)
    bus->gen = GEN_list_del(bus->gen,gen);
}

void BUS_add_load(Bus* bus, Load* load) {
  if (bus)
    bus->load = LOAD_list_add(bus->load,load);
}

void BUS_add_reg_gen(Bus* bus, Gen* reg_gen) {
  if (bus)
    bus->reg_gen = GEN_list_reg_add(bus->reg_gen,reg_gen);
}

void BUS_del_reg_gen(Bus* bus, Gen* reg_gen) {
  if (bus)
    bus->reg_gen = GEN_list_reg_del(bus->reg_gen,reg_gen);
}

void BUS_add_reg_tran(Bus* bus, Branch* reg_tran) {
  if (bus)
    bus->reg_tran = BRANCH_list_reg_add(bus->reg_tran,reg_tran);
}

void BUS_del_reg_tran(Bus* bus, Branch* reg_tran) {
  if (bus)
    bus->reg_tran = BRANCH_list_reg_del(bus->reg_tran,reg_tran);
}

void BUS_add_shunt(Bus* bus, Shunt* shunt) {
  if (bus)
    bus->shunt = SHUNT_list_add(bus->shunt,shunt);
}

void BUS_add_reg_shunt(Bus* bus, Shunt* reg_shunt) {
  if (bus)
    bus->reg_shunt = SHUNT_list_reg_add(bus->reg_shunt,reg_shunt);
}

void BUS_add_branch_k(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_k = BRANCH_list_k_add(bus->branch_k,branch);
}

void BUS_del_branch_k(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_k = BRANCH_list_k_del(bus->branch_k,branch);
}

void BUS_add_branch_m(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_m = BRANCH_list_m_add(bus->branch_m,branch);
}

void BUS_del_branch_m(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_m = BRANCH_list_m_del(bus->branch_m,branch);
}

void BUS_add_vargen(Bus* bus, Vargen* gen) {
  if (bus)
    bus->vargen = VARGEN_list_add(bus->vargen,gen);
}

void BUS_add_bat(Bus* bus, Bat* bat) {
  if (bus)
    bus->bat = BAT_list_add(bus->bat,bat);
}

BOOL BUS_array_check(Bus* bus_array, int size, BOOL verbose) {
  int i;
  BOOL bus_ok = TRUE;
  if (!bus_array)
    return bus_ok;
  for (i = 0; i < size; i++)
    bus_ok &= BUS_check(&(bus_array[i]),verbose);
  return bus_ok;
}

void BUS_array_del(Bus* bus_array, int size) {
  int i;
  Bus* bus;
  if (bus_array) {
    for (i = 0; i < size; i++) {
      bus = &(bus_array[i]);
      free(bus->v_mag);
      free(bus->v_ang);
      free(bus->v_set);
      free(bus->price);
      free(bus->index_v_mag);
      free(bus->index_v_ang);
      free(bus->sens_P_balance);
      free(bus->sens_Q_balance);
      free(bus->sens_v_mag_u_bound);
      free(bus->sens_v_mag_l_bound);
      free(bus->sens_v_ang_u_bound);
      free(bus->sens_v_ang_l_bound);
      free(bus->sens_v_reg_by_gen);
      free(bus->sens_v_reg_by_tran);
      free(bus->sens_v_reg_by_shunt);
      free(bus->P_mis);
      free(bus->Q_mis);
    }
    free(bus_array);
  }
}

void* BUS_array_get(void* bus_array, int index) {
  if (bus_array)
    return (void*)&(((Bus*)bus_array)[index]);
  else
    return NULL;
}

Bus* BUS_array_new(int size, int num_periods) {
  int i;
  if (num_periods > 0) {
    Bus* bus_array = (Bus*)malloc(sizeof(Bus)*size);
    for (i = 0; i < size; i++) {
      BUS_init(&(bus_array[i]),num_periods);
      BUS_set_index(&(bus_array[i]),i);
    }
    return bus_array;
  }
  else
    return NULL;
}

void BUS_array_show(Bus* bus_array, int size, int t) {
  int i;
  if (bus_array) {
    for (i = 0; i < size; i++)
      BUS_show(&(bus_array[i]),t);
  }
}

void BUS_array_get_max_mismatches(Bus* bus_array, int size, REAL* P, REAL* Q, int t) {
  int i;
  if (bus_array && t >= 0 && t < bus_array->num_periods) {
    for (i = 0; i < size; i++) {
      if (fabs(bus_array[i].P_mis[t]) > *P)
	*P = fabs(bus_array[i].P_mis[t]);
      if (fabs(bus_array[i].Q_mis[t]) > *Q)
	*Q = fabs(bus_array[i].Q_mis[t]);
    }
  }
}

BOOL BUS_check(Bus* bus, BOOL verbose) {

  // Local variables
  BOOL bus_ok = TRUE;

  // Null
  if (!bus) {
    if (verbose)
      fprintf(stderr,"NULL bus\n");
    return FALSE;
  }

  // Regulation voltage limits
  if (bus->v_min_reg > bus->v_max_reg) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"bad bus regulation voltage limits\n");
  }

  // Normal voltage violation limits
  if (bus->v_min_norm > bus->v_max_norm) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"bad bus normal voltage limits\n");
  }

  // Emergency voltage violation limits
  if (bus->v_min_emer > bus->v_max_emer) {
  bus_ok = FALSE;
  if (verbose)
    fprintf(stderr,"bad bus emergency voltage limits\n");
  }

  // Reg gen number
  if (BUS_is_regulated_by_gen(bus) && BUS_get_num_reg_gens(bus) < 1) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"reg-by-gen bus has no regulating generators\n");
  }

  // Reg tran
  if (BUS_is_regulated_by_tran(bus) && BUS_get_num_reg_trans(bus) < 1) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"reg-by-tran bus has no regulating transformer\n");
  }

  // Overall
  return bus_ok;

}

void BUS_clear_flags(Bus* bus, char flag_type) {
  if (bus) {
    if (flag_type == FLAG_VARS)
      bus->vars = 0x00;
    else if (flag_type == FLAG_BOUNDED)
      bus->bounded = 0x00;
    else if (flag_type == FLAG_FIXED)
      bus->fixed = 0x00;
    else if (flag_type == FLAG_SPARSE)
      bus->sparse = 0x00;
  }
}

void BUS_clear_sensitivities(Bus* bus) {
  int t;
  if (bus) {
    for (t = 0; t < bus->num_periods; t++) {
      bus->sens_P_balance[t] = 0;
      bus->sens_Q_balance[t] = 0;
      bus->sens_v_mag_u_bound[t] = 0;
      bus->sens_v_mag_l_bound[t] = 0;
      bus->sens_v_ang_u_bound[t] = 0;
      bus->sens_v_ang_l_bound[t] = 0;
      bus->sens_v_reg_by_gen[t] = 0;
      bus->sens_v_reg_by_tran[t] = 0;
      bus->sens_v_reg_by_shunt[t] = 0;
    }
  }
}

void BUS_clear_mismatches(Bus* bus) {
  int t;
  if (bus) {
    for (t = 0; t < bus->num_periods; t++) {
      bus->P_mis[t] = 0;
      bus->Q_mis[t] = 0;
    }
  }
}

void BUS_clear_vargen(Bus* bus) {
  if (bus)
    bus->vargen = NULL;
}

void BUS_clear_bat(Bus* bus) {
  if (bus)
    bus->bat = NULL;
}

char BUS_get_obj_type(void* bus) {
  if (bus)
    return OBJ_BUS;
  else
    return OBJ_UNKNOWN;
}

int BUS_get_degree(Bus* bus) {
  if (bus)
    return BRANCH_list_k_len(bus->branch_k)+BRANCH_list_m_len(bus->branch_m);
  else
    return 0;
}

REAL BUS_get_price(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->price[t];
  else
    return 0;
}

int BUS_get_index(Bus* bus) {
  if (bus)
    return bus->index;
  else
    return 0;
}

int BUS_get_index_v_mag(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->index_v_mag[t];
  else
    return 0;
}

int BUS_get_index_v_ang(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->index_v_ang[t];
  else
    return 0;
}

int BUS_get_index_P(Bus* bus) {
  if (bus)
    return 2*bus->index;
  else
    return 0;
}

int BUS_get_index_Q(Bus* bus) {
  if (bus)
    return 2*bus->index+1;
  else
    return 0;
}

Bus* BUS_get_next(Bus* bus) {
  if (bus)
    return bus->next;
  else
    return NULL;
}

int BUS_get_number(Bus* bus) {
  if (bus)
    return bus->number;
  else
    return 0;
}

char* BUS_get_name(Bus* bus) {
  if (bus)
    return bus->name;
  else
    return NULL;
}

int BUS_get_num_periods(Bus* bus) {
  if (bus)
    return bus->num_periods;
  else
    return 0;
}

int BUS_get_num_gens(Bus* bus) {
  if (bus)
    return GEN_list_len(bus->gen);
  else
    return 0;
}

int BUS_get_num_loads(Bus* bus) {
  if (bus)
    return LOAD_list_len(bus->load);
  else
    return 0;
}

int BUS_get_num_reg_gens(Bus* bus) {
  if (bus)
    return GEN_list_reg_len(bus->reg_gen);
  else
    return 0;
}

int BUS_get_num_reg_trans(Bus* bus) {
  if (bus)
    return BRANCH_list_reg_len(bus->reg_tran);
  else
    return 0;
}

int BUS_get_num_shunts(Bus* bus) {
  if (bus)
    return SHUNT_list_len(bus->shunt);
  else
    return 0;
}

int BUS_get_num_reg_shunts(Bus* bus) {
  if (bus)
    return SHUNT_list_reg_len(bus->reg_shunt);
  else
    return 0;
}

int BUS_get_num_vargens(Bus* bus) {
  if (bus)
    return VARGEN_list_len(bus->vargen);
  else
    return 0;
}

int BUS_get_num_bats(Bus* bus) {
  if (bus)
    return BAT_list_len(bus->bat);
  else
    return 0;
}

Gen* BUS_get_gen(Bus* bus) {
  if (bus)
    return bus->gen;
  else
    return NULL;
}

Load* BUS_get_load(Bus* bus) {
  if (bus)
    return bus->load;
  else
    return NULL;
}

Gen* BUS_get_reg_gen(Bus* bus) {
  if (bus)
    return bus->reg_gen;
  else
    return NULL;
}

Branch* BUS_get_reg_tran(Bus* bus) {
  if (bus)
    return bus->reg_tran;
  else
    return NULL;
}

Shunt* BUS_get_shunt(Bus* bus) {
  if (bus)
    return bus->shunt;
  else
    return NULL;
}

Shunt* BUS_get_reg_shunt(Bus* bus) {
  if (bus)
    return bus->reg_shunt;
  else
    return NULL;
}

Branch* BUS_get_branch_k(Bus* bus) {
  if (bus)
    return bus->branch_k;
  else
    return NULL;
}

Branch* BUS_get_branch_m(Bus* bus) {
  if (bus)
    return bus->branch_m;
  else
    return NULL;
}

Vargen* BUS_get_vargen(Bus* bus) {
  if (bus)
    return bus->vargen;
  else
    return NULL;
}

Bat* BUS_get_bat(Bus* bus) {
  if (bus)
    return bus->bat;
  else
    return NULL;
}

REAL BUS_get_P_mis(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->P_mis[t];
  else
    return 0;
}

REAL BUS_get_Q_mis(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->Q_mis[t];
  else
    return 0;
}

REAL BUS_get_total_gen_P(Bus* bus, int t) {
  Gen* gen;
  REAL P = 0;
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    P += GEN_get_P(gen,t);
  return P;
}

REAL BUS_get_total_gen_Q(Bus* bus, int t) {
  Gen* gen;
  REAL Q = 0;
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    Q += GEN_get_Q(gen,t);
  return Q;
}

REAL BUS_get_total_gen_Q_max(Bus* bus) {
  Gen* gen;
  REAL Qmax = 0;
  if (!bus)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    Qmax += GEN_get_Q_max(gen);
  return Qmax;
}

REAL BUS_get_total_gen_Q_min(Bus* bus) {
  Gen* gen;
  REAL Qmin = 0;
  if (!bus)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    Qmin += GEN_get_Q_min(gen);
  return Qmin;
}

REAL BUS_get_total_reg_gen_Qmax(Bus* bus) {
  Gen* gen;
  REAL Qmax = 0;
  if (!bus)
    return 0;
  for (gen = bus->reg_gen; gen != NULL; gen = GEN_get_reg_next(gen))
    Qmax += GEN_get_Q_max(gen);
  return Qmax;
}

REAL BUS_get_total_reg_gen_Qmin(Bus* bus) {
  Gen* gen;
  REAL Qmin = 0;
  if (!bus)
    return 0;
  for (gen = bus->reg_gen; gen != NULL; gen = GEN_get_reg_next(gen))
    Qmin += GEN_get_Q_min(gen);
  return Qmin;
}

REAL BUS_get_total_load_P(Bus* bus, int t) {
  Load* load;
  REAL P = 0;
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  for (load = bus->load; load != NULL; load = LOAD_get_next(load))
    P += LOAD_get_P(load,t);
  return P;
}

REAL BUS_get_total_load_Q(Bus* bus, int t) {
  Load* load;
  REAL Q = 0;
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  for (load = bus->load; load != NULL; load = LOAD_get_next(load))
    Q += LOAD_get_Q(load,t);
  return Q;
}

REAL BUS_get_total_shunt_g(Bus* bus) {
  Shunt* shunt;
  REAL g = 0;
  if (!bus)
    return 0;
  for (shunt = bus->shunt; shunt != NULL; shunt = SHUNT_get_next(shunt))
    g += SHUNT_get_g(shunt);
  return g;
}

REAL BUS_get_total_shunt_b(Bus* bus, int t) {
  Shunt* shunt;
  REAL b = 0;
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  for (shunt = bus->shunt; shunt != NULL; shunt = SHUNT_get_next(shunt))
    b += SHUNT_get_b(shunt,t);
  return b;
}

REAL BUS_get_v_mag(Bus* bus, int t) {
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  else
    return bus->v_mag[t];
}

REAL BUS_get_v_ang(Bus* bus, int t) {
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  else
    return bus->v_ang[t];
}

REAL BUS_get_v_set(Bus* bus, int t) {
  if (!bus || t < 0 || t >= bus->num_periods)
    return 0;
  else
    return bus->v_set[t];
}

REAL BUS_get_v_max_reg(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_max_reg;
}

REAL BUS_get_v_min_reg(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_min_reg;
}

REAL BUS_get_v_max_norm(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_max_norm;
}

REAL BUS_get_v_min_norm(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_min_norm;
}

REAL BUS_get_v_max_emer(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_max_emer;
}

REAL BUS_get_v_min_emer(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_min_emer;
}

void BUS_get_var_values(Bus* bus, Vec* values, int code) {

  // Local vars
  int t;

  // No bus
  if (!bus)
    return;

  for (t = 0; t < bus->num_periods; t++) {
    
    // Voltage magnitude
    if (bus->vars & BUS_VAR_VMAG) {
      switch (code) {

      case UPPER_LIMITS:
	if (bus->bounded & BUS_VAR_VMAG)
	  VEC_set(values,bus->index_v_mag[t],bus->v_max_norm);
	else
	  VEC_set(values,bus->index_v_mag[t],BUS_INF_V_MAG);
	break;

      case LOWER_LIMITS:
	if (bus->bounded & BUS_VAR_VMAG)
	  VEC_set(values,bus->index_v_mag[t],bus->v_min_norm);
	else
	  VEC_set(values,bus->index_v_mag[t],-BUS_INF_V_MAG);
	break;

      default:
	VEC_set(values,bus->index_v_mag[t],bus->v_mag[t]);
      }
    }

    // Voltage angle
    if (bus->vars & BUS_VAR_VANG) {
      switch(code) {
	
      case UPPER_LIMITS:
	VEC_set(values,bus->index_v_ang[t],BUS_INF_V_ANG);
	break;

      case LOWER_LIMITS:
	VEC_set(values,bus->index_v_ang[t],-BUS_INF_V_ANG);
	break;

      default:
	VEC_set(values,bus->index_v_ang[t],bus->v_ang[t]);
      }
    }
  }
}

int BUS_get_num_vars(void* vbus, unsigned char var, int t_start, int t_end) {

  // Local vars
  Bus* bus = (Bus*)vbus;
  int num_vars = 0;
  int dt;

  // Cheks
  if (!bus)
    return 0;
  if (t_start < 0)
    t_start = 0;
  if (t_end > bus->num_periods-1)
    t_end = bus->num_periods-1;

  // Num vars
  dt = t_end-t_start+1;
  if ((var & BUS_VAR_VMAG) && (bus->vars & BUS_VAR_VMAG))
    num_vars += dt;
  if ((var & BUS_VAR_VANG) && (bus->vars & BUS_VAR_VANG))
    num_vars += dt;
  return num_vars;
}

Vec* BUS_get_var_indices(void* vbus, unsigned char var, int t_start, int t_end) {

  // Local vars
  Bus* bus = (Bus*)vbus;
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
  indices = VEC_new(BUS_get_num_vars(vbus,var,t_start,t_end));
  if ((var & BUS_VAR_VMAG) && (bus->vars & BUS_VAR_VMAG)) { // v mag
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,bus->index_v_mag[t]);
      offset++;
    }
  }
  if ((var & BUS_VAR_VANG) && (bus->vars & BUS_VAR_VANG)) { // v ang
    for (t = t_start; t <= t_end; t++) {
      VEC_set(indices,offset,bus->index_v_ang[t]);
      offset++;
    }
  }
  return indices;
}

REAL BUS_get_sens_P_balance(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_P_balance[t];
  else
    return 0;
}

REAL BUS_get_sens_Q_balance(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_Q_balance[t];
  else
    return 0;
}

REAL BUS_get_sens_v_mag_u_bound(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_mag_u_bound[t];
  else
    return 0;
}

REAL BUS_get_sens_v_mag_l_bound(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_mag_l_bound[t];
  else
    return 0;
}

REAL BUS_get_sens_v_ang_u_bound(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_ang_u_bound[t];
  else
    return 0;
}

REAL BUS_get_sens_v_ang_l_bound(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_ang_l_bound[t];
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_gen(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_reg_by_gen[t];
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_tran(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_reg_by_tran[t];
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_shunt(Bus* bus, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    return bus->sens_v_reg_by_shunt[t];
  else
    return 0;
}

REAL BUS_get_largest_sens(Bus* bus, int t) {
  REAL sens = 0;
  if (bus && t >= 0 && t < bus->num_periods) {

    if (fabs(bus->sens_P_balance[t]) >= fabs(sens))
      sens = bus->sens_P_balance[t];

    if (fabs(bus->sens_Q_balance[t]) >= fabs(sens))
      sens = bus->sens_Q_balance[t];

    if (fabs(bus->sens_v_mag_u_bound[t]) >= fabs(sens))
      sens = bus->sens_v_mag_u_bound[t];

    if (fabs(bus->sens_v_mag_l_bound[t]) >= fabs(sens))
      sens = bus->sens_v_mag_l_bound[t];

    if (fabs(bus->sens_v_ang_u_bound[t]) >= fabs(sens))
      sens = bus->sens_v_ang_u_bound[t];

    if (fabs(bus->sens_v_ang_l_bound[t]) >= fabs(sens))
      sens = bus->sens_v_ang_l_bound[t];

    if (fabs(bus->sens_v_reg_by_gen[t]) >= fabs(sens))
      sens = bus->sens_v_reg_by_gen[t];

    if (fabs(bus->sens_v_reg_by_tran[t]) >= fabs(sens))
      sens = bus->sens_v_reg_by_tran[t];

    if (fabs(bus->sens_v_reg_by_shunt[t]) >= fabs(sens))
      sens = bus->sens_v_reg_by_shunt[t];
  }
  return sens;
}

int BUS_get_largest_sens_type(Bus* bus, int t) {
  REAL sens = 0;
  int type = BUS_SENS_P_BALANCE;
  if (bus && t >= 0 && t < bus->num_periods) {

    if (fabs(bus->sens_P_balance[t]) >= fabs(sens)) {
      sens = bus->sens_P_balance[t];
      type = BUS_SENS_P_BALANCE;
    }

    if (fabs(bus->sens_Q_balance[t]) >= fabs(sens)) {
      sens = bus->sens_Q_balance[t];
      type = BUS_SENS_Q_BALANCE;
    }

    if (fabs(bus->sens_v_mag_u_bound[t]) >= fabs(sens)) {
      sens = bus->sens_v_mag_u_bound[t];
      type = BUS_SENS_V_MAG_U_BOUND;
    }

    if (fabs(bus->sens_v_mag_l_bound[t]) >= fabs(sens)) {
      sens = bus->sens_v_mag_l_bound[t];
      type = BUS_SENS_V_MAG_L_BOUND;
    }

    if (fabs(bus->sens_v_ang_u_bound[t]) >= fabs(sens)) {
      sens = bus->sens_v_ang_u_bound[t];
      type = BUS_SENS_V_ANG_U_BOUND;
    }

    if (fabs(bus->sens_v_ang_l_bound[t]) >= fabs(sens)) {
      sens = bus->sens_v_ang_l_bound[t];
      type = BUS_SENS_V_ANG_L_BOUND;
    }

    if (fabs(bus->sens_v_reg_by_gen[t]) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_gen[t];
      type = BUS_SENS_V_REG_BY_GEN;
    }

    if (fabs(bus->sens_v_reg_by_tran[t]) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_tran[t];
      type = BUS_SENS_V_REG_BY_TRAN;
    }

    if (fabs(bus->sens_v_reg_by_shunt[t]) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_shunt[t];
      type = BUS_SENS_V_REG_BY_SHUNT;
    }
  }
  return type;
}

REAL BUS_get_largest_mis(Bus* bus, int t) {
  REAL mis = 0;
  if (bus && t >= 0 && t < bus->num_periods) {
    if (fabs(bus->P_mis[t]) >= fabs(mis))
      mis = bus->P_mis[t];
    if (fabs(bus->Q_mis[t]) >= fabs(mis))
      mis = bus->Q_mis[t];
  }
  return mis;
}

int BUS_get_largest_mis_type(Bus* bus, int t) {
  REAL mis = 0;
  int type = BUS_MIS_ACTIVE;
  if (bus && t >= 0 && t < bus->num_periods) {
    if (fabs(bus->P_mis[t]) >= fabs(mis)) {
      mis = bus->P_mis[t];
      type = BUS_MIS_ACTIVE;
    }
    if (fabs(bus->Q_mis[t]) >= fabs(mis)) {
      mis = bus->Q_mis[t];
      type = BUS_MIS_REACTIVE;
    }
  }
  return type;
}

REAL BUS_get_quantity(Bus* bus, int qtype, int t) {

  switch (qtype) {

  case BUS_SENS_LARGEST:
    return BUS_get_largest_sens(bus,t);

  case BUS_SENS_P_BALANCE:
    return BUS_get_sens_P_balance(bus,t);

  case BUS_SENS_Q_BALANCE:
    return BUS_get_sens_Q_balance(bus,t);

  case BUS_SENS_V_MAG_U_BOUND:
    return BUS_get_sens_v_mag_u_bound(bus,t);

  case BUS_SENS_V_MAG_L_BOUND:
    return BUS_get_sens_v_mag_l_bound(bus,t);

  case BUS_SENS_V_ANG_U_BOUND:
    return BUS_get_sens_v_ang_u_bound(bus,t);

  case BUS_SENS_V_ANG_L_BOUND:
    return BUS_get_sens_v_ang_l_bound(bus,t);

  case BUS_SENS_V_REG_BY_GEN:
    return BUS_get_sens_v_reg_by_gen(bus,t);

  case BUS_SENS_V_REG_BY_TRAN:
    return BUS_get_sens_v_reg_by_tran(bus,t);

  case BUS_SENS_V_REG_BY_SHUNT:
    return BUS_get_sens_v_reg_by_shunt(bus,t);

  case BUS_MIS_LARGEST:
    return BUS_get_largest_mis(bus,t);

  case BUS_MIS_ACTIVE:
    return BUS_get_P_mis(bus,t);

  case BUS_MIS_REACTIVE:
    return BUS_get_Q_mis(bus,t);

  default:
    return 0;
  }
}

char* BUS_get_json_string(Bus* bus, char* output) {
  
  // Local variables
  char temp[BUS_BUFFER_SIZE];
  char* output_start;
  BOOL resize;

  // No bus
  if (!bus)
    return NULL;

  // Output
  if (output)
    resize = FALSE;
  else {
    output = (char*)malloc(sizeof(char)*BUS_BUFFER_SIZE*BUS_NUM_JSON_FIELDS*bus->num_periods);
    resize = TRUE;
  }
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"number",bus->number,FALSE);
  JSON_str(temp,output,"name",bus->name,FALSE);
  JSON_int(temp,output,"num_periods",bus->num_periods,FALSE);
  JSON_array_float(temp,output,"v_mag",bus->v_mag,bus->num_periods,FALSE);
  JSON_array_float(temp,output,"v_ang",bus->v_ang,bus->num_periods,FALSE);
  JSON_array_float(temp,output,"v_set",bus->v_set,bus->num_periods,FALSE);
  JSON_float(temp,output,"v_max_reg",bus->v_max_reg,FALSE);
  JSON_float(temp,output,"v_min_reg",bus->v_min_reg,FALSE);
  JSON_float(temp,output,"v_max_norm",bus->v_max_norm,FALSE);
  JSON_float(temp,output,"v_min_norm",bus->v_min_norm,FALSE);
  JSON_float(temp,output,"v_max_emer",bus->v_max_emer,FALSE);
  JSON_float(temp,output,"v_min_emer",bus->v_min_emer,FALSE);
  JSON_bool(temp,output,"slack",bus->slack,FALSE);
  JSON_array_float(temp,output,"price",bus->price,bus->num_periods,FALSE);
  JSON_int(temp,output,"index",bus->index,FALSE);
  JSON_list_int(temp,output,"generators",bus,Gen,BUS_get_gen,GEN_get_index,GEN_get_next,FALSE);
  JSON_list_int(temp,output,"reg_generators",bus,Gen,BUS_get_reg_gen,GEN_get_index,GEN_get_reg_next,FALSE);
  JSON_list_int(temp,output,"loads",bus,Load,BUS_get_load,LOAD_get_index,LOAD_get_next,FALSE);
  JSON_list_int(temp,output,"shunts",bus,Shunt,BUS_get_shunt,SHUNT_get_index,SHUNT_get_next,FALSE);
  JSON_list_int(temp,output,"reg_shunts",bus,Shunt,BUS_get_reg_shunt,SHUNT_get_index,SHUNT_get_reg_next,FALSE);
  JSON_list_int(temp,output,"branches_k",bus,Branch,BUS_get_branch_k,BRANCH_get_index,BRANCH_get_next_k,FALSE);
  JSON_list_int(temp,output,"branches_m",bus,Branch,BUS_get_branch_m,BRANCH_get_index,BRANCH_get_next_m,FALSE);
  JSON_list_int(temp,output,"reg_transformers",bus,Branch,BUS_get_reg_tran,BRANCH_get_index,BRANCH_get_reg_next,FALSE);
  JSON_list_int(temp,output,"var_generators",bus,Vargen,BUS_get_vargen,VARGEN_get_index,VARGEN_get_next,FALSE);
  JSON_list_int(temp,output,"batteries",bus,Bat,BUS_get_bat,BAT_get_index,BAT_get_next,TRUE);
  JSON_end(output);
  
  // Resize
  if (resize)
    output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

BOOL BUS_has_flags(void* vbus, char flag_type, unsigned char mask) {
  Bus* bus = (Bus*)vbus;
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

BOOL BUS_has_properties(void* vbus, char prop) {
  Bus* bus = (Bus*)vbus;
  if (!bus)
    return FALSE;
  if ((prop & BUS_PROP_SLACK) && !BUS_is_slack(bus))
    return FALSE;
  if ((prop & BUS_PROP_REG_BY_GEN) && !BUS_is_regulated_by_gen(bus))
    return FALSE;
  if ((prop & BUS_PROP_REG_BY_TRAN) && !BUS_is_regulated_by_tran(bus))
    return FALSE;
  if ((prop & BUS_PROP_REG_BY_SHUNT) && !BUS_is_regulated_by_shunt(bus))
    return FALSE;
  if ((prop & BUS_PROP_NOT_REG_BY_GEN) && BUS_is_regulated_by_gen(bus))
    return FALSE;
  if ((prop & BUS_PROP_NOT_SLACK) && BUS_is_slack(bus))
    return FALSE;
  return TRUE;
}

Bus* BUS_hash_number_add(Bus* bus_hash,Bus* bus) {
  HASH_ADD(hh_number,bus_hash,number,sizeof(int),bus);
  return bus_hash;
}

void BUS_hash_number_del(Bus* bus_hash) {
  while (bus_hash != NULL)
    HASH_DELETE(hh_number,bus_hash,bus_hash);
}

Bus* BUS_hash_number_find(Bus* bus_hash,int number) {
  Bus* bus;
  HASH_FIND(hh_number,bus_hash,&number,sizeof(int),bus);
  return bus;
}

int BUS_hash_number_len(Bus* bus_hash) {
  return HASH_CNT(hh_number,bus_hash);
}

Bus* BUS_hash_name_add(Bus* bus_hash, Bus* bus) {
  HASH_ADD(hh_name,bus_hash,name[0],strlen(bus->name),bus);
  return bus_hash;
}

void BUS_hash_name_del(Bus* bus_hash) {
  while (bus_hash != NULL)
    HASH_DELETE(hh_name,bus_hash,bus_hash);
}

Bus* BUS_hash_name_find(Bus* bus_hash, char* name) {
  Bus* bus;
  HASH_FIND(hh_name,bus_hash,name,(unsigned)strlen(name),bus);
  return bus;
}

int BUS_hash_name_len(Bus* bus_hash) {
  return HASH_CNT(hh_name,bus_hash);
}

void BUS_init(Bus* bus, int num_periods) {

  // Local vars
  int i;
  int T;
  int t;

  // No bus
  if (!bus)
    return;

  T = num_periods;
  bus->num_periods = num_periods;

  bus->number = 0;
  for (i = 0; i < BUS_BUFFER_SIZE; i++)
    bus->name[i] = 0;

  bus->index = 0;

  bus->v_max_reg = BUS_DEFAULT_V_MAX;
  bus->v_min_reg = BUS_DEFAULT_V_MIN;
  bus->v_max_norm = BUS_DEFAULT_V_MAX;
  bus->v_min_norm = BUS_DEFAULT_V_MIN;
  bus->v_max_emer = BUS_DEFAULT_V_MAX;
  bus->v_min_emer = BUS_DEFAULT_V_MIN;

  bus->slack = FALSE;
  bus->fixed = 0x00;
  bus->bounded = 0x00;
  bus->sparse = 0x00;
  bus->vars = 0x00;

  bus->gen = NULL;
  bus->reg_gen = NULL;
  bus->reg_tran = NULL;
  bus->reg_shunt = NULL;
  bus->load = NULL;
  bus->shunt = NULL;
  bus->branch_k = NULL;
  bus->branch_m = NULL;
  bus->vargen = NULL;
  bus->bat = NULL;

  ARRAY_zalloc(bus->v_mag,REAL,T);
  ARRAY_zalloc(bus->v_ang,REAL,T);
  ARRAY_zalloc(bus->v_set,REAL,T);

  ARRAY_zalloc(bus->price,REAL,T);

  ARRAY_zalloc(bus->index_v_mag,int,T);
  ARRAY_zalloc(bus->index_v_ang,int,T);

  ARRAY_zalloc(bus->sens_P_balance,REAL,T);
  ARRAY_zalloc(bus->sens_Q_balance,REAL,T);
  ARRAY_zalloc(bus->sens_v_mag_u_bound,REAL,T);
  ARRAY_zalloc(bus->sens_v_mag_l_bound,REAL,T);
  ARRAY_zalloc(bus->sens_v_ang_u_bound,REAL,T);
  ARRAY_zalloc(bus->sens_v_ang_l_bound,REAL,T);
  ARRAY_zalloc(bus->sens_v_reg_by_gen,REAL,T);
  ARRAY_zalloc(bus->sens_v_reg_by_tran,REAL,T);
  ARRAY_zalloc(bus->sens_v_reg_by_shunt,REAL,T);

  ARRAY_zalloc(bus->P_mis,REAL,T);
  ARRAY_zalloc(bus->Q_mis,REAL,T);

  for (t = 0; t < bus->num_periods; t++) {
    bus->v_mag[t] = 1.;
    bus->v_set[t] = 1.;
  }

  bus->next = NULL;
}

void BUS_inject_P(Bus* bus, REAL P, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->P_mis[t] += P; // p.u.
}

void BUS_inject_Q(Bus* bus, REAL Q, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->Q_mis[t] += Q; // p.u.
}

BOOL BUS_is_equal(Bus* bus, Bus* other) {
  return bus == other;
}

BOOL BUS_is_regulated_by_gen(Bus* bus) {
  if (bus)
    return GEN_is_regulator(bus->reg_gen);
  else
    return FALSE;
}

BOOL BUS_is_regulated_by_tran(Bus* bus) {
  if (bus)
    return BRANCH_is_tap_changer_v(bus->reg_tran);
  else
    return FALSE;
}

BOOL BUS_is_regulated_by_shunt(Bus* bus) {
  if (bus)
    return SHUNT_is_switched_v(bus->reg_shunt);
  else
    return FALSE;
}

BOOL BUS_is_slack(Bus* bus) {
  if (bus)
    return bus->slack;
  else
    return FALSE;
}

Bus* BUS_list_add(Bus* bus_list, Bus* bus_new) {
  LIST_add(Bus,bus_list,bus_new,next);
  return bus_list;
}

Bus* BUS_list_add_sorting(Bus* bus_list, Bus* bus_new, int sort_by, int t) {

  // Local variables
  Bus* bus;
  Bus* bus_prev = NULL;
  REAL val_bus;
  REAL val_bus_new;

  for (bus = bus_list; bus != NULL; bus = bus->next) {

    // Get sorting quantities
    val_bus = BUS_get_quantity(bus,sort_by,t);
    val_bus_new = BUS_get_quantity(bus_new,sort_by,t);

    if (fabs(val_bus_new) >= fabs(val_bus)) // new bus comes before
      break;
    bus_prev = bus;
  }

  // New bus comes before bus
  bus_new->next = bus;
  if (bus_prev)
    bus_prev->next = bus_new;
  else
    bus_list = bus_new;
  return bus_list;
}

int BUS_list_len(Bus* bus_list) {
  int len;
  LIST_len(Bus,bus_list,next,len);
  return len;
}

Bus* BUS_new(int num_periods) {
  if (num_periods > 0) {
    Bus* bus = (Bus*)malloc(sizeof(Bus));
    BUS_init(bus,num_periods);
    return bus;
  }
  else
    return NULL;
}

void BUS_set_next(Bus* bus, Bus* next_bus) {
  if (bus) {
    bus->next = next_bus;
  }
}

void BUS_set_number(Bus* bus, int number) {
  if (bus)
    bus->number = number;
}

void BUS_set_name(Bus* bus, char* name) {
  if (bus)
    strncpy(bus->name,name,(size_t)(BUS_BUFFER_SIZE-1));
}

void BUS_set_price(Bus* bus, REAL price, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->price[t] = price;
}

void BUS_set_v_mag(Bus* bus, REAL v_mag, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->v_mag[t] = v_mag;
}

void BUS_set_v_ang(Bus* bus, REAL v_ang, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->v_ang[t] = v_ang;
}

void BUS_set_v_set(Bus* bus, REAL v_set, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->v_set[t] = v_set;
}

void BUS_set_v_max_reg(Bus* bus, REAL v_max_reg) {
  if (bus)
    bus->v_max_reg = v_max_reg;
}

void BUS_set_v_min_reg(Bus* bus, REAL v_min_reg) {
  if (bus)
    bus->v_min_reg = v_min_reg;
}

void BUS_set_v_max_norm(Bus* bus, REAL v_max_norm) {
  if (bus)
    bus->v_max_norm = v_max_norm;
}

void BUS_set_v_min_norm(Bus* bus, REAL v_min_norm) {
  if (bus)
    bus->v_min_norm = v_min_norm;
}

void BUS_set_v_max_emer(Bus* bus, REAL v_max_emer) {
  if (bus)
    bus->v_max_emer = v_max_emer;
}

void BUS_set_v_min_emer(Bus* bus, REAL v_min_emer) {
  if (bus)
    bus->v_min_emer = v_min_emer;
}

void BUS_set_slack(Bus* bus, BOOL slack) {
  if (bus)
    bus->slack = slack;
}

void BUS_set_index(Bus* bus, int index) {
  if (bus)
    bus->index = index;
}

int BUS_set_flags(void* vbus, char flag_type, unsigned char mask, int index) {

  // Local variables
  char* flags_ptr = NULL;
  Bus* bus = (Bus*)vbus;
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
  if (!((*flags_ptr) & BUS_VAR_VMAG) && (mask & BUS_VAR_VMAG)) { // voltage magnitude
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < bus->num_periods; t++)
	bus->index_v_mag[t] = index+t;
    }
    (*flags_ptr) |= BUS_VAR_VMAG;
    index += bus->num_periods;
  }
  if (!((*flags_ptr) & BUS_VAR_VANG) && (mask & BUS_VAR_VANG)) { // voltage angle
    if (flag_type == FLAG_VARS) {
      for (t = 0; t < bus->num_periods; t++)
	bus->index_v_ang[t] = index+t;
    }
    (*flags_ptr) |= BUS_VAR_VANG;
    index += bus->num_periods;
  }
  return index;
}

void BUS_set_var_values(Bus* bus, Vec* values) {

  // Local vars
  int t;

  // No bus
  if (!bus)
    return;

  // Time loop
  for (t = 0; t < bus->num_periods; t++) {

    // Voltage
    if (bus->vars & BUS_VAR_VMAG)      // voltage magnitude (p.u.)
      bus->v_mag[t] = VEC_get(values,bus->index_v_mag[t]);
    if (bus->vars & BUS_VAR_VANG)      // voltage angle (radians)
      bus->v_ang[t] = VEC_get(values,bus->index_v_ang[t]);
  }
}

void BUS_set_sens_P_balance(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_P_balance[t] = value;
}

void BUS_set_sens_Q_balance(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_Q_balance[t] = value;
}

void BUS_set_sens_v_mag_u_bound(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_mag_u_bound[t] = value;
}

void BUS_set_sens_v_mag_l_bound(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_mag_l_bound[t] = value;
}

void BUS_set_sens_v_ang_u_bound(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_ang_u_bound[t] = value;
}

void BUS_set_sens_v_ang_l_bound(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_ang_l_bound[t] = value;
}

void BUS_set_sens_v_reg_by_gen(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_reg_by_gen[t] = value;
}

void BUS_set_sens_v_reg_by_tran(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_reg_by_tran[t] = value;
}

void BUS_set_sens_v_reg_by_shunt(Bus* bus, REAL value, int t) {
  if (bus && t >= 0 && t < bus->num_periods)
    bus->sens_v_reg_by_shunt[t] = value;
}

void BUS_show(Bus* bus, int t) {
  printf("bus %d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	 BUS_get_number(bus),
	 BUS_is_slack(bus),
	 BUS_is_regulated_by_gen(bus),
	 BUS_get_num_gens(bus),
	 BUS_get_num_reg_gens(bus),
	 BUS_get_num_loads(bus),
	 BUS_get_num_shunts(bus));
}

void BUS_propagate_data_in_time(Bus* bus) {
  int t;
  if (bus) {
    for (t = 1; t < bus->num_periods; t++) {
      bus->v_mag[t] = bus->v_mag[0];
      bus->v_ang[t] = bus->v_ang[0];
      bus->v_set[t] = bus->v_set[0];
    }
  }
}
