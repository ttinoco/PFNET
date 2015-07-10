/** @file bus.c
 *  @brief This file defines the Bus data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/bus.h>

struct Bus {

  // Properties
  int number;        /**< @brief Bus number */

  // Voltage
  REAL v_mag;        /**< @brief Voltage magnitude (p.u.) */
  REAL v_ang;        /**< @brief Voltage angle (radians) */
  REAL v_set;        /**< @brief Voltage magnitude set point (p.u.) */
  REAL v_max;        /**< @brief Maximum voltage magnitude (p.u.) */
  REAL v_min;        /**< @brief Minimum voltage magnitude (p.u.) */

  // Flags
  BOOL slack;        /**< @brief Flag for indicating the the bus is a slack bus */
  char fixed;        /**< @brief Flags for indicating which quantities should be fixed to their current value */
  char bounded;      /**< @brief Flags for indicating which quantities should be bounded */
  char sparse;       /**< @brief Flags for indicating which control adjustments should be sparse */
  char vars;         /**< @brief Flags for indicating which quantities should be treated as variables */

  // Components
  Gen* gen;            /**< @brief List of generators connected to bus */
  Gen* reg_gen;        /**< @brief List of generators regulating the voltage magnitude of bus */
  Branch* reg_tran;    /**< @brief List of transformers regulating the voltage magnitude of bus */
  Load* load;          /**< @brief List of loads connected to bus */
  Shunt* shunt;        /**< @brief List of shunt devices connected to bus */
  Shunt* reg_shunt;    /**< @brief List of shunt devices regulating the voltage magnitude of bus */
  Branch* branch_from; /**< @brief List of branches having this bus on the "from" side */
  Branch* branch_to;   /**< @brief List of branches having this bus on the "to" side */
  
  // Indices
  int index;         /**< @brief Bus index */
  int index_v_mag;   /**< @brief Voltage magnitude index */
  int index_v_ang;   /**< @brief Voltage angle index */
  int index_y;       /**< @brief Voltage magnitude positive-deviation index */
  int index_z;       /**< @brief Voltage magnitude negative-deviation index */
  int index_vl;      /**< @brief Voltage magnitude violation of v_min */
  int index_vh;      /**< @brief Voltage magnitude violation of v_max */

  // Sensitivities
  REAL sens_P_balance;      /**< @brief Sensitivity of active power balance */
  REAL sens_Q_balance;      /**< @brief Sensitivity of reactive power balance */
  REAL sens_v_mag_u_bound;  /**< @brief Sensitivity of voltage magnitude upper bound */
  REAL sens_v_mag_l_bound;  /**< @brief Sensitivity of voltage magnitude lower bound */
  REAL sens_v_reg_by_gen;   /**< @brief Sensitivity of voltage regulation by generator */
  REAL sens_v_reg_by_tran;  /**< @brief Sensitivity of voltage regulation by transformer */
  REAL sens_v_reg_by_shunt; /**< @brief Sensitivity of voltage regulation by shunt device */

  // Mismatches
  REAL P_mis; /**< @brief Active power mismatch (p.u. system base power) */
  REAL Q_mis; /**< @brief Reactive power mismatch (p.u. system base power) */
  
  // Hash
  UT_hash_handle hh; /**< @brief Handle for bus hash table */

  // List
  Bus* next; /**< @brief List of buses */
};

void BUS_add_gen(Bus* bus, Gen* gen) {
  if (bus)
    bus->gen = GEN_list_add(bus->gen,gen);
}

void BUS_add_load(Bus* bus, Load* load) {
  if (bus)
    bus->load = LOAD_list_add(bus->load,load);
}

void BUS_add_reg_gen(Bus* bus, Gen* reg_gen) {
  if (bus)
    bus->reg_gen = GEN_list_reg_add(bus->reg_gen,reg_gen);
}

void BUS_add_reg_tran(Bus* bus, Branch* reg_tran) {
  if (bus)
    bus->reg_tran = BRANCH_list_reg_add(bus->reg_tran,reg_tran);
}

void BUS_add_shunt(Bus* bus, Shunt* shunt) {
  if (bus)
    bus->shunt = SHUNT_list_add(bus->shunt,shunt);
}

void BUS_add_reg_shunt(Bus* bus, Shunt* reg_shunt) {
  if (bus)
    bus->reg_shunt = SHUNT_list_reg_add(bus->reg_shunt,reg_shunt);
}

void BUS_add_branch_from(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_from = BRANCH_list_from_add(bus->branch_from,branch);
}

void BUS_add_branch_to(Bus* bus, Branch* branch) {
  if (bus)
    bus->branch_to = BRANCH_list_to_add(bus->branch_to,branch);
}

BOOL BUS_array_check(Bus* bus, int num, BOOL verbose) {
  int i;
  BOOL bus_ok = TRUE;
  for (i = 0; i < num; i++)
    bus_ok &= BUS_check(&(bus[i]),verbose);
  return bus_ok;
}

void* BUS_array_get(void* bus, int index) {
  if (bus)
    return (void*)&(((Bus*)bus)[index]);
  else
    return NULL;
}

Bus* BUS_array_new(int num) {
  int i;
  Bus* bus = (Bus*)malloc(sizeof(Bus)*num);
  for (i = 0; i < num; i++) {
    BUS_init(&(bus[i]));
    BUS_set_index(&(bus[i]),i);
  }
  return bus;
}

void BUS_array_show(Bus* bus, int num) {
  int i;
  if (bus) {
    for (i = 0; i < num; i++) 
      BUS_show(&(bus[i]));
  }
}

void BUS_array_get_max_mismatches(Bus* bus, int num, REAL* P, REAL* Q) {
  int i;
  if (bus) {
    for (i = 0; i < num; i++) {
      if (fabs(bus[i].P_mis) > *P)
	*P = fabs(bus[i].P_mis);
      if (fabs(bus[i].Q_mis) > *Q)
	*Q = fabs(bus[i].Q_mis);
    }
  }
}

BOOL BUS_check(Bus* bus, BOOL verbose) {

  // Local variables
  BOOL bus_ok = TRUE;
  REAL frac_sum;
  Gen* gen;

  // Null
  if (!bus) {
    if (verbose)
      fprintf(stderr,"NULL bus\n");
    return FALSE;    
  }

  // v limits
  if (bus->v_min > bus->v_max) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"bad bus v limits\n");
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
  if (bus) {
    bus->sens_P_balance = 0;
    bus->sens_Q_balance = 0;
    bus->sens_v_mag_u_bound = 0;
    bus->sens_v_mag_l_bound = 0;
    bus->sens_v_reg_by_gen = 0;
    bus->sens_v_reg_by_tran = 0;
    bus->sens_v_reg_by_shunt = 0;
  }
}

void BUS_clear_mismatches(Bus* bus) {
  if (bus) {
    bus->P_mis = 0;
    bus->Q_mis = 0;
  }
}

int BUS_get_degree(Bus* bus) {
  if (bus) {
    return (BRANCH_list_from_len(bus->branch_from) +
	    BRANCH_list_to_len(bus->branch_to));
  }
  else
    return 0;
}

int BUS_get_index(Bus* bus) {
  if (bus)
    return bus->index;
  else
    return 0;
}

int BUS_get_index_v_mag(Bus* bus) {
  if (bus)
    return bus->index_v_mag;
  else
    return 0;
}

int BUS_get_index_v_ang(Bus* bus) {
  if (bus)
    return bus->index_v_ang;
  else
    return 0;
}

int BUS_get_index_y(Bus* bus) {
  if (bus)
    return bus->index_y;
  else
    return 0;
}

int BUS_get_index_z(Bus* bus) {
  if (bus)
    return bus->index_z;
  else
    return 0;
}

int BUS_get_index_vl(Bus* bus) {
  if (bus)
    return bus->index_vl;
  else
    return 0;
}

int BUS_get_index_vh(Bus* bus) {
  if (bus)
    return bus->index_vh;
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

Branch* BUS_get_branch_from(Bus* bus) {
  if (bus)
    return bus->branch_from;
  else
    return NULL;
}

Branch* BUS_get_branch_to(Bus* bus) {
  if (bus)
    return bus->branch_to;
  else
    return NULL;
}

REAL BUS_get_P_mis(Bus* bus) {
  if (bus)
    return bus->P_mis;
  else
    return 0;
}

REAL BUS_get_Q_mis(Bus* bus) {
  if (bus)
    return bus->Q_mis;
  else
    return 0;
}

REAL BUS_get_total_gen_P(Bus* bus) {
  Gen* gen;
  REAL P = 0;
  if (!bus)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    P += GEN_get_P(gen);
  return P;
}

REAL BUS_get_total_gen_Q(Bus* bus) {
  Gen* gen;
  REAL Q = 0;
  if (!bus)
    return 0;
  for (gen = bus->gen; gen != NULL; gen = GEN_get_next(gen))
    Q += GEN_get_Q(gen);
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

REAL BUS_get_total_load_P(Bus* bus) {
  Load* load;
  REAL P = 0;
  if (!bus)
    return 0;
  for (load = bus->load; load != NULL; load = LOAD_get_next(load))
    P += LOAD_get_P(load);
  return P;
}

REAL BUS_get_total_load_Q(Bus* bus) {
  Load* load;
  REAL Q = 0;
  if (!bus)
    return 0;
  for (load = bus->load; load != NULL; load = LOAD_get_next(load))
    Q += LOAD_get_Q(load);
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

REAL BUS_get_total_shunt_b(Bus* bus) {
  Shunt* shunt;
  REAL b = 0;
  if (!bus)
    return 0;
  for (shunt = bus->shunt; shunt != NULL; shunt = SHUNT_get_next(shunt)) 
    b += SHUNT_get_b(shunt);
  return b;
}

REAL BUS_get_v_mag(Bus* bus) {
  if (!bus)
    return 0;
  else 
    return bus->v_mag;
}

REAL BUS_get_v_ang(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_ang;
}

REAL BUS_get_v_set(Bus* bus) {
  if (!bus)
    return 0;
  else
    return bus->v_set;
}

REAL BUS_get_v_max(Bus* bus) {
  if (!bus)
    return 0;
  else 
    return bus->v_max;
}

REAL BUS_get_v_min(Bus* bus) {
  if (!bus)
    return 0;
  else 
    return bus->v_min;
}

void BUS_get_var_values(Bus* bus, Vec* values) {

  // No bus
  if (!bus)
    return;

  // Voltage
  if (bus->vars & BUS_VAR_VMAG)      // voltage magnitude
    VEC_set(values,bus->index_v_mag,bus->v_mag);
  if (bus->vars & BUS_VAR_VANG)      // voltage angle
    VEC_set(values,bus->index_v_ang,bus->v_ang);
  if (bus->vars & BUS_VAR_VDEV) {
    if (bus->v_mag > bus->v_set) { // pos voltage mag deviation
      VEC_set(values,bus->index_y,bus->v_mag-bus->v_set); 
      VEC_set(values,bus->index_z,0.);
    }
    else {                         // neg voltage mag deviation
      VEC_set(values,bus->index_y,0.);
      VEC_set(values,bus->index_z,bus->v_set-bus->v_mag); 
    }    
  }
  if (bus->vars & BUS_VAR_VVIO) { // max min mag bound violations
    VEC_set(values,bus->index_vl,0.);
    VEC_set(values,bus->index_vh,0.);
  }
}

REAL BUS_get_sens_P_balance(Bus* bus) {
  if (bus)
    return bus->sens_P_balance;
  else
    return 0;
}

REAL BUS_get_sens_Q_balance(Bus* bus) {
  if (bus)
    return bus->sens_Q_balance;
  else
    return 0;
}

REAL BUS_get_sens_v_mag_u_bound(Bus* bus) {
  if (bus)
    return bus->sens_v_mag_u_bound;
  else
    return 0;
}

REAL BUS_get_sens_v_mag_l_bound(Bus* bus) {
  if (bus)
    return bus->sens_v_mag_l_bound;
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_gen(Bus* bus) {
  if (bus)
    return bus->sens_v_reg_by_gen;
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_tran(Bus* bus) {
  if (bus)
    return bus->sens_v_reg_by_tran;
  else
    return 0;
}

REAL BUS_get_sens_v_reg_by_shunt(Bus* bus) {
  if (bus)
    return bus->sens_v_reg_by_shunt;
  else
    return 0;
}

REAL BUS_get_largest_sens(Bus* bus) {
  REAL sens = 0;
  if (bus) {
    if (fabs(bus->sens_P_balance) >= fabs(sens))
      sens = bus->sens_P_balance;
    if (fabs(bus->sens_Q_balance) >= fabs(sens))
      sens = bus->sens_Q_balance;
    if (fabs(bus->sens_v_mag_u_bound) >= fabs(sens))
      sens = bus->sens_v_mag_u_bound;
    if (fabs(bus->sens_v_mag_l_bound) >= fabs(sens))
      sens = bus->sens_v_mag_l_bound;
    if (fabs(bus->sens_v_reg_by_gen) >= fabs(sens))
      sens = bus->sens_v_reg_by_gen;
    if (fabs(bus->sens_v_reg_by_tran) >= fabs(sens))
      sens = bus->sens_v_reg_by_tran;
    if (fabs(bus->sens_v_reg_by_shunt) >= fabs(sens))
      sens = bus->sens_v_reg_by_shunt;
  }
  return sens;
}

int BUS_get_largest_sens_type(Bus* bus) {
  REAL sens = 0;
  int type = BUS_SENS_P_BALANCE;
  if (bus) {
    if (fabs(bus->sens_P_balance) >= fabs(sens)) {
      sens = bus->sens_P_balance;
      type = BUS_SENS_P_BALANCE;
    }    
    if (fabs(bus->sens_Q_balance) >= fabs(sens)) {
      sens = bus->sens_Q_balance;
      type = BUS_SENS_Q_BALANCE;
    }
    if (fabs(bus->sens_v_mag_u_bound) >= fabs(sens)) {
      sens = bus->sens_v_mag_u_bound;
      type = BUS_SENS_V_MAG_U_BOUND;
    }
    if (fabs(bus->sens_v_mag_l_bound) >= fabs(sens)) {
      sens = bus->sens_v_mag_l_bound;
      type = BUS_SENS_V_MAG_L_BOUND;
    }
    if (fabs(bus->sens_v_reg_by_gen) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_gen;
      type = BUS_SENS_V_REG_BY_GEN;
    }
    if (fabs(bus->sens_v_reg_by_tran) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_tran;
      type = BUS_SENS_V_REG_BY_TRAN;
    }
    if (fabs(bus->sens_v_reg_by_shunt) >= fabs(sens)) {
      sens = bus->sens_v_reg_by_shunt;
      type = BUS_SENS_V_REG_BY_SHUNT;
    }
  }
  return type;
}

REAL BUS_get_largest_mis(Bus* bus) {
  REAL mis = 0;
  if (bus) {
    if (fabs(bus->P_mis) >= fabs(mis))
      mis = bus->P_mis;
    if (fabs(bus->Q_mis) >= fabs(mis))
      mis = bus->Q_mis;
  }
  return mis;
}

int BUS_get_largest_mis_type(Bus* bus) {
  REAL mis = 0;
  int type = BUS_MIS_ACTIVE;
  if (bus) {
    if (fabs(bus->P_mis) >= fabs(mis)) {
      mis = bus->P_mis;
      type = BUS_MIS_ACTIVE;
    }
    if (fabs(bus->Q_mis) >= fabs(mis)) {
      mis = bus->Q_mis;
      type = BUS_MIS_REACTIVE;
    }
  }
  return type;
}

REAL BUS_get_quantity(Bus* bus, int qtype) {
  
  switch (qtype) {
  case BUS_SENS_LARGEST:
    return BUS_get_largest_sens(bus);
  case BUS_SENS_P_BALANCE:
    return BUS_get_sens_P_balance(bus);
  case BUS_SENS_Q_BALANCE:
    return BUS_get_sens_Q_balance(bus);
  case BUS_SENS_V_MAG_U_BOUND:
    return BUS_get_sens_v_mag_u_bound(bus);
  case BUS_SENS_V_MAG_L_BOUND:
    return BUS_get_sens_v_mag_l_bound(bus);
  case BUS_SENS_V_REG_BY_GEN:
    return BUS_get_sens_v_reg_by_gen(bus);
  case BUS_SENS_V_REG_BY_TRAN:
    return BUS_get_sens_v_reg_by_tran(bus);
  case BUS_SENS_V_REG_BY_SHUNT:
    return BUS_get_sens_v_reg_by_shunt(bus);
  case BUS_MIS_LARGEST:
    return BUS_get_largest_mis(bus);
  case BUS_MIS_ACTIVE:
    return BUS_get_P_mis(bus);
  case BUS_MIS_REACTIVE:
    return BUS_get_Q_mis(bus);
  default:
    return 0;
  }
}

BOOL BUS_has_flags(Bus* bus, char flag_type, char mask) {
  if (bus) {
    if (flag_type == FLAG_VARS)
      return (bus->vars & mask);
    else if (flag_type == FLAG_BOUNDED)
      return (bus->bounded & mask);
    else if (flag_type == FLAG_FIXED)
      return (bus->fixed & mask);
    else if (flag_type == FLAG_SPARSE)
      return (bus->sparse & mask);
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

Bus* BUS_hash_add(Bus* bus_hash,Bus* bus) {
  HASH_ADD_INT(bus_hash,number,bus);
  return bus_hash;
}

void BUS_hash_del(Bus* bus_hash) {
  while (bus_hash != NULL)
    HASH_DEL(bus_hash,bus_hash);
}

Bus* BUS_hash_find(Bus* bus_hash,int number) {
  Bus* bus;
  HASH_FIND_INT(bus_hash,&number,bus);
  return bus;
}

int BUS_hash_len(Bus* bus_hash) {
  return HASH_COUNT(bus_hash);
}

void BUS_init(Bus* bus) {

  bus->number = 0;

  bus->v_mag = 1.;
  bus->v_ang = 0.;
  bus->v_set = 1.;
  bus->v_max = BUS_DEFAULT_V_MAX;
  bus->v_min = BUS_DEFAULT_V_MIN;

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
  bus->branch_from = NULL;
  bus->branch_to = NULL;

  bus->index = 0;
  bus->index_v_mag = 0;
  bus->index_v_ang = 0;
  bus->index_y = 0;
  bus->index_z = 0;
  bus->index_vl = 0;
  bus->index_vh = 0;

  bus->sens_P_balance = 0;
  bus->sens_Q_balance = 0;
  bus->sens_v_mag_u_bound = 0;
  bus->sens_v_mag_l_bound = 0;
  bus->sens_v_reg_by_gen = 0;
  bus->sens_v_reg_by_tran = 0;
  bus->sens_v_reg_by_shunt = 0;

  bus->P_mis = 0;
  bus->Q_mis = 0;

  bus->next = NULL;
}

void BUS_inject_P(Bus* bus, REAL P) {
  if (bus)
    bus->P_mis += P; // p.u.
}

void BUS_inject_Q(Bus* bus, REAL Q) {
  if (bus)
    bus->Q_mis += Q; // p.u.
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

Bus* BUS_list_add(Bus* bus_list, Bus* bus_new, int sort_by) {

  // Local variables
  Bus* bus;
  Bus* bus_prev = NULL;
  REAL val_bus;
  REAL val_bus_new;

  for (bus = bus_list; bus != NULL; bus = bus->next) {
    
    // Get sorting quantities
    val_bus = BUS_get_quantity(bus,sort_by);
    val_bus_new = BUS_get_quantity(bus_new,sort_by);

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

Bus* BUS_new(void) {
  Bus* bus = (Bus*)malloc(sizeof(Bus));
  BUS_init(bus);
  return bus;
}

void BUS_set_number(Bus* bus, int number) {
  if (bus)
    bus->number = number;
}

void BUS_set_v_mag(Bus* bus, REAL v_mag) {
  if (bus)
    bus->v_mag = v_mag;
}

void BUS_set_v_ang(Bus* bus, REAL v_ang) {
  if (bus)
    bus->v_ang = v_ang;
}

void BUS_set_v_set(Bus* bus, REAL v_set) {
  if (bus)
    bus->v_set = v_set;
}

void BUS_set_v_max(Bus* bus, REAL v_max) {
  if (bus)
    bus->v_max = v_max;
}

void BUS_set_v_min(Bus* bus, REAL v_min) {
  if (bus)
    bus->v_min = v_min;
}

void BUS_set_slack(Bus* bus, BOOL slack) {
  if (bus)
    bus->slack = slack;
}

void BUS_set_index(Bus* bus, int index) {
  if (bus)
    bus->index = index;
}

int BUS_set_flags(void* vbus, char flag_type, char mask, int index) {
  
  // Local variables
  char* flags_ptr = NULL;
  Bus* bus = (Bus*)vbus;
  
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
    if (flag_type == FLAG_VARS)
      bus->index_v_mag = index;
    (*flags_ptr) |= BUS_VAR_VMAG;
    index++;
  }
  if (!((*flags_ptr) & BUS_VAR_VANG) && (mask & BUS_VAR_VANG)) { // voltage angle
    if (flag_type == FLAG_VARS)
      bus->index_v_ang = index;
    (*flags_ptr) |= BUS_VAR_VANG;
    index++;
  }
  if (!((*flags_ptr) & BUS_VAR_VDEV) && (mask & BUS_VAR_VDEV)) { // voltage mag deviation
    if (flag_type == FLAG_VARS) {
      bus->index_y = index;
      bus->index_z = index+1;
    }
    (*flags_ptr) |= BUS_VAR_VDEV;
    index += 2;
  }
  if (!((*flags_ptr) & BUS_VAR_VVIO) && (mask & BUS_VAR_VVIO)) { // voltage mag max min bound violations
    if (flag_type == FLAG_VARS) {
      bus->index_vl = index;
      bus->index_vh = index+1;
    }
    (*flags_ptr) |= BUS_VAR_VVIO;
    index += 2;
  }
  return index;
}

void BUS_set_var_values(Bus* bus, Vec* values) {

  // No bus
  if (!bus)
    return;

  // Voltage
  if (bus->vars & BUS_VAR_VMAG)      // voltage magnitude (p.u.)
    bus->v_mag = VEC_get(values,bus->index_v_mag); 
  if (bus->vars & BUS_VAR_VANG)      // voltage angle (radians)
    bus->v_ang = VEC_get(values,bus->index_v_ang);
}

void BUS_set_sens_P_balance(Bus* bus, REAL value) {
  if (bus)
    bus->sens_P_balance = value;
}

void BUS_set_sens_Q_balance(Bus* bus, REAL value) {
  if (bus)
    bus->sens_Q_balance = value;
}

void BUS_set_sens_v_mag_u_bound(Bus* bus, REAL value) {
  if (bus)
    bus->sens_v_mag_u_bound = value;
}

void BUS_set_sens_v_mag_l_bound(Bus* bus, REAL value) {
  if (bus)
    bus->sens_v_mag_l_bound = value;
}

void BUS_set_sens_v_reg_by_gen(Bus* bus, REAL value) {
  if (bus)
    bus->sens_v_reg_by_gen = value;
}

void BUS_set_sens_v_reg_by_tran(Bus* bus, REAL value) {
  if (bus)
    bus->sens_v_reg_by_tran = value;
}

void BUS_set_sens_v_reg_by_shunt(Bus* bus, REAL value) {
  if (bus)
    bus->sens_v_reg_by_shunt = value;
}

void BUS_show(Bus *bus) {
  printf("bus %d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	 BUS_get_number(bus),
	 BUS_is_slack(bus),
	 BUS_is_regulated_by_gen(bus),
	 BUS_get_num_gens(bus),
	 BUS_get_num_reg_gens(bus),
	 BUS_get_num_loads(bus),
	 BUS_get_num_shunts(bus));
}
