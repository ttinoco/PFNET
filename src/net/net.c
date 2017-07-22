/** @file net.c
 *  @brief This file defines the Net data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/net.h>
#include <pfnet/array.h>

struct Net {

  // Error
  BOOL error_flag;                    /**< @brief Error flag. */
  char error_string[NET_BUFFER_SIZE]; /**< @brief Error string. */

  // Num periods
  int num_periods;   /**< @brief Number of time periods. */

  // Output
  char output_string[NET_BUFFER_SIZE]; /**< @brief Output string. */

  // Components
  Bus* bus;             /**< @brief Bus array. */
  Branch* branch;       /**< @brief Branch array. */
  Gen* gen;             /**< @brief Gen array. */
  Load* load;           /**< @brief Load array. */
  Shunt* shunt;         /**< @brief Shunt array. */
  Vargen* vargen;       /**< @brief Vargen array. */
  Bat* bat;             /**< @brief Bat array. */

  // Hash tables
  Bus* bus_hash_number;     /**< @brief Bus hash table indexed by bus numbers. */
  Bus* bus_hash_name;       /**< @brief Bus hash table indexed by bus names. */
  Vargen* vargen_hash_name; /**< @brief Vargen hash table indexed by vargen names. */

  // Number of components
  int num_buses;     /**< @brief Number of buses (size of Bus array). */
  int num_branches;  /**< @brief Number of branches (size of Branch array). */
  int num_gens;      /**< @brief Number of generators (size of Gen array). */
  int num_loads;     /**< @brief Number of loads (size of Load array). */
  int num_shunts;    /**< @brief Number of shunts (size of Shunt array). */
  int num_vargens;   /**< @brief Number of variable generators (size of Vargen array). */
  int num_bats;      /**< @brief Number of batteries (size of Bat array). */

  // Number of flags
  int num_vars;      /**< @brief Number of variable quantities. */
  int num_fixed;     /**< @brief Number of fixed quantities. */
  int num_bounded;   /**< @brief Number of bounded quantities. */
  int num_sparse;    /**< @brief Number of sparse control quantities. */

  // Base power
  REAL base_power; /**< @brief System base power (MVA). */

  // Properties
  REAL* bus_v_max;    /**< @brief Maximum bus voltage magnitude (p.u.). */
  REAL* bus_v_min;    /**< @brief Minimum bus volatge magnitude (p.u.). */
  REAL* bus_v_vio;    /**< @brief Maximum bus voltage magnitude limit violation (p.u.). */
  REAL* bus_P_mis;    /**< @brief Maximum bus active power mismatch (MW). */
  REAL* bus_Q_mis;    /**< @brief Maximum bus reactive power mismatch (MVAr). */

  REAL* gen_P_cost;   /**< @brief Total active power generation cost ($/hr). */
  REAL* gen_v_dev;    /**< @brief Maximum generator voltage setpoint deviation (p.u.). */
  REAL* gen_Q_vio;    /**< @brief Maximum generator reactive power limit violation (MVAr). */
  REAL* gen_P_vio;    /**< @brief Maximum generator active power limit violation (MW). */

  REAL* tran_v_vio;   /**< @brief Maximum transformer-controlled bus voltage magnitude band violation (p.u.). */
  REAL* tran_r_vio;   /**< @brief Maximum tap ratio limit violation of tap-changing transformer (unitless). */
  REAL* tran_p_vio;   /**< @brief Maximum phase shift limit violation of phase-shifting trasnformer (radians). */

  REAL* shunt_v_vio;  /**< @brief Maximum shunt-controlled bus voltage mangnitude band violation (p.u.). */
  REAL* shunt_b_vio;  /**< @brief Maximum susceptance limit volation of switched shunt device (p.u.). */

  REAL* load_P_util;  /**< @brief Total active power consumption utility ($/hr). */
  REAL* load_P_vio;   /**< @brief Maximum load active power limit violation (MW). */

  int* num_actions;   /**< @brief Number of control actions. */

  // Spatial correlation
  REAL vargen_corr_radius; /**< @brief Correlation radius for variable generators. **/
  REAL vargen_corr_value;  /**< @brief Correlation value for variable generators. **/

  // Utils
  char* bus_counted;  /**< @brief Flags for processing buses */
};

void NET_add_vargens(Net* net, Bus* bus_list, REAL power_capacity, REAL power_base, REAL power_std, REAL corr_radius, REAL corr_value) {
  
  // Local variables
  REAL total_load_P;
  REAL max_total_load_P;
  Vargen* vargen;
  int num;
  int i;
  int t;

  // Check
  if (!net)
    return;

  // Check
  if (power_capacity < 0 ||                 // percentage of max total load power
      power_base < 0 || power_base > 100 || // percentage of power capacity
      power_std < 0 ||                      // percentage of power capacity
      corr_radius < 0 ||                    // correlation radius
      corr_value < -1 || corr_value > 1) {  // correlation coefficient
    sprintf(net->error_string,"invalid arguments for adding variable generators");
    net->error_flag = TRUE;
    return;
  }

  // Clear
  VARGEN_hash_name_del(net->vargen_hash_name);
  VARGEN_array_del(net->vargen,net->num_vargens);
  net->vargen_hash_name = NULL;
  net->vargen = NULL;
  net->num_vargens = 0;

  // Save
  net->vargen_corr_radius = corr_radius;
  net->vargen_corr_value = corr_value;

  // Max total load power
  max_total_load_P = 0;
  for (t = 0; t < net->num_periods; t++) {
    total_load_P = 0;
    for (i = 0; i < net->num_loads; i++)
      total_load_P += LOAD_get_P(NET_get_load(net,i),t); // p.u.
    if (fabs(total_load_P) > max_total_load_P)
      max_total_load_P = fabs(total_load_P);
  }

  // Number
  num = BUS_list_len(bus_list);

  // Allocate
  NET_set_vargen_array(net,VARGEN_array_new(num,net->num_periods),num);

  // Set buses
  NET_set_vargen_buses(net,bus_list);

  // Set hash
  for (i = 0; i < net->num_vargens; i++) {
    vargen = NET_get_vargen(net,i);
    NET_vargen_hash_name_add(net,vargen);
  }

  // Set properties
  for (i = 0; i < net->num_vargens; i++) {
    vargen = NET_get_vargen(net,i);
    VARGEN_set_P_min(vargen,0.);
    VARGEN_set_P_max(vargen,(power_capacity/100.)*max_total_load_P/net->num_vargens);
    for (t = 0; t < net->num_periods; t++) {
      VARGEN_set_P_ava(vargen,(power_base/100.)*VARGEN_get_P_max(vargen),t);
      VARGEN_set_P(vargen,(power_base/100.)*VARGEN_get_P_max(vargen),t);
      VARGEN_set_P_std(vargen,(power_std/100.)*VARGEN_get_P_max(vargen),t);
    }
  }

  // Check hash
  if (VARGEN_hash_name_len(net->vargen_hash_name) != num) {
    sprintf(net->error_string,"unable to create vargen hash table");
    net->error_flag = TRUE;
    return;
  }
}

void NET_add_batteries(Net* net, Bus* bus_list, REAL power_capacity,  REAL energy_capacity, REAL eta_c, REAL eta_d) {
  
  // Local variables
  REAL total_load_P;
  REAL max_total_load_P;
  Bat* bat;
  int num;
  int i;
  int t;

  // Check
  if (!net)
    return;

  // Check
  if (power_capacity < 0 ||     // percentage of max total load power
      energy_capacity < 0 ||    // percentage of max total load energy during one interval
      eta_c <= 0 || eta_c > 1 || // charging efficiency in (0,1]
      eta_d <= 0 || eta_d > 1) { // discharging efficiency in (0,1]
    sprintf(net->error_string,"invalid arguments for adding batteries");
    net->error_flag = TRUE;
    return;
  }
  
  // Clear
  BAT_array_del(net->bat,net->num_bats);
  net->bat = NULL;
  net->num_bats = 0;

  // Max total load power
  max_total_load_P = 0;
  for (t = 0; t < net->num_periods; t++) {
    total_load_P = 0;
    for (i = 0; i < net->num_loads; i++)
      total_load_P += LOAD_get_P(NET_get_load(net,i),t); // p.u.
    if (fabs(total_load_P) > max_total_load_P)
      max_total_load_P = fabs(total_load_P);
  }

  // Number
  num = BUS_list_len(bus_list);

  // Allocate
  NET_set_bat_array(net,BAT_array_new(num,net->num_periods),num);

  // Set buses
  NET_set_bat_buses(net,bus_list);

  // Set properties
  for (i = 0; i < net->num_bats; i++) {
    bat = NET_get_bat(net,i);
    BAT_set_P_min(bat,-(power_capacity/100.)*max_total_load_P/net->num_bats);
    BAT_set_P_max(bat,(power_capacity/100.)*max_total_load_P/net->num_bats);
    BAT_set_eta_c(bat,eta_c);
    BAT_set_eta_d(bat,eta_d);
    BAT_set_E_max(bat,(energy_capacity/100.)*max_total_load_P/net->num_bats);
    BAT_set_E_init(bat,0.5*BAT_get_E_max(bat));
    BAT_set_E_final(bat,0.5*BAT_get_E_max(bat));
    for (t = 0; t < net->num_periods; t++) {
      BAT_set_P(bat,0,t);
      BAT_set_E(bat,BAT_get_E_init(bat),t);
    }
  }
}

void NET_adjust_generators(Net* net) {
  /** This function adjusts the powers of slack or regulator generators
   *  connected to the same bus or regulating the same bus voltage magnitude.
   *  The adjustment is done to obtain specific participations without affecting
   *  their total power. For active power, the participation is equal for every
   *  generator. For reactive power, the participaion is proportional to the generator
   *  reactive power resources.
   */

  // Local variables
  Bus* bus;
  Gen* gen;
  REAL num;
  REAL Ptot;
  REAL Qtot;
  REAL dQtot;
  REAL Q;
  REAL dQ;
  REAL Qmintot;
  REAL frac;
  int i;
  int t;

  // No net
  if (!net)
    return;

  for (i = 0; i < net->num_buses; i++) {

    bus = NET_get_bus(net,i);

    // Slack gens
    if (BUS_is_slack(bus)) {
      for (t = 0; t < net->num_periods; t++) {
	num = 0;
	Ptot = 0;
	for(gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	  Ptot += GEN_get_P(gen,t);
	  num += 1;
	}
	for(gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen))
	  GEN_set_P(gen,Ptot/num,t);
      }
    }

    // Regulating gens
    if (BUS_is_regulated_by_gen(bus)) {
      for (t = 0; t < net->num_periods; t++) {
	Qtot = 0;
	dQtot = 0;
	Qmintot = 0;
	for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	  Qtot += GEN_get_Q(gen,t);
	  dQtot += GEN_get_Q_max(gen)-GEN_get_Q_min(gen);
	  Qmintot += GEN_get_Q_min(gen);
	}
	gen = BUS_get_reg_gen(bus);
	dQ = GEN_get_Q_max(gen)-GEN_get_Q_min(gen);
	Q = GEN_get_Q_min(gen)+dQ*(Qtot-Qmintot)/dQtot;
	frac = (Q-GEN_get_Q_min(gen))/dQ;
	GEN_set_Q(gen,Q,t);
	for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen))
	  GEN_set_Q(gen,GEN_get_Q_min(gen)+frac*(GEN_get_Q_max(gen)-GEN_get_Q_min(gen)),t);
      }
    }
  }
}

void NET_bus_hash_number_add(Net* net, Bus* bus) {
  if (net)
    net->bus_hash_number = BUS_hash_number_add(net->bus_hash_number,bus);
}

Bus* NET_bus_hash_number_find(Net* net, int number) {
  if (net)
    return BUS_hash_number_find(net->bus_hash_number,number);
  else
    return NULL;
}

void NET_bus_hash_name_add(Net* net, Bus* bus) {
  if (net)
    net->bus_hash_name = BUS_hash_name_add(net->bus_hash_name,bus);
}

Bus* NET_bus_hash_name_find(Net* net, char* name) {
  if (net)
    return BUS_hash_name_find(net->bus_hash_name,name);
  else
    return NULL;
}

void NET_vargen_hash_name_add(Net* net, Vargen* vargen) {
  if (net)
    net->vargen_hash_name = VARGEN_hash_name_add(net->vargen_hash_name,vargen);
}

Vargen* NET_vargen_hash_name_find(Net* net, char* name) {
  if (net)
    return VARGEN_hash_name_find(net->vargen_hash_name,name);
  else
    return NULL;
}

BOOL NET_check(Net* net, BOOL verbose) {

  // Local variables
  BOOL bus_ok = TRUE;
  BOOL base_ok = TRUE;

  // Net
  if (!net) {
    fprintf(stderr,"NULL network\n");
  }

  // Base
  if (net->base_power > 0)
    base_ok = TRUE;
  else {
    base_ok = FALSE;
    if (verbose)
      fprintf(stderr,"non-positive base power\n");
  }

  // Buses
  if (!(net->bus)) {
    bus_ok = FALSE;
    if (verbose)
      fprintf(stderr,"NULL net bus array\n");
  }
  else
    bus_ok = BUS_array_check(net->bus,net->num_buses,verbose);

  // Branches

  // Gens

  // Loads

  // Shunt

  // Vargens

  // Bats

  // Overall
  return (base_ok & bus_ok);

}

void NET_clear_data(Net* net) {

  // No net
  if (!net)
    return;

  // Free hash tables
  BUS_hash_number_del(net->bus_hash_number);
  BUS_hash_name_del(net->bus_hash_name);
  VARGEN_hash_name_del(net->vargen_hash_name);

  // Free components
  BUS_array_del(net->bus,net->num_buses);
  BRANCH_array_del(net->branch,net->num_branches);
  GEN_array_del(net->gen,net->num_gens);
  SHUNT_array_del(net->shunt,net->num_shunts);
  LOAD_array_del(net->load,net->num_loads);
  VARGEN_array_del(net->vargen,net->num_vargens);
  BAT_array_del(net->bat,net->num_bats);

  // Free properties
  free(net->bus_v_max);
  free(net->bus_v_min);
  free(net->bus_v_vio);
  free(net->bus_P_mis);
  free(net->bus_Q_mis);
  free(net->gen_P_cost);
  free(net->gen_v_dev);
  free(net->gen_Q_vio);
  free(net->gen_P_vio);
  free(net->tran_v_vio);
  free(net->tran_r_vio);
  free(net->tran_p_vio);
  free(net->shunt_v_vio);
  free(net->shunt_b_vio);
  free(net->load_P_util);
  free(net->load_P_vio);
  free(net->num_actions);

  // Free utils
  free(net->bus_counted);

  // Re-initialize
  NET_init(net,net->num_periods);
}

void NET_clear_error(Net* net) {
  if (net) {
    net->error_flag = FALSE;
    strcpy(net->error_string,"");
  }
}

void NET_clear_flags(Net* net) {
  Branch* br;
  Gen* gen;
  Bus* bus;
  Shunt* shunt;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  int i;

  // Check
  if (!net)
    return;

  // Branches
  for (i = 0; i < net->num_branches; i++) {
    br = BRANCH_array_get(net->branch,i);
    BRANCH_clear_flags(br,FLAG_VARS);
    BRANCH_clear_flags(br,FLAG_FIXED);
    BRANCH_clear_flags(br,FLAG_BOUNDED);
    BRANCH_clear_flags(br,FLAG_SPARSE);
  }

  // Buses
  for (i = 0; i < net->num_buses; i++) {
    bus = BUS_array_get(net->bus,i);
    BUS_clear_flags(bus,FLAG_VARS);
    BUS_clear_flags(bus,FLAG_FIXED);
    BUS_clear_flags(bus,FLAG_BOUNDED);
    BUS_clear_flags(bus,FLAG_SPARSE);
  }

  // Gens
  for (i = 0; i < net->num_gens; i++) {
    gen = GEN_array_get(net->gen,i);
    GEN_clear_flags(gen,FLAG_VARS);
    GEN_clear_flags(gen,FLAG_FIXED);
    GEN_clear_flags(gen,FLAG_BOUNDED);
    GEN_clear_flags(gen,FLAG_SPARSE);
  }

  // Shunts
  for (i = 0; i < net->num_shunts; i++) {
    shunt = SHUNT_array_get(net->shunt,i);
    SHUNT_clear_flags(shunt,FLAG_VARS);
    SHUNT_clear_flags(shunt,FLAG_FIXED);
    SHUNT_clear_flags(shunt,FLAG_BOUNDED);
    SHUNT_clear_flags(shunt,FLAG_SPARSE);
  }

  // Loads
  for (i = 0; i < net->num_loads; i++) {
    load = LOAD_array_get(net->load,i);
    LOAD_clear_flags(load,FLAG_VARS);
    LOAD_clear_flags(load,FLAG_FIXED);
    LOAD_clear_flags(load,FLAG_BOUNDED);
    LOAD_clear_flags(load,FLAG_SPARSE);
  }

  // Vargens
  for (i = 0; i < net->num_vargens; i++) {
    vargen = VARGEN_array_get(net->vargen,i);
    VARGEN_clear_flags(vargen,FLAG_VARS);
    VARGEN_clear_flags(vargen,FLAG_FIXED);
    VARGEN_clear_flags(vargen,FLAG_BOUNDED);
    VARGEN_clear_flags(vargen,FLAG_SPARSE);
  }

  // Batteries
  for (i = 0; i < net->num_bats; i++) {
    bat = BAT_array_get(net->bat,i);
    BAT_clear_flags(bat,FLAG_VARS);
    BAT_clear_flags(bat,FLAG_FIXED);
    BAT_clear_flags(bat,FLAG_BOUNDED);
    BAT_clear_flags(bat,FLAG_SPARSE);
  }

  // Clear counters
  net->num_vars = 0;
  net->num_fixed = 0;
  net->num_bounded = 0;
  net->num_sparse = 0;
}

void NET_clear_outages(Net* net) {

  // Local vars
  int i;

  // No net
  if (!net)
    return;

  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_set_outage(NET_get_gen(net,i),FALSE);

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_set_outage(NET_get_branch(net,i),FALSE);
}

void NET_clear_properties(Net* net) {

  // Local variables
  int i;
  int t;

  // No net
  if (!net)
    return;

  // Time loop
  for (t = 0; t < net->num_periods; t++) {

    // Bus
    net->bus_v_max[t] = 0;
    net->bus_v_min[t] = 0;
    net->bus_v_vio[t] = 0;
    net->bus_P_mis[t] = 0;
    net->bus_Q_mis[t] = 0;

    // Gen
    net->gen_P_cost[t] = 0;
    net->gen_v_dev[t] = 0;
    net->gen_Q_vio[t] = 0;
    net->gen_P_vio[t] = 0;

    // Branch
    net->tran_v_vio[t] = 0;
    net->tran_r_vio[t] = 0;
    net->tran_p_vio[t] = 0;

    // Shunt
    net->shunt_v_vio[t] = 0;
    net->shunt_b_vio[t] = 0;

    // Load
    net->load_P_util[t] = 0;
    net->load_P_vio[t] = 0;

    // Battery

    // Actions
    net->num_actions[t] = 0;
  }

  // Counters
  for (i = 0; i < net->num_buses; i++)
    BUS_clear_mismatches(BUS_array_get(net->bus,i));

  // Bus counted
  ARRAY_clear(net->bus_counted,char,net->num_buses*net->num_periods);
}

void NET_clear_sensitivities(Net* net) {

  // Local variables
  int i;

  // No net
  if (!net)
    return;

  // Buses
  for (i = 0; i < net->num_buses; i++)
    BUS_clear_sensitivities(BUS_array_get(net->bus,i));

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_clear_sensitivities(BRANCH_array_get(net->branch,i));

  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_clear_sensitivities(GEN_array_get(net->gen,i));

  // Loads
  for (i = 0; i < net->num_loads; i++)
    LOAD_clear_sensitivities(LOAD_array_get(net->load,i));

  // Vargens

  // Shunts

  // Batteries
}

Bus* NET_create_sorted_bus_list(Net* net, int sort_by, int t) {

  // Local variables
  Bus* bus_list = NULL;
  int i;

  if (!net || t < 0 || t >= net->num_periods)
    return bus_list;

  for (i = 0; i < net->num_buses; i++)
    bus_list = BUS_list_add_sorting(bus_list,BUS_array_get(net->bus,i),sort_by,t);
  return bus_list;
}

int NET_get_bus_neighbors(Net* net, Bus* bus, int spread, int* neighbors, char* queued) {
  /** Returns number of neighbors including itself that are at most "spread"
   *  branches away.
   */

  // Local variables
  Bus* bus1;
  Bus* bus2;
  Branch* br;
  int neighbors_total;
  int neighbors_curr;
  int num_new;
  int i;

  // Check
  if (!neighbors || !queued)
    return -1;

  // Add self to be processed
  neighbors_total = 1;
  neighbors[0] = BUS_get_index(bus);
  queued[BUS_get_index(bus)] = TRUE;

  // Neighbors
  neighbors_curr = 0;
  for (i = 0; i < spread; i++) {
    num_new = 0;
    while (neighbors_curr < neighbors_total) {
      bus1 = NET_get_bus(net,neighbors[neighbors_curr]);
      for (br = BUS_get_branch_k(bus1); br != NULL; br = BRANCH_get_next_k(br)) {
	if (bus1 != BRANCH_get_bus_k(br)) {
	  sprintf(net->error_string,"unable to construct covariance matrix");
	  net->error_flag = TRUE;
	}
	bus2 = BRANCH_get_bus_m(br);
	if (!queued[BUS_get_index(bus2)]) {
	  neighbors[neighbors_total+num_new] = BUS_get_index(bus2);
	  queued[BUS_get_index(bus2)] = TRUE;
	  num_new++;
	}
      }
      for (br = BUS_get_branch_m(bus1); br != NULL; br = BRANCH_get_next_m(br)) {
	if (bus1 != BRANCH_get_bus_m(br)) {
	  sprintf(net->error_string,"unable to construct covariance matrix");
	  net->error_flag = TRUE;
	}
	bus2 = BRANCH_get_bus_k(br);
	if (!queued[BUS_get_index(bus2)]) {
	  neighbors[neighbors_total+num_new] = BUS_get_index(bus2);
	  queued[BUS_get_index(bus2)] = TRUE;
	  num_new++;
	}
      }
      neighbors_curr++;
    }
    neighbors_total += num_new;
    if (num_new == 0)
      break;
  }
  return neighbors_total;
}

Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr) {
  /** This function constructs a "spatial" covariance matrix for the active powers of
   *  variable generators. The matrix is constructed such that the correlation
   *  coefficients of the (variable) active powers of vargens that are less than
   *  "spread" branches away is equal to "corr". Only the lower triangular part
   *  of the covaraicen matrix is stored. The resulting matrix should be checked
   *  to make sure it is a valid covariance matrix.
   */

  // Local variables
  Mat* sigma;
  Bus* bus_main;
  Bus* bus;
  Vargen* vgen_main;
  Vargen* vg;
  char* queued;
  int* neighbors;
  int nnz_counter;
  int num_neighbors;
  int i;
  int j;
  int t;

  // Check
  if (!net)
    return NULL;

  // Allocate arrays
  ARRAY_alloc(queued,char,net->num_buses);
  ARRAY_alloc(neighbors,int,net->num_buses);

  // Count nnz
  //**********
  nnz_counter = 0;
  for (i = 0; i < net->num_vargens; i++) {

    // Clear arrays
    for (j = 0; j < net->num_buses; j++) {
      neighbors[j] = 0;
      queued[j] = FALSE;
    }

    // Main
    vgen_main = NET_get_vargen(net,i);
    bus_main = VARGEN_get_bus(vgen_main);

    // Check variable
    if (!VARGEN_has_flags(vgen_main,FLAG_VARS,VARGEN_VAR_P))
      continue;

    // Neighbors
    num_neighbors = NET_get_bus_neighbors(net,bus_main,spread,neighbors,queued);
    if (num_neighbors < 0) {
      sprintf(net->error_string,"unable to construct covariance matrix");
      net->error_flag = TRUE;
    }

    // Diagonals
    nnz_counter += net->num_periods;

    // Off diagonals
    for (j = 0; j < num_neighbors; j++) {
      bus = NET_get_bus(net,neighbors[j]);
      for (vg = BUS_get_vargen(bus); vg != NULL; vg = VARGEN_get_next(vg)) {
	for (t = 0; t < net->num_periods; t++) {
	  if (VARGEN_has_flags(vg,FLAG_VARS,VARGEN_VAR_P) &&
	      VARGEN_get_index_P(vgen_main,t) > VARGEN_get_index_P(vg,t)) {
	    nnz_counter++;
	  }
	}
      }
    }
  }

  // Allocate
  //*********
  sigma = MAT_new(net->num_vars,
		  net->num_vars,
		  nnz_counter);

  // Fill
  //*****
  nnz_counter = 0;
  for (i = 0; i < net->num_vargens; i++) {

    // Clear arrays
    for (j = 0; j < net->num_buses; j++) {
      neighbors[j] = 0;
      queued[j] = FALSE;
    }

    // Main
    vgen_main = NET_get_vargen(net,i);
    bus_main = VARGEN_get_bus(vgen_main);

    // Check variable
    if (!VARGEN_has_flags(vgen_main,FLAG_VARS,VARGEN_VAR_P))
      continue;

    // Neighbors
    num_neighbors = NET_get_bus_neighbors(net,bus_main,spread,neighbors,queued);
    if (num_neighbors < 0) {
      sprintf(net->error_string,"unable to construct covariance matrix");
      net->error_flag = TRUE;
    }

    // Diagonal
    for (t = 0; t < net->num_periods; t++) {
      MAT_set_i(sigma,nnz_counter,VARGEN_get_index_P(vgen_main,t));
      MAT_set_j(sigma,nnz_counter,VARGEN_get_index_P(vgen_main,t));
      MAT_set_d(sigma,nnz_counter,pow(VARGEN_get_P_std(vgen_main,t),2.));
      nnz_counter++;
    }

    // Off diagonals
    for (j = 0; j < num_neighbors; j++) {
      bus = NET_get_bus(net,neighbors[j]);
      for (vg = BUS_get_vargen(bus); vg != NULL; vg = VARGEN_get_next(vg)) {
	for (t = 0; t < net->num_periods; t++) {
	  if (VARGEN_has_flags(vg,FLAG_VARS,VARGEN_VAR_P) &&
	      VARGEN_get_index_P(vgen_main,t) > VARGEN_get_index_P(vg,t)) {
	    MAT_set_i(sigma,nnz_counter,VARGEN_get_index_P(vgen_main,t));
	    MAT_set_j(sigma,nnz_counter,VARGEN_get_index_P(vg,t));
	    MAT_set_d(sigma,nnz_counter,VARGEN_get_P_std(vgen_main,t)*VARGEN_get_P_std(vg,t)*corr);
	    nnz_counter++;
	  }
	}
      }
    }
  }

  // Check
  if (nnz_counter != MAT_get_nnz(sigma)) {
    sprintf(net->error_string,"unable to construct covariance matrix");
    net->error_flag = TRUE;
  }

  // Clean up
  free(queued);
  free(neighbors);

  // Return
  return sigma;
}

void NET_del(Net* net) {
  if (net) {
    NET_clear_data(net);
    free(net);
  }
}

void NET_init(Net* net, int num_periods) {

  // Local vars
  int T;

  // No net
  if (!net)
    return;

  // Number of periods
  T = num_periods;
  net->num_periods = num_periods;

  // Error
  net->error_flag = FALSE;
  strcpy(net->error_string,"");

  // Output
  strcpy(net->output_string,"");

  // Components
  net->bus = NULL;
  net->branch = NULL;
  net->gen = NULL;
  net->load = NULL;
  net->shunt = NULL;
  net->vargen = NULL;
  net->bat = NULL;

  // Hash tables
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;
  net->vargen_hash_name = NULL;

  // Number of components
  net->num_buses = 0;
  net->num_branches = 0;
  net->num_gens = 0;
  net->num_loads = 0;
  net->num_shunts = 0;
  net->num_vargens = 0;
  net->num_bats = 0;

  // Number flags
  net->num_vars = 0;
  net->num_fixed = 0;
  net->num_bounded = 0;
  net->num_sparse = 0;

  // Base power
  net->base_power = NET_BASE_POWER;

  // Spatial correlation
  net->vargen_corr_radius = 1;
  net->vargen_corr_value = 0;

  // Properties
  ARRAY_zalloc(net->bus_v_max,REAL,T);
  ARRAY_zalloc(net->bus_v_min,REAL,T);
  ARRAY_zalloc(net->bus_v_vio,REAL,T);
  ARRAY_zalloc(net->bus_P_mis,REAL,T);
  ARRAY_zalloc(net->bus_Q_mis,REAL,T);

  ARRAY_zalloc(net->gen_P_cost,REAL,T);
  ARRAY_zalloc(net->gen_v_dev,REAL,T);
  ARRAY_zalloc(net->gen_Q_vio,REAL,T);
  ARRAY_zalloc(net->gen_P_vio,REAL,T);

  ARRAY_zalloc(net->tran_v_vio,REAL,T);
  ARRAY_zalloc(net->tran_r_vio,REAL,T);
  ARRAY_zalloc(net->tran_p_vio,REAL,T);

  ARRAY_zalloc(net->shunt_v_vio,REAL,T);
  ARRAY_zalloc(net->shunt_b_vio,REAL,T);

  ARRAY_zalloc(net->load_P_util,REAL,T);
  ARRAY_zalloc(net->load_P_vio,REAL,T);

  ARRAY_zalloc(net->num_actions,int,T);

  // Utils
  net->bus_counted = NULL;
}

REAL NET_get_base_power(Net* net) {
  if (net)
    return net->base_power;
  else
    return NET_BASE_POWER;
}

Branch* NET_get_branch(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_branches)
    return NULL;
  else
    return BRANCH_array_get(net->branch,index);
}

Bus* NET_get_bus(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_buses)
    return NULL;
  else
    return BUS_array_get(net->bus,index);
}

Bus* NET_get_bus_hash_number(Net* net) {
  if (!net)
    return NULL;
  else
    return net->bus_hash_number;
}

Bus* NET_get_bus_hash_name(Net* net) {
  if (!net)
    return NULL;
  else
    return net->bus_hash_name;
}

Vargen* NET_get_vargen_hash_name(Net* net) {
  if (!net)
    return NULL;
  else
    return net->vargen_hash_name;
}

char* NET_get_error_string(Net* net) {
  if (!net)
    return NULL;
  else
    return net->error_string;
}

Gen* NET_get_gen(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_gens)
    return NULL;
  else
    return GEN_array_get(net->gen,index);
}

Load* NET_get_load(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_loads)
    return NULL;
  else
    return LOAD_array_get(net->load,index);
}

Shunt* NET_get_shunt(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_shunts)
    return NULL;
  else
    return SHUNT_array_get(net->shunt,index);
}

Vargen* NET_get_vargen(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_vargens)
    return NULL;
  else
    return VARGEN_array_get(net->vargen,index);
}

Bat* NET_get_bat(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_bats)
    return NULL;
  else
    return BAT_array_get(net->bat,index);
}

Bus* NET_get_gen_buses(Net* net) {

  Bus* bus_list = NULL;
  Bus* bus;
  int i;

  if (!net)
    return bus_list;

  for (i = 0; i < net->num_buses; i++) {
    bus = NET_get_bus(net,i);
    if (BUS_get_gen(bus))
      bus_list = BUS_list_add(bus_list,bus);
  }
  return bus_list;
}

Bus* NET_get_load_buses(Net* net) {

  Bus* bus_list = NULL;
  Bus* bus;
  int i;

  if (!net)
    return bus_list;

  for (i = 0; i < net->num_buses; i++) {
    bus = NET_get_bus(net,i);
    if (BUS_get_load(bus))
      bus_list = BUS_list_add(bus_list,bus);
  }
  return bus_list;
}

int NET_get_num_periods(Net* net) {
  if (net)
    return net->num_periods;
  else
    return 0;
}

int NET_get_num_buses(Net* net) {
  if (net)
    return net->num_buses;
  else
    return 0;
}

int NET_get_num_slack_buses(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_slack(BUS_array_get(net->bus,i)))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_gen(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_regulated_by_gen(BUS_array_get(net->bus,i)))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_tran(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_regulated_by_tran(BUS_array_get(net->bus,i)))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_tran_only(Net* net) {
  int i;
  int n = 0;
  Bus* bus;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    bus = BUS_array_get(net->bus,i);
    if (BUS_is_regulated_by_tran(bus) &&
	!BUS_is_regulated_by_gen(bus) &&
	!BUS_is_regulated_by_shunt(bus))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_shunt(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_regulated_by_shunt(BUS_array_get(net->bus,i)))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_shunt_only(Net* net) {
  int i;
  int n = 0;
  Bus* bus;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    bus = BUS_array_get(net->bus,i);
    if (BUS_is_regulated_by_shunt(bus) &&
	!BUS_is_regulated_by_gen(bus) &&
	!BUS_is_regulated_by_tran(bus))
      n++;
  }
  return n;
}

int NET_get_num_branches(Net* net) {
  if (net)
    return net->num_branches;
  else
    return 0;
}

int NET_get_num_branches_not_on_outage(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_branches; i++) {
    if (!BRANCH_is_on_outage(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_fixed_trans(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_fixed_tran(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_lines(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_line(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_phase_shifters(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_phase_shifter(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_tap_changers(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_tap_changer(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_tap_changers_v(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_tap_changer_v(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_tap_changers_Q(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_branches; i++) {
    if (BRANCH_is_tap_changer_Q(BRANCH_array_get(net->branch,i)))
      n++;
  }
  return n;
}

int NET_get_num_gens(Net* net) {
  if (net)
    return net->num_gens;
  else
    return 0;
}

int NET_get_num_gens_not_on_outage(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_gens; i++) {
    if (!GEN_is_on_outage(GEN_array_get(net->gen,i)))
      n++;
  }
  return n;
}

int NET_get_num_reg_gens(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_gens; i++) {
    if (GEN_is_regulator(GEN_array_get(net->gen,i)))
      n++;
  }
  return n;
}

int NET_get_num_slack_gens(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_gens; i++) {
    if (GEN_is_slack(GEN_array_get(net->gen,i)))
      n++;
  }
  return n;
}

int NET_get_num_P_adjust_gens(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_gens; i++) {
    if (GEN_is_P_adjustable(GEN_array_get(net->gen,i)))
      n++;
  }
  return n;
}

int NET_get_num_loads(Net* net) {
  if (net)
    return net->num_loads;
  else
    return 0;
}

int NET_get_num_P_adjust_loads(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_loads; i++) {
    if (LOAD_is_P_adjustable(LOAD_array_get(net->load,i)))
      n++;
  }
  return n;
}

int NET_get_num_shunts(Net* net) {
  if (net)
    return net->num_shunts;
  else
    return 0;
}

int NET_get_num_fixed_shunts(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_shunts; i++) {
    if (SHUNT_is_fixed(SHUNT_array_get(net->shunt,i)))
      n++;
  }
  return n;
}

int NET_get_num_switched_shunts(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_shunts; i++) {
    if (SHUNT_is_switched(SHUNT_array_get(net->shunt,i)))
      n++;
  }
  return n;
}

int NET_get_num_vargens(Net* net) {
  if (net)
    return net->num_vargens;
  else
    return 0;
}

int NET_get_num_bats(Net* net) {
  if (net)
    return net->num_bats;
  else
    return 0;
}

int NET_get_num_bounded(Net* net) {
  if (net)
    return net->num_bounded;
  else
    return 0;
}

int NET_get_num_sparse(Net* net) {
  if (net)
    return net->num_sparse;
  else
    return 0;
}

int NET_get_num_fixed(Net* net) {
  if (net)
    return net->num_fixed;
  else
    return 0;
}

int NET_get_num_vars(Net* net) {
  if (net)
    return net->num_vars;
  else
    return 0;
}

REAL NET_get_total_gen_P(Net* net, int t) {
  int i;
  REAL P = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_gens; i++) {
    P += GEN_get_P(GEN_array_get(net->gen,i),t); // p.u.
  }
  return P*net->base_power; // MW
}

REAL NET_get_total_gen_Q(Net* net, int t) {
  int i;
  REAL Q = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_gens; i++) {
    Q += GEN_get_Q(GEN_array_get(net->gen,i),t); // p.u.
  }
  return Q*net->base_power; // MVAr
}

REAL NET_get_total_load_P(Net* net, int t) {
  int i;
  REAL P = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_loads; i++) {
    P += LOAD_get_P(LOAD_array_get(net->load,i),t); // p.u.
  }
  return P*net->base_power; // MW
}

REAL NET_get_total_load_Q(Net* net, int t) {
  int i;
  REAL Q = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_loads; i++) {
    Q += LOAD_get_Q(LOAD_array_get(net->load,i),t); // p.u.
  }
  return Q*net->base_power; // MVAr
}

Vec* NET_get_var_values(Net* net, int code) {

  // Local variables
  int i;
  Vec* values;

  if (!net)
    return NULL;

  // Vector of var values
  values = VEC_new(net->num_vars);

  // Buses
  for (i = 0; i < net->num_buses; i++)
    BUS_get_var_values(BUS_array_get(net->bus,i),values,code);

  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_get_var_values(GEN_array_get(net->gen,i),values,code);

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_get_var_values(BRANCH_array_get(net->branch,i),values,code);

  // Shunts
  for (i = 0; i < net->num_shunts; i++)
    SHUNT_get_var_values(SHUNT_array_get(net->shunt,i),values,code);

  // Loads
  for (i = 0; i < net->num_loads; i++)
    LOAD_get_var_values(LOAD_array_get(net->load,i),values,code);

  // Variable generators
  for (i = 0; i < net->num_vargens; i++)
    VARGEN_get_var_values(VARGEN_array_get(net->vargen,i),values,code);

  // Batteries
  for (i = 0; i < net->num_bats; i++)
    BAT_get_var_values(BAT_array_get(net->bat,i),values,code);

  // Return
  return values;
}

Mat* NET_get_var_projection(Net* net, char obj_type, unsigned char var, int t_start, int t_end) {

  // Local variables
  int num_subvars;
  Vec* indices;
  Mat* proj;
  int i;
  int j;

  // Check
  if (!net)
    return NULL;

  // Check
  if (t_start < 0)
    t_start = 0;
  if (t_end > net->num_periods-1)
    t_end = net->num_periods-1;

  // Check
  if ((obj_type == OBJ_ALL) && (var != 0xFF)) {
    sprintf(net->error_string,"component-specific flag cannot be used on all components");
    net->error_flag = TRUE;
    return NULL;
  }

  // Count
  num_subvars = 0;
  if ((obj_type == OBJ_BUS) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_buses; i++)
      num_subvars += BUS_get_num_vars(NET_get_bus(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_GEN) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_gens; i++)
      num_subvars += GEN_get_num_vars(NET_get_gen(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_LOAD) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_loads; i++)
      num_subvars += LOAD_get_num_vars(NET_get_load(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_BRANCH) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_branches; i++)
      num_subvars += BRANCH_get_num_vars(NET_get_branch(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_SHUNT) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_shunts; i++)
      num_subvars += SHUNT_get_num_vars(NET_get_shunt(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_VARGEN) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_vargens; i++)
      num_subvars += VARGEN_get_num_vars(NET_get_vargen(net,i),var,t_start,t_end);
  }
  if ((obj_type == OBJ_BAT) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_bats; i++)
      num_subvars += BAT_get_num_vars(NET_get_bat(net,i),var,t_start,t_end);
  }

  // Allocate
  proj = MAT_new(num_subvars,
		 net->num_vars,
		 num_subvars);

  // Fill
  num_subvars = 0;
  if ((obj_type == OBJ_BUS) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_buses; i++) {
      indices = BUS_get_var_indices(NET_get_bus(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_GEN) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_gens; i++) {
      indices = GEN_get_var_indices(NET_get_gen(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_LOAD) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_loads; i++) {
      indices = LOAD_get_var_indices(NET_get_load(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_BRANCH) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_branches; i++) {
      indices = BRANCH_get_var_indices(NET_get_branch(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_SHUNT) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_shunts; i++) {
      indices = SHUNT_get_var_indices(NET_get_shunt(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_VARGEN) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_vargens; i++) {
      indices = VARGEN_get_var_indices(NET_get_vargen(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }
  if ((obj_type == OBJ_BAT) || (obj_type == OBJ_ALL)) {
    for (i = 0; i < net->num_bats; i++) {
      indices = BAT_get_var_indices(NET_get_bat(net,i),var,t_start,t_end);
      for (j = 0; j < VEC_get_size(indices); j++) {
	MAT_set_i(proj,num_subvars,num_subvars);
	MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
	MAT_set_d(proj,num_subvars,1.);
	num_subvars++;
      }
      VEC_del(indices);
    }
  }

  // Return
  return proj;
}

REAL NET_get_bus_v_max(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->bus_v_max[t];
  else
    return 0;
}

REAL NET_get_bus_v_min(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->bus_v_min[t];
  else
    return 0;
}

REAL NET_get_bus_v_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->bus_v_vio[t];
  else
    return 0;
}

REAL NET_get_bus_P_mis(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->bus_P_mis[t];
  else
    return 0;
}

REAL NET_get_bus_Q_mis(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->bus_Q_mis[t];
  else
    return 0;
}

REAL NET_get_gen_P_cost(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->gen_P_cost[t];
  else
    return 0;
}

REAL NET_get_gen_v_dev(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->gen_v_dev[t];
  else
    return 0;
}

REAL NET_get_gen_Q_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->gen_Q_vio[t];
  else
    return 0;
}

REAL NET_get_gen_P_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->gen_P_vio[t];
  else
    return 0;
}

REAL NET_get_tran_v_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->tran_v_vio[t];
  else
    return 0;
}

REAL NET_get_tran_r_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->tran_r_vio[t];
  else
    return 0;
}

REAL NET_get_tran_p_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->tran_p_vio[t];
  else
    return 0;
}

REAL NET_get_shunt_v_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->shunt_v_vio[t];
  else
    return 0;
}

REAL NET_get_shunt_b_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->shunt_b_vio[t];
  else
    return 0;
}

REAL NET_get_load_P_util(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->load_P_util[t];
  else
    return 0;
}

REAL NET_get_load_P_vio(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->load_P_vio[t];
  else
    return 0;
}

int NET_get_num_actions(Net* net, int t) {
  if (net && t >= 0 && t < net->num_periods)
    return net->num_actions[t];
  else
    return 0;
}

REAL NET_get_vargen_corr_radius(Net* net) {
  if (net)
    return net->vargen_corr_radius;
  else
    return 0;
}

REAL NET_get_vargen_corr_value(Net* net) {
  if (net)
    return net->vargen_corr_value;
  else
    return 0;
}

char* NET_get_json_string(Net* net) {

   // Local variables
  char temp[NET_BUFFER_SIZE];
  char* output;
  char* element_output;
  int max_size;
  int i;

  // No network
  if (!net)
    return NULL;

  // Max size
  max_size = (2*NET_BUFFER_SIZE +
	      BUS_BUFFER_SIZE*BUS_NUM_JSON_FIELDS*net->num_buses +
	      BRANCH_BUFFER_SIZE*BRANCH_NUM_JSON_FIELDS*net->num_branches +
	      GEN_BUFFER_SIZE*GEN_NUM_JSON_FIELDS*net->num_gens +
	      LOAD_BUFFER_SIZE*LOAD_NUM_JSON_FIELDS*net->num_loads +
	      SHUNT_BUFFER_SIZE*SHUNT_NUM_JSON_FIELDS*net->num_shunts +
	      VARGEN_BUFFER_SIZE*VARGEN_NUM_JSON_FIELDS*net->num_vargens +
	      BAT_BUFFER_SIZE*BAT_NUM_JSON_FIELDS*net->num_bats)*net->num_periods;

  // Alloc
  output = (char*)malloc(sizeof(char)*max_size);
  
  // Start
  strcpy(output,"{ ");
  
  // Num periods
  sprintf(temp,"\"num_periods\" : %d", net->num_periods);
  strcat(output,temp);
  strcat(output,", ");

  // Base power
  sprintf(temp,"\"base_power\" : %.10e", net->base_power);
  strcat(output,temp);
  strcat(output,", ");

  // Buses
  element_output = (char*)malloc(sizeof(char)*BUS_BUFFER_SIZE*BUS_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"buses\" : [ ");
  for (i = 0; i < net->num_buses; i++) {
    BUS_get_json_string(BUS_array_get(net->bus,i),element_output);
    strcat(output,element_output);
    if (i < net->num_buses-1)
      strcat(output,", ");
  }
  strcat(output," ]");
  free(element_output);

  /*
  // Branches
  element_output = (char*)malloc(sizeof(char)*BRANCH_BUFFER_SIZE*BRANCH_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"branches\" : [ ");
  for (i = 0; i < net->num_branches; i++) {
    BRANCH_get_json_string(BRANCH_array_get(net->branch,i),element_output);
    strcat(output,element_output);
    if (i < net->num_branches-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  free(element_output);

  // Generators
  element_output = (char*)malloc(sizeof(char)*GEN_BUFFER_SIZE*GEN_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"generators\" : [ ");
  for (i = 0; i < net->num_gens; i++) {
    GEN_get_json_string(GEN_array_get(net->gen,i),element_output);
    strcat(output,element_output);
    if (i < net->num_gens-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  free(element_output);

  // Loads
  element_output = (char*)malloc(sizeof(char)*LOAD_BUFFER_SIZE*LOAD_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"loads\" : [ ");
  for (i = 0; i < net->num_loads; i++) {
    LOAD_get_json_string(LOAD_array_get(net->load,i),element_output);
    strcat(output,element_output);
    if (i < net->num_loads-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  free(element_output);

  // Shunts
  element_output = (char*)malloc(sizeof(char)*SHUNT_BUFFER_SIZE*SHUNT_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"shunts\" : [ ");
  for (i = 0; i < net->num_shunts; i++) {
    SHUNT_get_json_string(SHUNT_array_get(net->shunt,i),element_output);
    strcat(output,element_output);
    if (i < net->num_shunts-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  free(element_output);

  // Var generators
  element_output = (char*)malloc(sizeof(char)*VARGEN_BUFFER_SIZE*VARGEN_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"var_generators\" : [ ");
  for (i = 0; i < net->num_vargens; i++) {
    VARGEN_get_json_string(VARGEN_array_get(net->vargen,i),element_output);
    strcat(output,element_output);
    if (i < net->num_vargens-1)
      strcat(output,", ");
  }
  strcat(output," ], ");
  free(element_output);

  // Batteries
  element_output = (char*)malloc(sizeof(char)*BAT_BUFFER_SIZE*BAT_NUM_JSON_FIELDS*net->num_periods);
  strcat(output,"\"batteries\" : [ ");
  for (i = 0; i < net->num_bats; i++) {
    BAT_get_json_string(BAT_array_get(net->bat,i),element_output);
    strcat(output,element_output);
    if (i < net->num_bats-1)
      strcat(output,", ");
  }
  strcat(output," ]");
  free(element_output);
  */

  // End
  strcat(output," }");
  
  // Resize
  output = (char*)realloc(output,sizeof(char)*(strlen(output)+1)); // +1 important!

  // Return
  return output;
}

BOOL NET_has_error(Net* net) {
  if (net)
    return net->error_flag;
  else
    return FALSE;
}

Net* NET_new(int num_periods) {
  if (num_periods > 0) {
    Net* net = (Net*)malloc(sizeof(Net));
    NET_init(net,num_periods);
    return net;
  }
  else
    return NULL;
}

void NET_set_base_power(Net* net, REAL base_power) {
  if (net)
    net->base_power = base_power;
}

void NET_set_branch_array(Net* net, Branch* branch, int num) {
  if (net) {
    net->branch = branch;
    net->num_branches = num;
  }
}

void NET_set_load_array(Net* net, Load* load, int num) {
  if (net) {
    net->load = load;
    net->num_loads = num;
  }
}

void NET_set_shunt_array(Net* net, Shunt* shunt, int num) {
  if (net) {
    net->shunt = shunt;
    net->num_shunts = num;
  }
}

void NET_set_bus_array(Net* net, Bus* bus, int num) {
  if (net) {
    net->bus = bus;
    net->num_buses = num;
    ARRAY_zalloc(net->bus_counted,char,net->num_buses*net->num_periods);
  }
}

void NET_set_gen_array(Net* net, Gen* gen, int num) {
  if (net) {
    net->gen = gen;
    net->num_gens = num;
  }
}

void NET_set_vargen_array(Net* net, Vargen* gen, int num) {

  // Local variables
  int i;

  if (net) {

    // Clear vargens connections to buses
    for (i = 0; i < net->num_buses; i++)
      BUS_clear_vargen(NET_get_bus(net,i));

    // Clear hash
    VARGEN_hash_name_del(net->vargen_hash_name);
    net->vargen_hash_name = NULL;

    // Clear array
    VARGEN_array_del(net->vargen,net->num_vargens);
    net->vargen = NULL;
    net->num_vargens = 0;

    // Check hash length
    if (VARGEN_hash_name_len(net->vargen_hash_name) != 0) {
      sprintf(net->error_string,"unable to clear vargen hash table");
      net->error_flag = TRUE;
      return;
    }

    // Set
    net->vargen = gen;         // array
    net->num_vargens = num;    // number
  }
}

void NET_set_bat_array(Net* net, Bat* bat, int num) {
  if (net) {
    net->bat = bat;
    net->num_bats = num;
  }
}

void NET_set_vargen_buses(Net* net, Bus* bus_list) {

  // Local vars
  int i;
  Bus* bus;
  Vargen* gen;

  // No net
  if (!net)
    return;

  // Clear connections
  for (i = 0; i < net->num_buses; i++)
    BUS_clear_vargen(NET_get_bus(net,i));

  i = 0;
  bus = bus_list;
  while (i < net->num_vargens && bus) {
    gen = VARGEN_array_get(net->vargen,i);
    VARGEN_set_bus(gen,bus);
    BUS_add_vargen(bus,gen);
    bus = BUS_get_next(bus);
    i++;
  }
}

void NET_set_bat_buses(Net* net, Bus* bus_list) {

  // Local vars
  int i;
  Bus* bus;
  Bat* bat;

  // No net
  if (!net)
    return;

  // Clear connections
  for (i = 0; i < net->num_buses; i++)
    BUS_clear_bat(NET_get_bus(net,i));

  i = 0;
  bus = bus_list;
  while (i < net->num_bats && bus) {
    bat = BAT_array_get(net->bat,i);
    BAT_set_bus(bat,bus);
    BUS_add_bat(bus,bat);
    bus = BUS_get_next(bus);
    i++;
  }
}

void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, unsigned char val_mask) {

  // Local variables
  int i;
  int num;
  void* obj;
  void* array;
  void* (*get_element)(void* array, int index);
  int (*set_flags)(void*,char,unsigned char,int);
  BOOL (*has_properties)(void*,char);

  // Check
  if (!net)
    return;

  // Set pointers
  switch (obj_type) {
  case OBJ_BUS:
    num = net->num_buses;
    array = net->bus;
    get_element = &BUS_array_get;
    set_flags = &BUS_set_flags;
    has_properties = &BUS_has_properties;
    break;
  case OBJ_GEN:
    num = net->num_gens;
    array = net->gen;
    get_element = &GEN_array_get;
    set_flags = &GEN_set_flags;
    has_properties = &GEN_has_properties;
    break;
  case OBJ_LOAD:
    num = net->num_loads;
    array = net->load;
    get_element = &LOAD_array_get;
    set_flags = &LOAD_set_flags;
    has_properties = &LOAD_has_properties;
    break;
  case OBJ_BRANCH:
    num = net->num_branches;
    array = net->branch;
    get_element = &BRANCH_array_get;
    set_flags = &BRANCH_set_flags;
    has_properties = &BRANCH_has_properties;
    break;
  case OBJ_SHUNT:
    num = net->num_shunts;
    array = net->shunt;
    get_element = &SHUNT_array_get;
    set_flags = &SHUNT_set_flags;
    has_properties = &SHUNT_has_properties;
    break;
  case OBJ_VARGEN:
    num = net->num_vargens;
    array = net->vargen;
    get_element = &VARGEN_array_get;
    set_flags = &VARGEN_set_flags;
    has_properties = &VARGEN_has_properties;
    break;
  case OBJ_BAT:
    num = net->num_bats;
    array = net->bat;
    get_element = &BAT_array_get;
    set_flags = &BAT_set_flags;
    has_properties = &BAT_has_properties;
    break;
  default:
    sprintf(net->error_string,"invalid object type");
    net->error_flag = TRUE;
    return;
  }

  // Set flags
  for (i = 0; i < num; i++) {
    obj = get_element(array,i);
    if (has_properties(obj,prop_mask)) {
      if (flag_mask & FLAG_VARS)
	net->num_vars = set_flags(obj,FLAG_VARS,val_mask,net->num_vars);
      if (flag_mask & FLAG_FIXED)
	net->num_fixed = set_flags(obj,FLAG_FIXED,val_mask,net->num_fixed);
      if (flag_mask & FLAG_BOUNDED)
	net->num_bounded = set_flags(obj,FLAG_BOUNDED,val_mask,net->num_bounded);
      if (flag_mask & FLAG_SPARSE)
	net->num_sparse = set_flags(obj,FLAG_SPARSE,val_mask,net->num_sparse);
    }
  }
}

void NET_set_flags_of_component(Net* net, void* obj, char obj_type, char flag_mask, unsigned char val_mask) {

  // Local variables
  int (*set_flags)(void*,char,unsigned char,int);
  char (*get_obj_type)(void*);

  // Check
  if (!net)
    return;

  // Set pointers
  switch (obj_type) {
  case OBJ_BUS:
    set_flags = &BUS_set_flags;
    get_obj_type = &BUS_get_obj_type;
    break;
  case OBJ_GEN:
    set_flags = &GEN_set_flags;
    get_obj_type = &GEN_get_obj_type;
    break;
  case OBJ_LOAD:
    set_flags = &LOAD_set_flags;
    get_obj_type = &LOAD_get_obj_type;
    break;
  case OBJ_BRANCH:
    set_flags = &BRANCH_set_flags;
    get_obj_type = &BRANCH_get_obj_type;
    break;
  case OBJ_SHUNT:
    set_flags = &SHUNT_set_flags;
    get_obj_type = &SHUNT_get_obj_type;
    break;
  case OBJ_VARGEN:
    set_flags = &VARGEN_set_flags;
    get_obj_type = &VARGEN_get_obj_type;
    break;
  case OBJ_BAT:
    set_flags = &BAT_set_flags;
    get_obj_type = &BAT_get_obj_type;
    break;
  default:
    sprintf(net->error_string,"invalid object type");
    net->error_flag = TRUE;
    return;
  }

  // Check type
  if (obj_type != get_obj_type(obj)) {
    sprintf(net->error_string,"object type mismatch");
    net->error_flag = TRUE;
    return;
  }

  // Set flags
  if (flag_mask & FLAG_VARS)
    net->num_vars = set_flags(obj,FLAG_VARS,val_mask,net->num_vars);
  if (flag_mask & FLAG_FIXED)
    net->num_fixed = set_flags(obj,FLAG_FIXED,val_mask,net->num_fixed);
  if (flag_mask & FLAG_BOUNDED)
    net->num_bounded = set_flags(obj,FLAG_BOUNDED,val_mask,net->num_bounded);
  if (flag_mask & FLAG_SPARSE)
    net->num_sparse = set_flags(obj,FLAG_SPARSE,val_mask,net->num_sparse);
}

void NET_set_var_values(Net* net, Vec* values) {

  // Local variables
  int i;

  if (!net)
    return;

  // Buses
  for (i = 0; i < net->num_buses; i++)
    BUS_set_var_values(BUS_array_get(net->bus,i),values);

  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_set_var_values(GEN_array_get(net->gen,i),values);

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_set_var_values(BRANCH_array_get(net->branch,i),values);

  // Shunts
  for (i = 0; i < net->num_shunts; i++)
    SHUNT_set_var_values(SHUNT_array_get(net->shunt,i),values);

  // Loads
  for (i = 0; i < net->num_loads; i++)
    LOAD_set_var_values(LOAD_array_get(net->load,i),values);

  // Vargens
  for (i = 0; i < net->num_vargens; i++)
    VARGEN_set_var_values(VARGEN_array_get(net->vargen,i),values);

  // Batteries
  for (i = 0; i < net->num_bats; i++)
    BAT_set_var_values(BAT_array_get(net->bat,i),values);
}

char* NET_get_show_components_str(Net* net) {

  char* out;

  if (!net)
    return NULL;

  out = net->output_string;
  strcpy(out,"");

  sprintf(out+strlen(out),"\nNetwork Components\n");
  sprintf(out+strlen(out),"------------------\n");
  sprintf(out+strlen(out),"buses            : %d\n",NET_get_num_buses(net));
  sprintf(out+strlen(out),"  slack          : %d\n",NET_get_num_slack_buses(net));
  sprintf(out+strlen(out),"  reg by gen     : %d\n",NET_get_num_buses_reg_by_gen(net));
  sprintf(out+strlen(out),"  reg by tran    : %d\n",NET_get_num_buses_reg_by_tran(net));
  sprintf(out+strlen(out),"  reg by shunt   : %d\n",NET_get_num_buses_reg_by_shunt(net));
  sprintf(out+strlen(out),"shunts           : %d\n",NET_get_num_shunts(net));
  sprintf(out+strlen(out),"  fixed          : %d\n",NET_get_num_fixed_shunts(net));
  sprintf(out+strlen(out),"  switched v     : %d\n",NET_get_num_switched_shunts(net));
  sprintf(out+strlen(out),"branches         : %d\n",NET_get_num_branches(net));
  sprintf(out+strlen(out),"  lines          : %d\n",NET_get_num_lines(net));
  sprintf(out+strlen(out),"  fixed trans    : %d\n",NET_get_num_fixed_trans(net));
  sprintf(out+strlen(out),"  phase shifters : %d\n",NET_get_num_phase_shifters(net));
  sprintf(out+strlen(out),"  tap changers v : %d\n",NET_get_num_tap_changers_v(net));
  sprintf(out+strlen(out),"  tap changers Q : %d\n",NET_get_num_tap_changers_Q(net));
  sprintf(out+strlen(out),"generators       : %d\n",NET_get_num_gens(net));
  sprintf(out+strlen(out),"  slack          : %d\n",NET_get_num_slack_gens(net));
  sprintf(out+strlen(out),"  reg            : %d\n",NET_get_num_reg_gens(net));
  sprintf(out+strlen(out),"  P adjust       : %d\n",NET_get_num_P_adjust_gens(net));
  sprintf(out+strlen(out),"loads            : %d\n",NET_get_num_loads(net));
  sprintf(out+strlen(out),"  P adjust       : %d\n",NET_get_num_P_adjust_loads(net));
  sprintf(out+strlen(out),"vargens          : %d\n",NET_get_num_vargens(net));
  sprintf(out+strlen(out),"batteries        : %d\n",NET_get_num_bats(net));

  return out;
}

void NET_show_components(Net* net) {

  printf("%s",NET_get_show_components_str(net));
}

char* NET_get_show_properties_str(Net* net, int t) {

  char* out;

  if (!net || t < 0 || t >= net->num_periods)
    return NULL;

  out = net->output_string;
  strcpy(out,"");

  sprintf(out+strlen(out),"\nNetwork Properties (t = %d)\n",NET_get_num_periods(net));
  sprintf(out+strlen(out),"------------------\n");
  sprintf(out+strlen(out),"bus v max   : %.2f     (p.u.)\n",NET_get_bus_v_max(net,t));
  sprintf(out+strlen(out),"bus v min   : %.2f     (p.u.)\n",NET_get_bus_v_min(net,t));
  sprintf(out+strlen(out),"bus v vio   : %.2f     (p.u.)\n",NET_get_bus_v_vio(net,t));
  sprintf(out+strlen(out),"bus P mis   : %.2e (MW)\n",NET_get_bus_P_mis(net,t));
  sprintf(out+strlen(out),"bus Q mis   : %.2e (MVAr)\n",NET_get_bus_Q_mis(net,t));
  sprintf(out+strlen(out),"gen P cost  : %.2e ($/hr)\n",NET_get_gen_P_cost(net,t));
  sprintf(out+strlen(out),"gen v dev   : %.2e (p.u.)\n",NET_get_gen_v_dev(net,t));
  sprintf(out+strlen(out),"gen Q vio   : %.2e (MVAr)\n",NET_get_gen_Q_vio(net,t));
  sprintf(out+strlen(out),"gen P vio   : %.2e (MW)\n",NET_get_gen_P_vio(net,t));
  sprintf(out+strlen(out),"tran v vio  : %.2e (p.u.)\n",NET_get_tran_v_vio(net,t));
  sprintf(out+strlen(out),"tran r vio  : %.2e       \n",NET_get_tran_r_vio(net,t));
  sprintf(out+strlen(out),"tran p vio  : %.2e (rad)\n",NET_get_tran_p_vio(net,t));
  sprintf(out+strlen(out),"shunt v vio : %.2e (p.u.)\n",NET_get_shunt_v_vio(net,t));
  sprintf(out+strlen(out),"shunt b vio : %.2e (p.u.)\n",NET_get_shunt_b_vio(net,t));
  sprintf(out+strlen(out),"load P util : %.2e ($/hr)\n",NET_get_load_P_util(net,t));
  sprintf(out+strlen(out),"load P vio  : %.2e (MW)\n",NET_get_load_P_vio(net,t));
  sprintf(out+strlen(out),"num actions : %d\n",NET_get_num_actions(net,t));

  return out;
}

void NET_show_properties(Net* net, int t) {

  printf("%s",NET_get_show_properties_str(net,t));
}

void NET_show_buses(Net* net, int number, int sort_by, int t) {

  // Local variables
  Bus* bus;
  int type;
  REAL value;
  int counter;
  char type_string[100];
  char units[100];

  printf("\nTop Buses\n");
  printf("---------\n");

  if (BUS_SENS_LARGEST <= sort_by && sort_by <= BUS_SENS_V_REG_BY_SHUNT) { // sensitivity
    strcpy(type_string,"  sensitivity  ");
    strcpy(units,"  normalized  ");
    printf("%7s %7s %15s %9s %14s\n","  index"," number",type_string,"   value ","     units    ");
  }
  else if (BUS_MIS_LARGEST <= sort_by && sort_by <= BUS_MIS_REACTIVE) { // mismatch
    strcpy(type_string,"   mismatch  ");
    printf("%7s %7s %13s %9s %10s\n","  index"," number",type_string,"   value ","   units  ");
  }
  else {
    printf("invalid sort option\n");
    return;
  }

  counter = 0;
  bus = NET_create_sorted_bus_list(net,sort_by,t);
  while (bus != NULL && counter < number) {

    printf("%7d ",BUS_get_index(bus));
    printf("%7d ",BUS_get_number(bus));

    if (BUS_SENS_LARGEST <= sort_by && sort_by <= BUS_SENS_V_REG_BY_SHUNT) { // sensitivity

      value = BUS_get_quantity(bus,sort_by,t);
      if (sort_by == BUS_SENS_LARGEST)
	type = BUS_get_largest_sens_type(bus,t);
      else
	type = sort_by;

      switch (type) {
      case BUS_SENS_P_BALANCE:
	printf("%15s ","   P_balance   ");
	break;
      case BUS_SENS_Q_BALANCE:
	printf("%15s ","   Q_balance   ");
	break;
      case BUS_SENS_V_MAG_U_BOUND:
	printf("%15s "," v_mag_u_bound ");
	break;
      case BUS_SENS_V_MAG_L_BOUND:
	printf("%15s "," v_mag_l_bound ");
	break;
      case BUS_SENS_V_ANG_U_BOUND:
	printf("%15s "," v_ang_u_bound ");
	break;
      case BUS_SENS_V_ANG_L_BOUND:
	printf("%15s "," v_ang_l_bound ");
	break;
      case BUS_SENS_V_REG_BY_GEN:
	printf("%15s ","  v_reg_by_gen ");
	break;
      case BUS_SENS_V_REG_BY_TRAN:
	printf("%15s "," v_reg_by_tran ");
	break;
      case BUS_SENS_V_REG_BY_SHUNT:
	printf("%15s "," v_reg_by_shunt");
	break;
      default:
	printf("%15s ","    unknown    ");
      }
      printf("% 9.2e ",value);
      printf("%14s\n",units);
    }

    else if (BUS_MIS_LARGEST <= sort_by && sort_by <= BUS_MIS_REACTIVE) { // mismatch

      value = BUS_get_quantity(bus,sort_by,t);
      if (sort_by == BUS_MIS_LARGEST)
	type = BUS_get_largest_mis_type(bus,t);
      else
	type = sort_by;

      switch (type) {
      case BUS_MIS_ACTIVE:
	printf("%13s ","    active   ");
	strcpy(units,"    MW    ");
	break;
      case BUS_MIS_REACTIVE:
	printf("%13s ","   reactive  ");
	strcpy(units,"   MVAr   ");
	break;
      default:
	printf("%13s ","   unknown   ");
	strcpy(units,"  unknown ");
      }
      printf("% 9.2e ",value);
      printf("%10s\n",units);
    }

    bus = BUS_get_next(bus);
    counter++;
  }
}

void NET_update_properties(Net* net, Vec* values) {

  // Local variables
  int i;
  int t;

  // Clear
  NET_clear_properties(net);

  // Update
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      NET_update_properties_step(net,NET_get_branch(net,i),t,values);
  }
}

void NET_update_properties_step(Net* net, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Shunt* shunt;
  Bat* bat;

  REAL P;
  REAL Q;
  REAL dQ;
  REAL dP;

  REAL v[2];
  REAL dv;

  REAL a;
  REAL da;
  REAL phi;
  REAL dphi;

  REAL shunt_b;
  REAL shunt_db;
  REAL shunt_g;

  int k;
  int T;

  // Check pointers
  if (!net || !br)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Periods
  T = net->num_periods;

  // Voltage magnitudes
  for (k = 0; k < 2; k++) {
    bus = buses[k];
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && var_values)
      v[k] = VEC_get(var_values,BUS_get_index_v_mag(bus,t));
    else
      v[k] = BUS_get_v_mag(bus,t);
  }

  // Branch data
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) && var_values)
    a = VEC_get(var_values,BRANCH_get_index_ratio(br,t));
  else
    a = BRANCH_get_ratio(br,t);
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) && var_values)
    phi = VEC_get(var_values,BRANCH_get_index_phase(br,t));
  else
    phi = BRANCH_get_phase(br,t);

  // Tap ratios
  if (BRANCH_is_tap_changer(br)) {

    // Tap ratio limit violations
    //***************************
    da = 0;
    if (a > BRANCH_get_ratio_max(br))
      da = (a-BRANCH_get_ratio_max(br));
    if (a < BRANCH_get_ratio_min(br))
      da = (BRANCH_get_ratio_min(br)-a);
    if (da > net->tran_r_vio[t])
      net->tran_r_vio[t] = da;

    // Tap ratio actions
    //******************
    da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
    if (da < NET_CONTROL_EPS)
      da = NET_CONTROL_EPS;
    if (100.*fabs(a-BRANCH_get_ratio(br,t))/da > NET_CONTROL_ACTION_PCT)
      net->num_actions[t]++;
  }

  // Phase shifts
  if (BRANCH_is_phase_shifter(br)) {

    // Phase shift limit violations
    //*****************************
    dphi = 0;
    if (phi > BRANCH_get_phase_max(br))
      dphi = (phi-BRANCH_get_phase_max(br));
    if (phi < BRANCH_get_phase_min(br))
      dphi = (BRANCH_get_phase_min(br)-phi);
    if (dphi > net->tran_p_vio[t])
      net->tran_p_vio[t] = dphi;

    // Phase shift actions
    //********************
    dphi = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
    if (dphi < NET_CONTROL_EPS)
      dphi = NET_CONTROL_EPS;
    if (100.*fabs(phi-BRANCH_get_phase(br,t))/dphi > NET_CONTROL_ACTION_PCT)
      net->num_actions[t]++;
  }

  // Branch flows
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    // Update injected P,Q at buses k and m
    if (k == 0) {
      BUS_inject_P(bus,-BRANCH_get_P_km(br,var_values,t),t);
      BUS_inject_Q(bus,-BRANCH_get_Q_km(br,var_values,t),t);
    }
    else {
      BUS_inject_P(bus,-BRANCH_get_P_mk(br,var_values,t),t);
      BUS_inject_Q(bus,-BRANCH_get_Q_mk(br,var_values,t),t);
    }
  }

  // Other flows
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    // Skip if already counted
    if (net->bus_counted[BUS_get_index(bus)*T+t])
      continue;
    else
      net->bus_counted[BUS_get_index(bus)*T+t] = TRUE;

    // Maximum and minimum voltage magnitudes
    //***************************************
    if (net->bus_v_max[t] == 0 && net->bus_v_min[t] == 0) {
      net->bus_v_max[t] = v[k];
      net->bus_v_min[t] = v[k];
    }
    else {
      if (v[k] > net->bus_v_max[t])
	net->bus_v_max[t] = v[k];
      if (v[k] < net->bus_v_min[t])
	net->bus_v_min[t] = v[k];
    }

    // Normal voltage magnitude limit violations
    //************************************
    dv = 0;
    if (v[k] > BUS_get_v_max_norm(bus))
      dv = v[k]-BUS_get_v_max_norm(bus);
    if (v[k] < BUS_get_v_min_norm(bus))
      dv = BUS_get_v_min_norm(bus)-v[k];
    if (dv > net->bus_v_vio[t])
      net->bus_v_vio[t] = dv;

    // Regulation voltage magntiude limit violations
    //**********************************************
    dv = 0;
    if (v[k] > BUS_get_v_max_reg(bus))
      dv = v[k]-BUS_get_v_max_reg(bus);
    if (v[k] < BUS_get_v_min_reg(bus))
      dv = BUS_get_v_min_reg(bus)-v[k];
    if (BUS_is_regulated_by_tran(bus)) {
      if (dv > net->tran_v_vio[t])
	net->tran_v_vio[t] = dv;
    }
    if (BUS_is_regulated_by_shunt(bus)) {
      if (dv > net->shunt_v_vio[t])
	net->shunt_v_vio[t] = dv;
    }

    // Bus regulated by gen
    if (BUS_is_regulated_by_gen(bus)) {

      // Voltage set point deviation
      //****************************
      if (fabs(v[k]-BUS_get_v_set(bus,t)) > net->gen_v_dev[t])
	net->gen_v_dev[t] = fabs(v[k]-BUS_get_v_set(bus,t));

      // Voltage set point action
      //*************************
      dv = BUS_get_v_max_reg(bus)-BUS_get_v_min_reg(bus);
      if (dv < NET_CONTROL_EPS)
	dv = NET_CONTROL_EPS;
      if (100.*fabs(v[k]-BUS_get_v_set(bus,t))/dv > NET_CONTROL_ACTION_PCT)
	net->num_actions[t]++;
    }

    // Generators
    for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) && var_values)
	P = VEC_get(var_values,GEN_get_index_P(gen,t));
      else
	P = GEN_get_P(gen,t);
      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q) && var_values)
	Q = VEC_get(var_values,GEN_get_index_Q(gen,t));
      else
	Q = GEN_get_Q(gen,t);

      // Injections
      BUS_inject_P(bus,P,t);
      BUS_inject_Q(bus,Q,t);

      // Active power generation cost
      //*****************************
      net->gen_P_cost[t] += GEN_get_P_cost_for(gen,P);

      // Reacive power
      if (GEN_is_regulator(gen)) { // Should this be done for all generators?

	// Reactive power limit violations
	//********************************
	dQ = 0;
	if (Q > GEN_get_Q_max(gen))
	  dQ = (Q-GEN_get_Q_max(gen))*net->base_power; // MVAr
	if (Q < GEN_get_Q_min(gen))
	  dQ = (GEN_get_Q_min(gen)-Q)*net->base_power; // MVAr
	if (dQ > net->gen_Q_vio[t])
	  net->gen_Q_vio[t] = dQ;
      }

      // Active power limit violations
      //******************************
      dP = 0;
      if (P > GEN_get_P_max(gen))
	dP = (P-GEN_get_P_max(gen))*net->base_power; // MW
      if (P < GEN_get_P_min(gen))
	dP = (GEN_get_P_min(gen)-P)*net->base_power; // MW
      if (dP > net->gen_P_vio[t])
	net->gen_P_vio[t] = dP;

      // Active power actions
      //*********************
      dP = GEN_get_P_max(gen)-GEN_get_P_min(gen);
      if (dP < NET_CONTROL_EPS)
	dP = NET_CONTROL_EPS;
      if (100.*fabs(P-GEN_get_P(gen,t))/dP > NET_CONTROL_ACTION_PCT)
	net->num_actions[t]++;
    }

    // Loads
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

      if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P) && var_values)
	P = VEC_get(var_values,LOAD_get_index_P(load,t));
      else
	P = LOAD_get_P(load,t);
      if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q) && var_values)
	Q = VEC_get(var_values,LOAD_get_index_Q(load,t));
      else
	Q = LOAD_get_Q(load,t);

      // Injections
      BUS_inject_P(bus,-P,t);
      BUS_inject_Q(bus,-Q,t);

      // Active power consumption utility
      //*********************************
      net->load_P_util[t] += LOAD_get_P_util_for(load,P);

      // Active power limit violations
      //******************************
      dP = 0;
      if (P > LOAD_get_P_max(load,t))
	dP = (P-LOAD_get_P_max(load,t))*net->base_power; // MW
      if (P < LOAD_get_P_min(load,t))
	dP = (LOAD_get_P_min(load,t)-P)*net->base_power; // MW
      if (dP > net->load_P_vio[t])
	net->load_P_vio[t] = dP;

      // Active power actions
      //*********************
      dP = LOAD_get_P_max(load,t)-LOAD_get_P_min(load,t);
      if (dP < NET_CONTROL_EPS)
	dP = NET_CONTROL_EPS;
      if (100.*fabs(P-LOAD_get_P(load,t))/dP > NET_CONTROL_ACTION_PCT)
	net->num_actions[t]++;
    }

    // Batteries
    for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

      if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P) && var_values)
	P = VEC_get(var_values,BAT_get_index_Pc(bat,t))-VEC_get(var_values,BAT_get_index_Pd(bat,t));
      else
	P = BAT_get_P(bat,t);

      // Injections
      BUS_inject_P(bus,-P,t);
    }

    // Variable generators
    for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

      if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P) && var_values)
	P = VEC_get(var_values,VARGEN_get_index_P(vargen,t));
      else
	P = VARGEN_get_P(vargen,t);
      if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q) && var_values)
	Q = VEC_get(var_values,VARGEN_get_index_Q(vargen,t));
      else
	Q = VARGEN_get_Q(vargen,t);

      // Injections
      BUS_inject_P(bus,P,t);
      BUS_inject_Q(bus,Q,t);
    }

    // Shunts
    for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

      shunt_g = SHUNT_get_g(shunt);
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) && var_values)
	shunt_b = VEC_get(var_values,SHUNT_get_index_b(shunt,t));
      else
	shunt_b = SHUNT_get_b(shunt,t);

      // Flows
      BUS_inject_P(bus,-shunt_g*v[k]*v[k],t);
      BUS_inject_Q(bus,shunt_b*v[k]*v[k],t);

      // Switched shunts
      if (SHUNT_is_switched_v(shunt)) {

	// Switched shunt susceptance violations
	//**************************************
	shunt_db = 0;
	if (shunt_b > SHUNT_get_b_max(shunt))
	  shunt_db = (shunt_b-SHUNT_get_b_max(shunt));
	if (shunt_b < SHUNT_get_b_min(shunt))
	  shunt_db = (SHUNT_get_b_min(shunt)-shunt_b);
	if (shunt_db > net->shunt_b_vio[t])
	  net->shunt_b_vio[t] = shunt_db;

	// Swtiched shunt susceptance actions
	//***********************************
	shunt_db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt);
	if (shunt_db < NET_CONTROL_EPS)
	  shunt_db = NET_CONTROL_EPS;
	if (100.*fabs(shunt_b-SHUNT_get_b(shunt,t))/shunt_db > NET_CONTROL_ACTION_PCT)
	  net->num_actions[t]++;
      }
    }
  }

  // Power mismatches
  if (BRANCH_get_index(br) == net->num_branches-1) {
    BUS_array_get_max_mismatches(net->bus,
				 net->num_buses,
				 &(net->bus_P_mis[t]),
				 &(net->bus_Q_mis[t]),
				 t);
    net->bus_P_mis[t] *= net->base_power;
    net->bus_Q_mis[t] *= net->base_power;
  }
}

void NET_update_set_points(Net* net) {

  // Local variables
  Bus* bus;
  int i;
  int t;

  // No net
  if (!net)
    return;

  // Update
  for (i = 0; i < net->num_buses; i++) {
    bus = BUS_array_get(net->bus,i);
    if (BUS_is_regulated_by_gen(bus)) {
      for (t = 0; t < net->num_periods; t++)
	BUS_set_v_set(bus,BUS_get_v_mag(bus,t),t);
    }
  }
}

void NET_propagate_data_in_time(Net* net) {

  // Local variables
  int i;

  // No net
  if (!net)
    return;

  // Buses
  for (i = 0; i < net->num_buses; i++)
    BUS_propagate_data_in_time(BUS_array_get(net->bus,i));

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_propagate_data_in_time(BRANCH_array_get(net->branch,i));
  
  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_propagate_data_in_time(GEN_array_get(net->gen,i));

  // Loads
  for (i = 0; i < net->num_loads; i++)
    LOAD_propagate_data_in_time(LOAD_array_get(net->load,i));

  // Vargens
  for (i = 0; i < net->num_vargens; i++)
    VARGEN_propagate_data_in_time(VARGEN_array_get(net->vargen,i));

  // Shunts
  for (i = 0; i < net->num_shunts; i++)
    SHUNT_propagate_data_in_time(SHUNT_array_get(net->shunt,i));

  // Batteries
  for (i = 0; i < net->num_bats; i++)
    BAT_propagate_data_in_time(BAT_array_get(net->bat,i));
}
