/** @file net.c
 *  @brief This file defines the Net data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/net.h>
#include <pfnet/parser_RAW.h>
#include <pfnet/parser_MAT.h>
#include <pfnet/parser_ART.h>

struct Net {
  
  // Error
  BOOL error_flag;                    /**< @brief Error flag */
  char error_string[NET_BUFFER_SIZE]; /**< @brief Error string */
  
  // Components
  Bus* bus;             /**< @brief Bus array */
  Branch* branch;       /**< @brief Branch array */
  Gen* gen;             /**< @brief Gen array */
  Load* load;           /**< @brief Load array */
  Shunt* shunt;         /**< @brief Shunt array */
  Vargen* vargen;       /**< @brief Vargen array */

  // Hash tables
  Bus* bus_hash_number; /**< @brief Bus hash table indexed by bus numbers */
  Bus* bus_hash_name;   /**< @brief Bus hash table indexed by bus names */

  // Number of components
  int num_buses;     /**< @brief Number of buses (size of Bus array) */
  int num_branches;  /**< @brief Number of branches (size of Branch array) */
  int num_gens;      /**< @brief Number of generators (size of Gen array) */
  int num_loads;     /**< @brief Number of loads (size of Load array) */
  int num_shunts;    /**< @brief Number of shunts (size of Shunt array) */
  int num_vargens;   /**< @brief Number of variable generators (size of Vargen array) */

  // Number of flags
  int num_vars;      /**< @brief Number of variable quantities. */
  int num_fixed;     /**< @brief Number of fixed quantities. */
  int num_bounded;   /**< @brief Number of bounded quantities. */
  int num_sparse;    /**< @brief Number of sparse control quantities. */

  // Base power
  REAL base_power; /**< @brief System base power (MVA) */

  // Properties
  REAL bus_v_max;    /**< @brief Maximum bus voltage magnitude (p.u.). */
  REAL bus_v_min;    /**< @brief Minimum bus volatge magnitude (p.u.). */
  REAL bus_v_vio;    /**< @brief Maximum bus voltage limits violation (p.u.). */
  REAL bus_P_mis;    /**< @brief Maximum bus active power mismatch (MW). */
  REAL bus_Q_mis;    /**< @brief Maximum bus reactive power mismatch (MVAr). */
  REAL gen_v_dev;    /**< @brief Maximum generator voltage setpoint deviation (p.u.). */
  REAL gen_Q_vio;    /**< @brief Maximum generator reactive power limit violation (MVAr). */
  REAL gen_P_vio;    /**< @brief Maximum generator active power limit violation (MW). */
  REAL tran_v_vio;   /**< @brief Maximum transformer-controlled bus voltage magnitude band violation (p.u.). */
  REAL tran_r_vio;   /**< @brief Maximum tap ratio limit violation of tap-changing transformer (unitless). */
  REAL tran_p_vio;   /**< @brief Maximum phase shift limit violation of phase-shifting trasnformer (radians). */
  REAL shunt_v_vio;  /**< @brief Maximum shunt-controlled bus voltage mangnitude band violation (p.u.). */
  REAL shunt_b_vio;  /**< @brief Maximum susceptance limit volation of switched shunt device (p.u.). */
  int num_actions;   /**< @brief Number of control actions. */

  // Spatial correlation
  REAL vargen_corr_radius; /**< @brief Correlation radius for variable generators. **/
  REAL vargen_corr_value;  /**< @brief Correlation value for variable generators. **/
  
  // Utils
  char* bus_counted;  /**< @brief Flags for processing buses */
  int branch_counter; /**< @brief Counter for processing branches */
};

void NET_add_vargens(Net* net, Bus* bus_list, REAL penetration, REAL uncertainty, REAL corr_radius, REAL corr_value) {
  
  // Local variables
  REAL total_load_P;
  Vargen* vargen;
  int i;

  // Check
  if (!net)
    return;

  // Check
  if (penetration < 0 ||
      uncertainty < 0 ||
      corr_radius < 0 ||
      corr_value < -1 ||
      corr_value > 1) {
    sprintf(net->error_string,"invalid arguments for adding variable generators");
    net->error_flag = TRUE;
    return;
  }

  // Clear
  free(net->vargen);
  net->vargen = NULL;
  net->num_vargens = 0;
  
  // Save
  net->vargen_corr_radius = corr_radius;
  net->vargen_corr_value = corr_value;

  // Total load
  total_load_P = 0;
  for (i = 0; i < net->num_loads; i++)
    total_load_P += LOAD_get_P(NET_get_load(net,i));

  // Number
  net->num_vargens = BUS_list_len(bus_list);

  // Allocate
  net->vargen = VARGEN_array_new(net->num_vargens);

  // Set buses
  NET_set_vargen_buses(net,bus_list);

  // Set properties
  for (i = 0; i < net->num_vargens; i++) {
    vargen = NET_get_vargen(net,i);
    VARGEN_set_P_min(vargen,0.);
    VARGEN_set_P_max(vargen,total_load_P/net->num_vargens);
    VARGEN_set_P_std(vargen,(uncertainty/100.)*VARGEN_get_P_max(vargen));
    VARGEN_set_P(vargen,(penetration/100.)*VARGEN_get_P_max(vargen));
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
  int i;
  REAL num;
  REAL Ptot;
  REAL Qtot;
  REAL dQtot;
  REAL Q;
  REAL dQ;
  REAL Qmintot;
  REAL frac;

  if (!net)
    return;
  
  for (i = 0; i < net->num_buses; i++) {

    bus = NET_get_bus(net,i);

    // Slack gens
    if (BUS_is_slack(bus)) {
      num = 0;
      Ptot = 0;
      for(gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	Ptot += GEN_get_P(gen);
	num += 1;
      }
      for(gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen))
	GEN_set_P(gen,Ptot/num);
    }

    // Regulating gens
    if (BUS_is_regulated_by_gen(bus)) {
      Qtot = 0;
      dQtot = 0;
      Qmintot = 0;
      for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	Qtot += GEN_get_Q(gen);
	dQtot += GEN_get_Q_max(gen)-GEN_get_Q_min(gen);
	Qmintot += GEN_get_Q_min(gen);
      }
      gen = BUS_get_reg_gen(bus);
      dQ = GEN_get_Q_max(gen)-GEN_get_Q_min(gen);
      Q = GEN_get_Q_min(gen)+dQ*(Qtot-Qmintot)/dQtot;
      frac = (Q-GEN_get_Q_min(gen))/dQ;
      GEN_set_Q(gen,Q);
      for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen))
	GEN_set_Q(gen,GEN_get_Q_min(gen)+frac*(GEN_get_Q_max(gen)-GEN_get_Q_min(gen)));
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

  // Overall
  return (base_ok & bus_ok);

}

void NET_clear_data(Net* net) {
  if (net) {

    // Free data
    BUS_hash_number_del(net->bus_hash_number);
    BUS_hash_name_del(net->bus_hash_name);
    free(net->bus);
    free(net->branch);
    free(net->gen);
    free(net->load);
    SHUNT_array_free(net->shunt,net->num_shunts);
    free(net->vargen);
    free(net->bus_counted);

    // Re-initialize
    NET_init(net);
  }
}

void NET_clear_error(Net* net) {
  if (net) {
    net->error_flag = FALSE;
    strcpy(net->error_string,"");
  }
}

void NET_clear_flags(Net* net) {
  int i;
  Branch* br;
  Gen* gen;
  Bus* bus;
  Shunt* shunt;
  Vargen* vargen;

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

  // Vargens
  for (i = 0; i < net->num_vargens; i++) {
    vargen = VARGEN_array_get(net->vargen,i);
    VARGEN_clear_flags(vargen,FLAG_VARS);
    VARGEN_clear_flags(vargen,FLAG_FIXED);
    VARGEN_clear_flags(vargen,FLAG_BOUNDED);
    VARGEN_clear_flags(vargen,FLAG_SPARSE);
  }

  // Clear counters
  net->num_vars = 0;
  net->num_fixed = 0;
  net->num_bounded = 0;
  net->num_sparse = 0;
}

void NET_clear_properties(Net* net) {
  int i;
  if (net) {
    net->bus_v_max = 0;
    net->bus_v_min = 0;
    net->bus_v_vio = 0;
    net->bus_P_mis = 0;
    net->bus_Q_mis = 0;
    net->gen_v_dev = 0;
    net->gen_Q_vio = 0;
    net->gen_P_vio = 0;
    net->tran_v_vio = 0;
    net->tran_r_vio = 0;
    net->tran_p_vio = 0;
    net->shunt_v_vio = 0;
    net->shunt_b_vio = 0;
    net->num_actions = 0;
    if (net->bus_counted && net->bus) {
      for (i = 0; i < net->num_buses; i++) {
	BUS_clear_mismatches(BUS_array_get(net->bus,i));
	net->bus_counted[i] = 0;
      }
    }
    net->branch_counter = 0;
  }
}

void NET_clear_sensitivities(Net* net) {
  
  // Local variables
  Bus* bus;
  int i;

  if (!net)
    return;

  // Buses
  for (i = 0; i < net->num_buses; i++) 
    BUS_clear_sensitivities(BUS_array_get(net->bus,i));

  // Others (TO DO)
}

Bus* NET_create_sorted_bus_list(Net* net, int sort_by) {
  
  // Local variables
  Bus* bus_list = NULL;
  int i;
  
  if (!net)
    return bus_list;
  
  for (i = 0; i < net->num_buses; i++)
    bus_list = BUS_list_add_sorting(bus_list,BUS_array_get(net->bus,i),sort_by);
  return bus_list;
}
 
Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr) {
  /* This function constructs a covariance matrix for the active powers of
   * variable generators. The matrix is constructed such that the correlation 
   * coefficients of the (variable) active powers of vargens that are less than 
   * "spread" branches away is equal to "corr". Only the lower triangular part 
   * of the covaraicen matrix is stored. The resulting matrix should be checked 
   * to make sure it is a valid covariance matrix.
   */

  // Local variables
  Mat* sigma;
  Bus* bus_main;
  Bus* bus1;
  Bus* bus2;
  Branch* br;
  Vargen* vgen_main;
  Vargen* vg1;
  
  char* queued;
  int* neighbors;
  int neighbors_total;
  int neighbors_curr;
  int num_new;
  
  int nnz_counter;
  int i;
  int j;

  // Check
  if (!net)
    return NULL;

  // Allocate arrays
  queued = (char*)malloc(net->num_buses*sizeof(char));
  neighbors = (int*)malloc(net->num_buses*sizeof(int));

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

    // Add self to be processed
    neighbors_total = 1;
    neighbors[0] = BUS_get_index(bus_main);
    queued[BUS_get_index(bus_main)] = TRUE;

    // Neighbors
    neighbors_curr = 0;
    for (j = 0; j < spread; j++) {
      num_new = 0;
      while (neighbors_curr < neighbors_total) {
	bus1 = NET_get_bus(net,neighbors[neighbors_curr]);
	for (br = BUS_get_branch_from(bus1); br != NULL; br = BRANCH_get_from_next(br)) {
	  if (bus1 != BRANCH_get_bus_from(br)) {
	    sprintf(net->error_string,"unable to construct covariance matrix");
	    net->error_flag = TRUE;
	  }
	  bus2 = BRANCH_get_bus_to(br);
	  if (!queued[BUS_get_index(bus2)]) {
	    neighbors[neighbors_total+num_new] = BUS_get_index(bus2);
	    queued[BUS_get_index(bus2)] = TRUE;
	    num_new++;
	  }
	}
	for (br = BUS_get_branch_to(bus1); br != NULL; br = BRANCH_get_to_next(br)) {
	  if (bus1 != BRANCH_get_bus_to(br)) {
	    sprintf(net->error_string,"unable to construct covariance matrix");
	    net->error_flag = TRUE;
	  }
	  bus2 = BRANCH_get_bus_from(br);
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

    // Diagonal
    nnz_counter++;

    // Off diagonals
    for (j = 0; j < neighbors_total; j++) {
      bus1 = NET_get_bus(net,neighbors[j]);
      for (vg1 = BUS_get_vargen(bus1); vg1 != NULL; vg1 = VARGEN_get_next(vg1)) {
	if (VARGEN_has_flags(vg1,FLAG_VARS,VARGEN_VAR_P) &&
	    VARGEN_get_index_P(vgen_main) > VARGEN_get_index_P(vg1)) {
	  nnz_counter++;
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

    // Add self to be processed
    neighbors_total = 1;
    neighbors[0] = BUS_get_index(bus_main);
    queued[BUS_get_index(bus_main)] = TRUE;

    // Neighbors
    neighbors_curr = 0;
    for (j = 0; j < spread; j++) {
      num_new = 0;
      while (neighbors_curr < neighbors_total) {
	bus1 = NET_get_bus(net,neighbors[neighbors_curr]);
	for (br = BUS_get_branch_from(bus1); br != NULL; br = BRANCH_get_from_next(br)) {
	  if (bus1 != BRANCH_get_bus_from(br)) {
	    sprintf(net->error_string,"unable to construct covariance matrix");
	    net->error_flag = TRUE;
	  }
	  bus2 = BRANCH_get_bus_to(br);
	  if (!queued[BUS_get_index(bus2)]) {
	    neighbors[neighbors_total+num_new] = BUS_get_index(bus2);
	    queued[BUS_get_index(bus2)] = TRUE;
	    num_new++;
	  }
	}
	for (br = BUS_get_branch_to(bus1); br != NULL; br = BRANCH_get_to_next(br)) {
	  if (bus1 != BRANCH_get_bus_to(br)) {
	    sprintf(net->error_string,"unable to construct covariance matrix");
	    net->error_flag = TRUE;
	  }
	  bus2 = BRANCH_get_bus_from(br);
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

    // Diagonal
    MAT_set_i(sigma,nnz_counter,VARGEN_get_index_P(vgen_main));
    MAT_set_j(sigma,nnz_counter,VARGEN_get_index_P(vgen_main));
    MAT_set_d(sigma,nnz_counter,pow(VARGEN_get_P_std(vgen_main),2.));
    nnz_counter++;

    // Off diagonals
    for (j = 0; j < neighbors_total; j++) {
      bus1 = NET_get_bus(net,neighbors[j]);
      for (vg1 = BUS_get_vargen(bus1); vg1 != NULL; vg1 = VARGEN_get_next(vg1)) {
	if (VARGEN_has_flags(vg1,FLAG_VARS,VARGEN_VAR_P) &&
	    VARGEN_get_index_P(vgen_main) > VARGEN_get_index_P(vg1)) {
	  MAT_set_i(sigma,nnz_counter,VARGEN_get_index_P(vgen_main));
	  MAT_set_j(sigma,nnz_counter,VARGEN_get_index_P(vg1));
	  MAT_set_d(sigma,nnz_counter,
		    VARGEN_get_P_std(vgen_main)*VARGEN_get_P_std(vg1)*corr);
	  nnz_counter++;
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
 
  return sigma;
}

void NET_del(Net* net) {
  if (net) {    
    NET_clear_data(net);
    free(net);    
    net = NULL;
  } 
}

void NET_init(Net* net) {

  if (!net)
    return;

  // Error
  net->error_flag = FALSE;
  strcpy(net->error_string,"");

  // Components
  net->bus = NULL;
  net->branch = NULL;
  net->gen = NULL;
  net->load = NULL;
  net->shunt = NULL;
  net->vargen = NULL;

  // Hash tables
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;

  // Number components
  net->num_buses = 0;
  net->num_branches = 0;
  net->num_gens = 0;
  net->num_loads = 0;
  net->num_shunts = 0;
  net->num_vargens = 0;

  // Number flags 
  net->num_vars = 0;
  net->num_fixed = 0;
  net->num_bounded = 0;
  net->num_sparse = 0;

  // Base
  net->base_power = NET_BASE_POWER;

  // Spatial correlation
  net->vargen_corr_radius = 1;
  net->vargen_corr_value = 0;

  // Properties
  net->bus_v_max = 0;
  net->bus_v_min = 0;
  net->bus_v_vio = 0;
  net->bus_P_mis = 0;
  net->bus_Q_mis = 0;
  net->gen_v_dev = 0;
  net->gen_Q_vio = 0;
  net->gen_P_vio = 0;
  net->tran_v_vio = 0;
  net->tran_r_vio = 0;
  net->tran_p_vio = 0;
  net->shunt_v_vio = 0;
  net->shunt_b_vio = 0;
  net->num_actions = 0;

  // Utils
  net->bus_counted = NULL;
  net->branch_counter = 0;
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

int NET_get_num_buses(Net* net) {
  if (!net)
    return 0;
  else
    return net->num_buses;
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

REAL NET_get_total_gen_P(Net* net) {
  int i;
  REAL P = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_gens; i++) {
    P += GEN_get_P(GEN_array_get(net->gen,i)); // p.u.
  }
  return P*net->base_power; // MW
}

REAL NET_get_total_gen_Q(Net* net) {
  int i;
  REAL Q = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_gens; i++) {
    Q += GEN_get_Q(GEN_array_get(net->gen,i)); // p.u.
  }
  return Q*net->base_power; // MVAr
}

REAL NET_get_total_load_P(Net* net) {
  int i;
  REAL P = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_loads; i++) {
    P += LOAD_get_P(LOAD_array_get(net->load,i)); // p.u.
  }
  return P*net->base_power; // MW
}

REAL NET_get_total_load_Q(Net* net) {
  int i;
  REAL Q = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_loads; i++) {
    Q += LOAD_get_Q(LOAD_array_get(net->load,i)); // p.u.
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

  // Variable generators
  for (i = 0; i < net->num_vargens; i++) 
    VARGEN_get_var_values(VARGEN_array_get(net->vargen,i),values,code);

  // Return
  return values;  
}

Mat* NET_get_var_projection(Net* net, char obj_type, char var) {

  // Local variables
  int num_subvars;
  Mat* proj;
  int i;

  int num;
  void* obj;
  void* array;
  void* (*get_element)(void* array, int index);
  BOOL (*has_flags)(void*,char,char);
  int (*get_var_index)(void*,char);

  // Check
  if (!net)
    return NULL;

  // Set pointers
  switch (obj_type) {
  case OBJ_BUS:
    num = net->num_buses;
    array = net->bus;
    get_element = &BUS_array_get;
    has_flags = &BUS_has_flags;
    get_var_index = &BUS_get_var_index;
    break;
  case OBJ_GEN:
    num = net->num_gens;
    array = net->gen;
    get_element = &GEN_array_get;
    has_flags = &GEN_has_flags;
    get_var_index = &GEN_get_var_index;
    break;
  case OBJ_BRANCH:
    num = net->num_branches;
    array = net->branch;
    get_element = &BRANCH_array_get;
    has_flags = &BRANCH_has_flags;
    get_var_index = &BRANCH_get_var_index;
    break;
  case OBJ_SHUNT:
    num = net->num_shunts;
    array = net->shunt;
    get_element = &SHUNT_array_get;
    has_flags = &SHUNT_has_flags;
    get_var_index = &SHUNT_get_var_index;
    break;
  case OBJ_VARGEN:
    num = net->num_vargens;
    array = net->vargen;
    get_element = &VARGEN_array_get;
    has_flags = &VARGEN_has_flags;
    get_var_index = &VARGEN_get_var_index;
    break;
  default:
    sprintf(net->error_string,"invalid object type");
    net->error_flag = TRUE;
    return NULL;
  }
    
  // Count
  num_subvars = 0;
  for (i = 0; i < num; i++) {
    obj = get_element(array,i);
    if (has_flags(obj,FLAG_VARS,var))
      num_subvars++;
  }

  // Allocate
  proj = MAT_new(num_subvars,
		 net->num_vars,
		 num_subvars);
  
  // Fill
  num_subvars = 0;
  for (i = 0; i < num; i++) {
    obj = get_element(array,i);
    if (has_flags(obj,FLAG_VARS,var)) {
      MAT_set_i(proj,num_subvars,num_subvars);
      MAT_set_j(proj,num_subvars,get_var_index(obj,var));
      MAT_set_d(proj,num_subvars,1.);
      num_subvars++;
    }
  }
       
  // Return
  return proj;
}

REAL NET_get_bus_v_max(Net* net) {
  if (net)
    return net->bus_v_max;
  else
    return 0;
}

REAL NET_get_bus_v_min(Net* net) {
  if (net)
    return net->bus_v_min;
  else
    return 0;
}

REAL NET_get_bus_v_vio(Net* net) {
  if (net)
    return net->bus_v_vio;
  else
    return 0;
}

REAL NET_get_bus_P_mis(Net* net) {
  if (net)
    return net->bus_P_mis;
  else
    return 0;
}

REAL NET_get_bus_Q_mis(Net* net) {
  if (net)
    return net->bus_Q_mis;
  else
    return 0;
}

REAL NET_get_gen_v_dev(Net* net) {
  if (net)
    return net->gen_v_dev;
  else
    return 0;
}

REAL NET_get_gen_Q_vio(Net* net) {
  if (net)
    return net->gen_Q_vio;
  else
    return 0;
}

REAL NET_get_gen_P_vio(Net* net) {
  if (net)
    return net->gen_P_vio;
  else
    return 0;
}

REAL NET_get_tran_v_vio(Net* net) {
  if (net)
    return net->tran_v_vio;
  else
    return 0;
}

REAL NET_get_tran_r_vio(Net* net) {
  if (net)
    return net->tran_r_vio;
  else
    return 0;
}

REAL NET_get_tran_p_vio(Net* net) {
  if (net)
    return net->tran_p_vio;
  else
    return 0;
}

REAL NET_get_shunt_v_vio(Net* net) {
  if (net)
    return net->shunt_v_vio;
  else
    return 0;
}

REAL NET_get_shunt_b_vio(Net* net) {
  if (net)
    return net->shunt_b_vio;
  else
    return 0;
}

int NET_get_num_actions(Net* net) {
  if (net)
    return net->num_actions;
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

BOOL NET_has_error(Net* net) {
  if (!net)
    return FALSE;
  else
    return net->error_flag;
}

void NET_load(Net* net, char* filename) {

  // Local variables
  char* ext;

  // No network
  if (!net)
    return;

  // Free data
  NET_clear_data(net);

  // Extension
  ext = strrchr(filename,'.');

  // Parse
  if (!ext) {
    sprintf(net->error_string,"unable to find extension in %s",filename);
    net->error_flag = TRUE;
  }
  else if (strcmp(ext+1,"raw") == 0) {
  
    // PSSE raw
    RAW_Parser* parser = RAW_PARSER_new();
    RAW_PARSER_read(parser,filename);
    #ifdef DEBUG
      RAW_PARSER_show(parser);
    #endif
    if (!RAW_PARSER_has_error(parser))
      RAW_PARSER_load(parser,net);
    if (RAW_PARSER_has_error(parser)) {
      strcpy(net->error_string,RAW_PARSER_get_error_string(parser));
      net->error_flag = TRUE;
    }
    RAW_PARSER_del(parser);
  }
  else if (strcmp(ext+1,"mat") == 0) {
  
    // MATPOWER MAT
    MAT_Parser* parser = MAT_PARSER_new();
    MAT_PARSER_read(parser,filename);
    #ifdef DEBUG
      MAT_PARSER_show(parser);
    #endif
    if (!MAT_PARSER_has_error(parser))
      MAT_PARSER_load(parser,net);
    if (MAT_PARSER_has_error(parser)) {
      strcpy(net->error_string,MAT_PARSER_get_error_string(parser));
      net->error_flag = TRUE;
    }
    MAT_PARSER_del(parser);
  }
  else if (strcmp(ext+1,"art") == 0) {
  
    // ARTERE ART
    ART_Parser* parser = ART_PARSER_new();
    ART_PARSER_read(parser,filename);
    #ifdef DEBUG
      ART_PARSER_show(parser);
    #endif
    if (!ART_PARSER_has_error(parser))
      ART_PARSER_load(parser,net);
    if (ART_PARSER_has_error(parser)) {
      strcpy(net->error_string,ART_PARSER_get_error_string(parser));
      net->error_flag = TRUE;
    }
    ART_PARSER_del(parser);
  }
  else {
    sprintf(net->error_string,"invalid file type (%s)",ext+1);
    net->error_flag = TRUE;
  }

  // Set up utilities
  net->bus_counted = (char*)calloc(net->num_buses,sizeof(char));

  // Properties
  NET_update_properties(net,NULL);
}

Net* NET_new(void) {
  Net* net = (Net*)malloc(sizeof(Net));
  NET_init(net);
  return net;
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
  }
}

void NET_set_gen_array(Net* net, Gen* gen, int num) {
  if (net) {
    net->gen = gen;
    net->num_gens = num;
  }
}

void NET_set_vargen_array(Net* net, Vargen* gen, int num) {
  if (net) {
    net->vargen = gen;
    net->num_vargens = num;
  }
}

void NET_set_vargen_buses(Net* net, Bus* bus_list) {
  
  int i;
  Bus* bus;
  Vargen* gen;

  if (!net)
    return;

  // Clear
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

void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, char val_mask) {

  // Local variables
  int i;
  int num;
  void* obj;
  void* array;
  void* (*get_element)(void* array, int index);
  int (*set_flags)(void*,char,char,int);
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

  // Vargens
  for (i = 0; i < net->num_vargens; i++) 
    VARGEN_set_var_values(VARGEN_array_get(net->vargen,i),values);
}

void NET_show_components(Net *net) {
 
  printf("\nNetwork Components\n");
  printf("------------------\n");
  printf("buses            : %d\n",NET_get_num_buses(net));
  printf("  slack          : %d\n",NET_get_num_slack_buses(net));
  printf("  reg by gen     : %d\n",NET_get_num_buses_reg_by_gen(net));
  printf("  reg by tran    : %d\n",NET_get_num_buses_reg_by_tran(net));
  printf("  reg by shunt   : %d\n",NET_get_num_buses_reg_by_shunt(net));
  printf("shunts           : %d\n",NET_get_num_shunts(net));
  printf("  fixed          : %d\n",NET_get_num_fixed_shunts(net));
  printf("  switched v     : %d\n",NET_get_num_switched_shunts(net));
  printf("branches         : %d\n",NET_get_num_branches(net));
  printf("  lines          : %d\n",NET_get_num_lines(net));
  printf("  fixed trans    : %d\n",NET_get_num_fixed_trans(net));
  printf("  phase shifters : %d\n",NET_get_num_phase_shifters(net));
  printf("  tap changers v : %d\n",NET_get_num_tap_changers_v(net));
  printf("  tap changers Q : %d\n",NET_get_num_tap_changers_Q(net));
  printf("generators       : %d\n",NET_get_num_gens(net));
  printf("  slack          : %d\n",NET_get_num_slack_gens(net));
  printf("  reg            : %d\n",NET_get_num_reg_gens(net));
  printf("loads            : %d\n",NET_get_num_loads(net));
  printf("vargens          : %d\n",NET_get_num_vargens(net));
}

void NET_show_properties(Net* net) {
  
  printf("\nNetwork Properties\n");
  printf("------------------\n");
  printf("bus v max   : %.2f     (p.u.)\n",NET_get_bus_v_max(net));
  printf("bus v min   : %.2f     (p.u.)\n",NET_get_bus_v_min(net));
  printf("bus v vio   : %.2f     (p.u.)\n",NET_get_bus_v_vio(net));
  printf("bus P mis   : %.2e (MW)\n",NET_get_bus_P_mis(net));
  printf("bus Q mis   : %.2e (MVAr)\n",NET_get_bus_Q_mis(net));
  printf("gen v dev   : %.2e (p.u.)\n",NET_get_gen_v_dev(net));
  printf("gen Q vio   : %.2e (MVAr)\n",NET_get_gen_Q_vio(net));
  printf("gen P vio   : %.2e (MW)\n",NET_get_gen_P_vio(net));
  printf("tran v vio  : %.2e (p.u.)\n",NET_get_tran_v_vio(net));
  printf("tran r vio  : %.2e       \n",NET_get_tran_r_vio(net));
  printf("tran p vio  : %.2e (rad)\n",NET_get_tran_p_vio(net));
  printf("shunt v vio : %.2e (p.u.)\n",NET_get_shunt_v_vio(net));
  printf("shunt b vio : %.2e (p.u.)\n",NET_get_shunt_b_vio(net));
  printf("num actions : %d\n",NET_get_num_actions(net));
}

void NET_show_buses(Net* net, int number, int sort_by) {
  
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
  bus = NET_create_sorted_bus_list(net,sort_by);
  while (bus != NULL && counter < number) {
  
    printf("%7d ",BUS_get_index(bus));
    printf("%7d ",BUS_get_number(bus));
    
    if (BUS_SENS_LARGEST <= sort_by && sort_by <= BUS_SENS_V_REG_BY_SHUNT) { // sensitivity
      
      value = BUS_get_quantity(bus,sort_by);
      if (sort_by == BUS_SENS_LARGEST)
	type = BUS_get_largest_sens_type(bus);
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
      
      value = BUS_get_quantity(bus,sort_by);
      if (sort_by == BUS_MIS_LARGEST)
	type = BUS_get_largest_mis_type(bus);
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

  // Clear
  NET_clear_properties(net);

  // Update
  for (i = 0; i < NET_get_num_branches(net); i++)
    NET_update_properties_branch(net,NET_get_branch(net,i),values);
 }

void NET_update_properties_branch(Net* net, Branch* br, Vec* var_values) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Shunt* shunt;
  
  REAL P;
  REAL Q;
  REAL dQ;
  REAL dP;

  REAL v[2];
  REAL w[2];
  REAL dv;

  REAL a;
  REAL da;
  REAL a_temp;
  REAL phi;
  REAL dphi;
  REAL phi_temp;
  
  REAL b;
  REAL b_sh[2];

  REAL g;
  REAL g_sh[2];

  REAL flowP[2];
  REAL flowP_sh[2];
  REAL flowQ[2];
  REAL flowQ_sh[2];

  REAL shunt_b;
  REAL shunt_db;
  REAL shunt_g;

  int k;
  int m;
  
  // Check data
  if (!net || !br)
    return;

  // Branch counter
  net->branch_counter++;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++) {
    bus = buses[k];
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG) && var_values)
      w[k] = VEC_get(var_values,BUS_get_index_v_ang(bus));
    else
      w[k] = BUS_get_v_ang(bus);
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && var_values)
      v[k] = VEC_get(var_values,BUS_get_index_v_mag(bus));
    else
      v[k] = BUS_get_v_mag(bus);
  }
  
  // Branch data
  b = BRANCH_get_b(br);
  b_sh[0] = BRANCH_get_b_from(br);
  b_sh[1] = BRANCH_get_b_to(br);
  g = BRANCH_get_g(br);
  g_sh[0] = BRANCH_get_g_from(br);
  g_sh[1] = BRANCH_get_g_to(br);
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) && var_values)
    a = VEC_get(var_values,BRANCH_get_index_ratio(br));
  else
    a = BRANCH_get_ratio(br);
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) && var_values)
    phi = VEC_get(var_values,BRANCH_get_index_phase(br));
  else
    phi = BRANCH_get_phase(br);

  // Tap ratios
  if (BRANCH_is_tap_changer(br)) {

    // Tap ratio limit violations
    //***************************
    da = 0;
    if (a > BRANCH_get_ratio_max(br))
      da = (a-BRANCH_get_ratio_max(br));
    if (a < BRANCH_get_ratio_min(br)) 
      da = (BRANCH_get_ratio_min(br)-a);
    if (da > net->tran_r_vio)
      net->tran_r_vio = da;

    // Tap ratio actions
    //******************
    da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
    if (da < NET_CONTROL_EPS)
      da = NET_CONTROL_EPS;
    if (100.*fabs(a-BRANCH_get_ratio(br))/da > NET_CONTROL_ACTION_PCT)
      net->num_actions++;
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
    if (dphi > net->tran_p_vio)
      net->tran_p_vio = dphi;

    // Phase shift actions
    //********************
    dphi = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
    if (dphi < NET_CONTROL_EPS)
      dphi = NET_CONTROL_EPS;
    if (100.*fabs(phi-BRANCH_get_phase(br))/dphi > NET_CONTROL_ACTION_PCT)
      net->num_actions++;
  }

  // Branch flows
  for (k = 0; k < 2; k++) {
    bus = buses[k];
    if (k == 0) {
      m = 1;
      a_temp = a;
      phi_temp = phi;
    }
    else {
      m = 0;
      a_temp = 1;
      phi_temp = -phi;
    }
    flowP[k] = -a*v[k]*v[m]*(g*cos(w[k]-w[m]-phi_temp)+b*sin(w[k]-w[m]-phi_temp));
    flowQ[k] = -a*v[k]*v[m]*(g*sin(w[k]-w[m]-phi_temp)-b*cos(w[k]-w[m]-phi_temp));
    flowP_sh[k] =  a_temp*a_temp*(g_sh[k]+g)*v[k]*v[k];
    flowQ_sh[k] = -a_temp*a_temp*(b_sh[k]+b)*v[k]*v[k];

    // Flows
    BUS_inject_P(bus,-flowP_sh[k]-flowP[k]); 
    BUS_inject_Q(bus,-flowQ_sh[k]-flowQ[k]);
  }
  
  // Other flows
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];
    
    // Skip if already counted
    if (net->bus_counted[BUS_get_index(bus)])
      continue;
    else
      net->bus_counted[BUS_get_index(bus)] = TRUE;
    
    // Maximum and minimum voltage magnitudes
    //***************************************
    if (net->bus_v_max == 0 && net->bus_v_min == 0) {
      net->bus_v_max = v[k];
      net->bus_v_min = v[k];
    }
    else {
      if (v[k] > net->bus_v_max)
	net->bus_v_max = v[k];
      if (v[k] < net->bus_v_min)
	net->bus_v_min = v[k];
    }

    // Voltage limit violations
    //*************************
    dv = 0;
    if (v[k] > BUS_get_v_max(bus))
      dv = v[k]-BUS_get_v_max(bus);
    if (v[k] < BUS_get_v_min(bus))
      dv = BUS_get_v_min(bus)-v[k];
    if (dv > net->bus_v_vio)
      net->bus_v_vio = dv;
    
    // Tran-controlled
    if (BUS_is_regulated_by_tran(bus)) {
      if (dv > net->tran_v_vio)
	net->tran_v_vio = dv;
    }
    
    // Shunt-controlled
    if (BUS_is_regulated_by_shunt(bus)) {
      if (dv > net->shunt_v_vio)
	net->shunt_v_vio = dv;
    }

    // Bus regulated by gen
    if (BUS_is_regulated_by_gen(bus)) {

      // Voltage set point deviation
      //****************************
      if (fabs(v[k]-BUS_get_v_set(bus)) > net->gen_v_dev)
	net->gen_v_dev = fabs(v[k]-BUS_get_v_set(bus));
      
      // Voltage set point action
      //*************************
      dv = BUS_get_v_max(bus)-BUS_get_v_min(bus);
      if (dv < NET_CONTROL_EPS)
	dv = NET_CONTROL_EPS;
      if (100.*fabs(v[k]-BUS_get_v_set(bus))/dv > NET_CONTROL_ACTION_PCT)
	net->num_actions++;
    }

    // Generators
    for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) && var_values)
	P = VEC_get(var_values,GEN_get_index_P(gen));
      else
	P = GEN_get_P(gen);
      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q) && var_values)
	Q = VEC_get(var_values,GEN_get_index_Q(gen));
      else
	Q = GEN_get_Q(gen);

      // Injections
      BUS_inject_P(bus,P);
      BUS_inject_Q(bus,Q);
      
      // Reacive power
      if (GEN_is_regulator(gen)) { // Should this be done for all generators?

	// Reactive power limit violations
	//********************************
	dQ = 0;
	if (Q > GEN_get_Q_max(gen))
	  dQ = (Q-GEN_get_Q_max(gen))*net->base_power; // MVAr
	if (Q < GEN_get_Q_min(gen)) 
	  dQ = (GEN_get_Q_min(gen)-Q)*net->base_power; // MVAr
	if (dQ > net->gen_Q_vio)
	  net->gen_Q_vio = dQ;
      }

      // Active power limit violations
      //******************************
      dP = 0;
      if (P > GEN_get_P_max(gen))
	dP = (P-GEN_get_P_max(gen))*net->base_power; // MW
      if (P < GEN_get_P_min(gen)) 
	dP = (GEN_get_P_min(gen)-P)*net->base_power; // MW
      if (dP > net->gen_P_vio)
	net->gen_P_vio = dP;

      // Active power actions
      //*********************
      dP = GEN_get_P_max(gen)-GEN_get_P_min(gen);
      if (dP < NET_CONTROL_EPS)
	dP = NET_CONTROL_EPS;
      if (100.*fabs(P-GEN_get_P(gen))/dP > NET_CONTROL_ACTION_PCT)
	net->num_actions++;
    }

    // Variable generators
    for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

      if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P) && var_values)
	P = VEC_get(var_values,VARGEN_get_index_P(vargen));
      else
	P = VARGEN_get_P(vargen);
      if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q) && var_values)
	Q = VEC_get(var_values,VARGEN_get_index_Q(vargen));
      else
	Q = VARGEN_get_Q(vargen);

      // Injections
      BUS_inject_P(bus,P);
      BUS_inject_Q(bus,Q);
    }

    // Loads
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

      // Injections
      BUS_inject_P(bus,-LOAD_get_P(load));
      BUS_inject_Q(bus,-LOAD_get_Q(load));
    }

    // Shunts
    for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

      shunt_g = SHUNT_get_g(shunt);
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) && var_values)
	shunt_b = VEC_get(var_values,SHUNT_get_index_b(shunt));
      else
	shunt_b = SHUNT_get_b(shunt);

      // Flows
      BUS_inject_P(bus,-shunt_g*v[k]*v[k]);
      BUS_inject_Q(bus,shunt_b*v[k]*v[k]);
      

      // Switched shunts
      if (SHUNT_is_switched_v(shunt)) {

	// Switched shunt susceptance violations
	//**************************************
	shunt_db = 0;
	if (shunt_b > SHUNT_get_b_max(shunt))
	  shunt_db = (shunt_b-SHUNT_get_b_max(shunt));
	if (shunt_b < SHUNT_get_b_min(shunt)) 
	  shunt_db = (SHUNT_get_b_min(shunt)-shunt_b);
	if (shunt_db > net->shunt_b_vio)
	  net->shunt_b_vio = shunt_db;

	// Swtiched shunt susceptance actions
	//***********************************
	shunt_db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt);
	if (shunt_db < NET_CONTROL_EPS)
	  shunt_db = NET_CONTROL_EPS;
	if (100.*fabs(shunt_b-SHUNT_get_b(shunt))/shunt_db > NET_CONTROL_ACTION_PCT)
	  net->num_actions++;
      }
    }    
  }

  // Power mismatches
  if (net->branch_counter == net->num_branches) {
    BUS_array_get_max_mismatches(net->bus,
				 net->num_buses,
				 &(net->bus_P_mis),
				 &(net->bus_Q_mis));
    net->bus_P_mis *= net->base_power;
    net->bus_Q_mis *= net->base_power;
  }
}

void NET_update_set_points(Net* net) {
  
  // Local variables
  Bus* bus;
  int i;

  // Check
  if (!net)
    return;
  
  // Update
  for (i = 0; i < net->num_buses; i++) {
    bus = BUS_array_get(net->bus,i);
    if (BUS_is_regulated_by_gen(bus))
      BUS_set_v_set(bus,BUS_get_v_mag(bus));
  }
}
