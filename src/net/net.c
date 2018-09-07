/** @file net.c
 *  @brief This file defines the Net data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/net.h>
#include <pfnet/array.h>
#include <pfnet/json_macros.h>
#include <pfnet/pfnet_config.h>

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
  ConvCSC* csc_conv;    /**< @brief CSC converter array */
  ConvVSC* vsc_conv;    /**< @brief VSC converter array */
  BusDC* dc_bus;        /**< @brief DC bus array */
  BranchDC* dc_branch;  /**< @brief DC branch array */
  Facts* facts;         /**< @brief Facts array */

  // Hash tables
  Bus* bus_hash_number;      /**< @brief Bus hash table indexed by bus numbers. */
  Bus* bus_hash_name;        /**< @brief Bus hash table indexed by bus names. */
  BusDC* dc_bus_hash_number; /**< @brief DC bus hash table indexed by DC bus numbers. */
  BusDC* dc_bus_hash_name;   /**< @brief DC bus hash table indexed by DC bus names. */

  // Number of components
  int num_buses;     /**< @brief Number of buses (size of Bus array). */
  int num_branches;  /**< @brief Number of branches (size of Branch array). */
  int num_gens;      /**< @brief Number of generators (size of Gen array). */
  int num_loads;     /**< @brief Number of loads (size of Load array). */
  int num_shunts;    /**< @brief Number of shunts (size of Shunt array). */
  int num_vargens;   /**< @brief Number of variable generators (size of Vargen array). */
  int num_bats;      /**< @brief Number of batteries (size of Bat array). */
  int num_csc_convs;   /**< @brief Number of CSC converters (size of ConvCSC array) */
  int num_vsc_convs;   /**< @brief Number of VSC converters (size of ConvVSC array) */
  int num_dc_buses;    /**< @brief Number of DC buses (size of BusDC array) */
  int num_dc_branches; /**< @brief Number of DC branches (size of BranchDC array) */
  int num_facts;       /**< @brief Number of FACTS devices (size of Facts array) */

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
  REAL vargen_corr_radius; /**< @brief Correlation radius for variable generators. */
  REAL vargen_corr_value;  /**< @brief Correlation value for variable generators. */

  // State tag
  unsigned long int state_tag; /**< @brief State tag. */
};

void NET_inc_state_tag(Net* net) {
  if (net)
    net->state_tag++;
}

unsigned long int NET_get_state_tag(Net* net) {
  if (net)
    return net->state_tag;
  else
    return 0;
}

void NET_add_buses(Net* net, Bus** bus_ptr_array, int size) {
  /** Adds buses to the network. The entire bus array is
   * relocated, the data is copied (except flags for the new buses), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Bus* bus_src;
  Bus* bus_dst;
  Bus* bus;
  Bus* old_bus_array; 
  int old_num_buses;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bus_ptr_array)
    return;

  // Old buses
  old_bus_array = net->bus;
  old_num_buses = net->num_buses;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    bus = bus_ptr_array[i];
    if (bus != NET_get_bus(net,BUS_get_index(bus))) // not in the network
      num++;
    else
      bus_ptr_array[i] = NULL;                      // clear to ignore below
  }

  // New buses
  NET_set_bus_array(net,BUS_array_new(net->num_buses+num,net->num_periods),net->num_buses+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_buses+size; i++) {

    if (i < old_num_buses)
      bus_src = BUS_array_get(old_bus_array,i);
    else
      bus_src = bus_ptr_array[i-old_num_buses];
    bus_dst = NET_get_bus(net,index);

    // Check
    if (!bus_src)
      continue;

    // Copy data (except index, hash info, and connections)
    BUS_copy_from_bus(bus_dst,bus_src);
    
    // Clear flags
    if (i >= old_num_buses) {
      BUS_clear_flags(bus_dst,FLAG_VARS);
      BUS_clear_flags(bus_dst,FLAG_FIXED);
      BUS_clear_flags(bus_dst,FLAG_BOUNDED);
      BUS_clear_flags(bus_dst,FLAG_SPARSE);
    }

    // Connections - gen
    while (BUS_get_gen(bus_src))
      GEN_set_bus(BUS_get_gen(bus_src),bus_dst);

    // Connections - reg gen
    while (BUS_get_reg_gen(bus_src))
      GEN_set_reg_bus(BUS_get_reg_gen(bus_src),bus_dst);
    
    // Connection - branch k
    while (BUS_get_branch_k(bus_src))
      BRANCH_set_bus_k(BUS_get_branch_k(bus_src),bus_dst);

    // Connection - branch m
    while (BUS_get_branch_m(bus_src))
      BRANCH_set_bus_m(BUS_get_branch_m(bus_src),bus_dst);

    // Connection - reg tran
    while (BUS_get_reg_tran(bus_src))
      BRANCH_set_reg_bus(BUS_get_reg_tran(bus_src),bus_dst);

    // Connections - shunt
    while (BUS_get_shunt(bus_src))
      SHUNT_set_bus(BUS_get_shunt(bus_src),bus_dst);

    // Connections - reg shunt
    while (BUS_get_reg_shunt(bus_src))
      SHUNT_set_reg_bus(BUS_get_reg_shunt(bus_src),bus_dst);

    // Connections - Load
    while (BUS_get_load(bus_src))
      LOAD_set_bus(BUS_get_load(bus_src),bus_dst);

    // Connections - vargen
    while (BUS_get_vargen(bus_src))
      VARGEN_set_bus(BUS_get_vargen(bus_src),bus_dst);

    // Connections - bat
    while (BUS_get_bat(bus_src))
      BAT_set_bus(BUS_get_bat(bus_src),bus_dst);

    // Connections - csc conv
    while (BUS_get_csc_conv(bus_src))
      CONVCSC_set_ac_bus(BUS_get_csc_conv(bus_src),bus_dst);
    
    // Connections - vsc conv
    while (BUS_get_vsc_conv(bus_src))
      CONVVSC_set_ac_bus(BUS_get_vsc_conv(bus_src),bus_dst);

    // Connections - reg vsc conv
    while (BUS_get_reg_vsc_conv(bus_src))
      CONVVSC_set_reg_bus(BUS_get_reg_vsc_conv(bus_src),bus_dst);

    // Connection - facts k
    while (BUS_get_facts_k(bus_src))
      FACTS_set_bus_k(BUS_get_facts_k(bus_src),bus_dst);

    // Connection - facts m
    while (BUS_get_facts_m(bus_src))
      FACTS_set_bus_m(BUS_get_facts_m(bus_src),bus_dst);

    // Connection - reg facts
    while (BUS_get_reg_facts(bus_src))
      FACTS_set_reg_bus(BUS_get_reg_facts(bus_src),bus_dst);

    index++;
  }

  // Update hash
  NET_update_hash_tables(net);

  // Delete old buses
  BUS_array_del(old_bus_array,old_num_buses);
}

void NET_del_buses(Net* net, Bus** bus_ptr_array, int size) {
  /** Removes buses from the network. The entire bus array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Bus* bus_src;
  Bus* bus_dst;
  Bus* bus;
  Bus* old_bus_array;
  int old_num_buses;
  char* delete;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bus_ptr_array)
    return;

  // Old buses
  old_bus_array = net->bus;
  old_num_buses = net->num_buses;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_buses);
  for (i = 0; i < size; i++) {
    bus = bus_ptr_array[i];
    if (bus) {
      if (bus == NET_get_bus(net,BUS_get_index(bus))) { // bus present in the network
        if (!delete[BUS_get_index(bus)]) {
          delete[BUS_get_index(bus)] = 1;
          num++;
        }
      }
      else
        bus_ptr_array[i] = NULL;                        // bus not present in the network
    }
  }

  // New buses
  NET_set_bus_array(net,BUS_array_new(net->num_buses-num,net->num_periods),net->num_buses-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_buses; i++) {

    // Delete
    if (delete[i]) {

      bus = BUS_array_get(old_bus_array,i);

      // Clear connections
      while (BUS_get_gen(bus))
        GEN_set_bus(BUS_get_gen(bus),NULL);
      while (BUS_get_reg_gen(bus))
        GEN_set_reg_bus(BUS_get_reg_gen(bus),NULL);
      while (BUS_get_branch_k(bus))
        BRANCH_set_bus_k(BUS_get_branch_k(bus),NULL);
      while (BUS_get_branch_m(bus))
        BRANCH_set_bus_m(BUS_get_branch_m(bus),NULL);
      while (BUS_get_reg_tran(bus))
        BRANCH_set_reg_bus(BUS_get_reg_tran(bus),NULL);
      while (BUS_get_shunt(bus))
        SHUNT_set_bus(BUS_get_shunt(bus),NULL);
      while (BUS_get_reg_shunt(bus))
        SHUNT_set_reg_bus(BUS_get_reg_shunt(bus),NULL);
      while (BUS_get_load(bus))
        LOAD_set_bus(BUS_get_load(bus),NULL);
      while (BUS_get_vargen(bus))
        VARGEN_set_bus(BUS_get_vargen(bus),NULL);
      while (BUS_get_bat(bus))
        BAT_set_bus(BUS_get_bat(bus),NULL);
      while (BUS_get_csc_conv(bus))
        CONVCSC_set_ac_bus(BUS_get_csc_conv(bus),NULL);
      while (BUS_get_vsc_conv(bus))
        CONVVSC_set_ac_bus(BUS_get_vsc_conv(bus),NULL);
      while (BUS_get_reg_vsc_conv(bus))
        CONVVSC_set_reg_bus(BUS_get_reg_vsc_conv(bus),NULL);
      while (BUS_get_facts_k(bus))
        FACTS_set_bus_k(BUS_get_facts_k(bus),NULL);
      while (BUS_get_facts_m(bus))
        FACTS_set_bus_m(BUS_get_facts_m(bus),NULL);
      while (BUS_get_reg_facts(bus))
        FACTS_set_reg_bus(BUS_get_reg_facts(bus),NULL);
    }

    // Keep
    else {

      bus_src = BUS_array_get(old_bus_array,i);
      bus_dst = NET_get_bus(net,index);

      // Copy data (except index, hash info, and connections)
      BUS_copy_from_bus(bus_dst,bus_src);

      // Connections - gen
      while (BUS_get_gen(bus_src))
        GEN_set_bus(BUS_get_gen(bus_src),bus_dst);

      // Connections - reg gen
      while (BUS_get_reg_gen(bus_src))
        GEN_set_reg_bus(BUS_get_reg_gen(bus_src),bus_dst);
      
      // Connection - branch k
      while (BUS_get_branch_k(bus_src))
        BRANCH_set_bus_k(BUS_get_branch_k(bus_src),bus_dst);
      
      // Connection - branch m
      while (BUS_get_branch_m(bus_src))
        BRANCH_set_bus_m(BUS_get_branch_m(bus_src),bus_dst);
      
      // Connection - reg tran
      while (BUS_get_reg_tran(bus_src))
        BRANCH_set_reg_bus(BUS_get_reg_tran(bus_src),bus_dst);
      
      // Connections - shunt
      while (BUS_get_shunt(bus_src))
        SHUNT_set_bus(BUS_get_shunt(bus_src),bus_dst);
      
      // Connections - reg shunt
      while (BUS_get_reg_shunt(bus_src))
        SHUNT_set_reg_bus(BUS_get_reg_shunt(bus_src),bus_dst);
      
      // Connections - Load
      while (BUS_get_load(bus_src))
        LOAD_set_bus(BUS_get_load(bus_src),bus_dst);
      
      // Connections - vargen
      while (BUS_get_vargen(bus_src))
        VARGEN_set_bus(BUS_get_vargen(bus_src),bus_dst);
      
      // Connections - bat
      while (BUS_get_bat(bus_src))
        BAT_set_bus(BUS_get_bat(bus_src),bus_dst);

      // Connections - csc conv
      while (BUS_get_csc_conv(bus_src))
        CONVCSC_set_ac_bus(BUS_get_csc_conv(bus_src),bus_dst);

      // Connections - vsc conv
      while (BUS_get_vsc_conv(bus_src))
        CONVVSC_set_ac_bus(BUS_get_vsc_conv(bus_src),bus_dst);
      
      // Connections - reg vsc conv
      while (BUS_get_reg_vsc_conv(bus_src))
        CONVVSC_set_reg_bus(BUS_get_reg_vsc_conv(bus_src),bus_dst);

      // Connection - facts k
      while (BUS_get_facts_k(bus_src))
        FACTS_set_bus_k(BUS_get_facts_k(bus_src),bus_dst);
      
      // Connection - facts m
      while (BUS_get_facts_m(bus_src))
        FACTS_set_bus_m(BUS_get_facts_m(bus_src),bus_dst);
      
      // Connection - reg facts
      while (BUS_get_reg_facts(bus_src))
        FACTS_set_reg_bus(BUS_get_reg_facts(bus_src),bus_dst);
      
      index++;
    }
  }

  // Update hash
  NET_update_hash_tables(net);

  // Delete old buses
  BUS_array_del(old_bus_array,old_num_buses);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_branches(Net* net, Branch** br_ptr_array, int size) {
  /** Adds branches to the network. The entire branch array is
   * relocated, the data is copied (except flags for the new branches), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Branch* br_src;
  Branch* br_dst;
  Branch* br;
  Branch* old_br_array;
  int old_num_branches;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !br_ptr_array)
    return;

  // Old branches
  old_br_array = net->branch;
  old_num_branches = net->num_branches;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    br = br_ptr_array[i];
    if (br != NET_get_branch(net,BRANCH_get_index(br))) // not in the network
      num++;
    else
      br_ptr_array[i] = NULL;                      // clear to ignore below
  }
  
  // New branches
  NET_set_branch_array(net,BRANCH_array_new(net->num_branches+num,net->num_periods),net->num_branches+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_branches+size; i++) {

    if (i < old_num_branches)
      br_src = BRANCH_array_get(old_br_array,i);
    else
      br_src = br_ptr_array[i-old_num_branches];
    br_dst = NET_get_branch(net,index);

    // Check
    if (!br_src)
      continue;

    // Copy data
    BRANCH_copy_from_branch(br_dst,br_src);

    // Clear flags
    if (i >= old_num_branches) {
      BRANCH_clear_flags(br_dst,FLAG_VARS);
      BRANCH_clear_flags(br_dst,FLAG_FIXED);
      BRANCH_clear_flags(br_dst,FLAG_BOUNDED);
      BRANCH_clear_flags(br_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus_k = BRANCH_get_bus_k(br_src);
    bus_m = BRANCH_get_bus_m(br_src);
    reg_bus = BRANCH_get_reg_bus(br_src);

    // Clear connections bus - old branch 
    BRANCH_set_bus_k(br_src,NULL);      // also removes branch from bus->branches_k list
    BRANCH_set_bus_m(br_src,NULL);      // also removes branch from bus->branches_m list
    BRANCH_set_reg_bus(br_src,NULL);    // also removes branch from bus->reg_trans list

    // Add connections bus - new branch
    BRANCH_set_bus_k(br_dst,bus_k);      // also adds branch to bus->branches_k list
    BRANCH_set_bus_m(br_dst,bus_m);      // also adds branch to bus->branches_m list
    BRANCH_set_reg_bus(br_dst,reg_bus);  // also adds branch to bus->reg_trans list

    index++;
  }

  // Delete old branches
  BRANCH_array_del(old_br_array,old_num_branches);
}

void NET_del_branches(Net* net, Branch** br_ptr_array, int size) {
  /** Removes branches from the network. The entire branch array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Branch* br_src;
  Branch* br_dst;
  Branch* br;
  Branch* old_br_array;
  int old_num_branches;
  char* delete;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !br_ptr_array)
    return;

  // Old branches
  old_br_array = net->branch;
  old_num_branches = net->num_branches;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_branches);
  for (i = 0; i < size; i++) {
    br = br_ptr_array[i];
    if (br) {
      if (br == NET_get_branch(net,BRANCH_get_index(br))) { // branch to delete is present in network
        if (!delete[BRANCH_get_index(br)]) {
          delete[BRANCH_get_index(br)] = 1;
          num++;
        }
      }
      else
        br_ptr_array[i] = NULL;
    }
  }

  // New branches
  NET_set_branch_array(net,BRANCH_array_new(net->num_branches-num,net->num_periods),net->num_branches-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_branches; i++) {

    // Delete
    if (delete[i]) {

      br = BRANCH_array_get(old_br_array,i);

      // Clear connections bus - "branch to be deleted"
      BRANCH_set_bus_k(br,NULL);   // also removes branch from bus->branches_k list
      BRANCH_set_bus_m(br,NULL);   // also removes branch from bus->branches_m list
      BRANCH_set_reg_bus(br,NULL); // also removes branch from bus->reg_trans list
    }

    // Keep
    else {

      br_src = BRANCH_array_get(old_br_array,i);
      br_dst = NET_get_branch(net,index);

      // Copy data
      BRANCH_copy_from_branch(br_dst,br_src);

      // Save old connections
      bus_k = BRANCH_get_bus_k(br_src);
      bus_m = BRANCH_get_bus_m(br_src);
      reg_bus = BRANCH_get_reg_bus(br_src);
      
      // Clear connections bus - old branch 
      BRANCH_set_bus_k(br_src,NULL);   // also removes branch from bus->branches_k list
      BRANCH_set_bus_m(br_src,NULL);   // also removes branch from bus->branches_m list
      BRANCH_set_reg_bus(br_src,NULL); // also removes branch from bus->reg_trans list
      
      // Add connections bus - new branch
      BRANCH_set_bus_k(br_dst,bus_k);      // also adds branch to bus->branches_k list
      BRANCH_set_bus_m(br_dst,bus_m);      // also adds branch to bus->branches_m list
      BRANCH_set_reg_bus(br_dst,reg_bus);  // also adds branch to bus->reg_trans list
      
      index++;
    }
  }

  // Delete old branches
  BRANCH_array_del(old_br_array,old_num_branches);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_gens(Net* net, Gen** gen_ptr_array, int size) {
  /** Adds generators to the network. The entire generator array is
   * relocated, the data is copied (except flags for the new gens), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Gen* gen_src;
  Gen* gen_dst;
  Gen* gen;
  Gen* old_gen_array;
  int old_num_gens;
  Bus* bus;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !gen_ptr_array)
    return;

  // Old gens
  old_gen_array = net->gen;
  old_num_gens = net->num_gens;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    gen = gen_ptr_array[i];
    if (gen != NET_get_gen(net,GEN_get_index(gen))) // not in the network
      num++;
    else
      gen_ptr_array[i] = NULL;                      // clear to ignore below
  }    
  
  // New gens
  NET_set_gen_array(net,GEN_array_new(net->num_gens+num,net->num_periods),net->num_gens+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_gens+size; i++) {
    
    if (i < old_num_gens)
      gen_src = GEN_array_get(old_gen_array,i);
    else
      gen_src = gen_ptr_array[i-old_num_gens];
    gen_dst = NET_get_gen(net,index);

    // Check
    if (!gen_src)
      continue;

    // Copy data
    GEN_copy_from_gen(gen_dst,gen_src);

    // Clear flags
    if (i >= old_num_gens) {
      GEN_clear_flags(gen_dst,FLAG_VARS);
      GEN_clear_flags(gen_dst,FLAG_FIXED);
      GEN_clear_flags(gen_dst,FLAG_BOUNDED);
      GEN_clear_flags(gen_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = GEN_get_bus(gen_src);
    reg_bus = GEN_get_reg_bus(gen_src);

    // Clear connections bus - old gen 
    GEN_set_bus(gen_src,NULL);         // also removes gen from bus->gens list
    GEN_set_reg_bus(gen_src,NULL);     // also removes gen from bus->reg_gens list

    // Add connections bus - new gen
    GEN_set_bus(gen_dst,bus);          // also adds gen to bus->gens list
    GEN_set_reg_bus(gen_dst,reg_bus);  // also adds gen to bus->reg_gens list

    index++;
  }

  // Delete old gens
  GEN_array_del(old_gen_array,old_num_gens);
}

void NET_del_gens(Net* net, Gen** gen_ptr_array, int size) {
  /** Removes generators from the network. The entire generator array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Gen* gen_src;
  Gen* gen_dst;
  Gen* gen;
  Gen* old_gen_array;
  int old_num_gens;
  char* delete;
  Bus* bus;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !gen_ptr_array)
    return;

  // Old gens
  old_gen_array = net->gen;
  old_num_gens = net->num_gens;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_gens);
  for (i = 0; i < size; i++) {
    gen = gen_ptr_array[i];
    if (gen) {
      if (gen == NET_get_gen(net,GEN_get_index(gen))) { // gen to delete is present in network
        if (!delete[GEN_get_index(gen)]) {
          delete[GEN_get_index(gen)] = 1;
          num++;
        }
      }
      else
        gen_ptr_array[i] = NULL;
    }
  }

  // New gens
  NET_set_gen_array(net,GEN_array_new(net->num_gens-num,net->num_periods),net->num_gens-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_gens; i++) {

    // Delete
    if (delete[i]) {

      gen = GEN_array_get(old_gen_array,i);

      // Clear connections bus - "gen to be deleted"
      GEN_set_bus(gen,NULL);     // also removes gen from bus->gens list
      GEN_set_reg_bus(gen,NULL); // also removes gen from bus->reg_gens list
    }

    // Keep
    else {

      gen_src = GEN_array_get(old_gen_array,i);
      gen_dst = NET_get_gen(net,index);

      // Copy data
      GEN_copy_from_gen(gen_dst,gen_src);

      // Save old connections
      bus = GEN_get_bus(gen_src);        
      reg_bus = GEN_get_reg_bus(gen_src);
      
      // Clear connections bus - old gen 
      GEN_set_bus(gen_src,NULL);         // also removes gen from bus->gens list
      GEN_set_reg_bus(gen_src,NULL);     // also removes gen from bus->reg_gens list
      
      // Add connections bus - new gen
      GEN_set_bus(gen_dst,bus);          // also adds gen to bus->gens list
      GEN_set_reg_bus(gen_dst,reg_bus);  // also adds gen to bus->reg_gens list
      
      index++;
    }
  }

  // Delete old gens
  GEN_array_del(old_gen_array,old_num_gens);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_loads(Net* net, Load** load_ptr_array, int size) {
  /** Adds loads to the network. The entire load array is
   * relocated, the data is copied (except flags for the new loads), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Load* load_src;
  Load* load_dst;
  Load* load;
  Load* old_load_array;
  int old_num_loads;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !load_ptr_array)
    return;

  // Old loads
  old_load_array = net->load;
  old_num_loads = net->num_loads;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    load = load_ptr_array[i];
    if (load != NET_get_load(net,LOAD_get_index(load))) // not in the network
      num++;
    else
      load_ptr_array[i] = NULL;                      // clear to ignore below
  }

  // New loads
  NET_set_load_array(net,LOAD_array_new(net->num_loads+num,net->num_periods),net->num_loads+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_loads+size; i++) {

    if (i < old_num_loads)
      load_src = LOAD_array_get(old_load_array,i);
    else
      load_src = load_ptr_array[i-old_num_loads];
    load_dst = NET_get_load(net,index);

    // Check
    if (!load_src)
      continue;

    // Copy data
    LOAD_copy_from_load(load_dst,load_src);

    // Clear flags
    if (i >= old_num_loads) {
      LOAD_clear_flags(load_dst,FLAG_VARS);
      LOAD_clear_flags(load_dst,FLAG_FIXED);
      LOAD_clear_flags(load_dst,FLAG_BOUNDED);
      LOAD_clear_flags(load_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = LOAD_get_bus(load_src);

    // Clear connections bus - old load
    LOAD_set_bus(load_src,NULL);         // also removes load from bus->loads list

    // Add connections bus - new load
    LOAD_set_bus(load_dst,bus);          // also adds load to bus->loads list

    index++;
  }

  // Delete old loads
  LOAD_array_del(old_load_array,old_num_loads);
}

void NET_del_loads(Net* net, Load** load_ptr_array, int size) {
  /** Removes loads from the network. The entire load array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Load* load_src;
  Load* load_dst;
  Load* load;
  Load* old_load_array;
  int old_num_loads;
  char* delete;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !load_ptr_array)
    return;

  // Old loads
  old_load_array = net->load;
  old_num_loads = net->num_loads;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_loads);
  for (i = 0; i < size; i++) {
    load = load_ptr_array[i];
    if (load) {
      if (load == NET_get_load(net,LOAD_get_index(load))) { // load to delete is present in network
        if (!delete[LOAD_get_index(load)]) {
          delete[LOAD_get_index(load)] = 1;
          num++;
        }
      }
      else
        load_ptr_array[i] = NULL;
    }
  }

  // New loads
  NET_set_load_array(net,LOAD_array_new(net->num_loads-num,net->num_periods),net->num_loads-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_loads; i++) {

    // Delete
    if (delete[i]) {

      load = LOAD_array_get(old_load_array,i);
      
      // Clear connections bus - "load to be deleted"
      LOAD_set_bus(load,NULL); // also removes load from bus->loads list
    }

    // Keep
    else {

      load_src = LOAD_array_get(old_load_array,i);
      load_dst = NET_get_load(net,index);

      // Copy data
      LOAD_copy_from_load(load_dst,load_src);

      // Save old connections
      bus = LOAD_get_bus(load_src);
      
      // Clear connections bus - old load 
      LOAD_set_bus(load_src,NULL);       // also removes load from bus->loads list
      
      // Add connections bus - new load
      LOAD_set_bus(load_dst,bus);        // also adds load to bus->loads list
      
      index++;
    }
  }

  // Delete old loads
  LOAD_array_del(old_load_array,old_num_loads);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_shunts(Net* net, Shunt** shunt_ptr_array, int size) {
  /** Adds shunts to the network. The entire shunt array is
   * relocated, the data is copied (except flags for the new shunts), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Shunt* shunt_src;
  Shunt* shunt_dst;
  Shunt* shunt;
  Shunt* old_shunt_array;
  int old_num_shunts;
  Bus* bus;
  Bus* reg_bus;
  int index;
  int num;
  int i;

  // Check
  if (!net || !shunt_ptr_array)
    return;

  // Old shunts
  old_shunt_array = net->shunt;
  old_num_shunts = net->num_shunts;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    shunt = shunt_ptr_array[i];
    if (shunt != NET_get_shunt(net,SHUNT_get_index(shunt))) // not in the network
      num++;
    else
      shunt_ptr_array[i] = NULL;                      // clear to ignore below
  }

  // New shunts
  NET_set_shunt_array(net,SHUNT_array_new(net->num_shunts+num,net->num_periods),net->num_shunts+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_shunts+size; i++) {

    if (i < old_num_shunts)
      shunt_src = SHUNT_array_get(old_shunt_array,i);
    else
      shunt_src = shunt_ptr_array[i-old_num_shunts];
    shunt_dst = NET_get_shunt(net,index);

    // Check
    if (!shunt_src)
      continue;

    // Copy data
    SHUNT_copy_from_shunt(shunt_dst,shunt_src);

    // Clear flags
    if (i >= old_num_shunts) {
      SHUNT_clear_flags(shunt_dst,FLAG_VARS);
      SHUNT_clear_flags(shunt_dst,FLAG_FIXED);
      SHUNT_clear_flags(shunt_dst,FLAG_BOUNDED);
      SHUNT_clear_flags(shunt_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = SHUNT_get_bus(shunt_src);        
    reg_bus = SHUNT_get_reg_bus(shunt_src);

    // Clear connections bus - old shunt 
    SHUNT_set_bus(shunt_src,NULL);         // also removes shunt from bus->shunts list
    SHUNT_set_reg_bus(shunt_src,NULL);     // also removes shunt from bus->reg_shunts list

    // Add connections bus - new shunt
    SHUNT_set_bus(shunt_dst,bus);          // also adds shunt to bus->shunts list
    SHUNT_set_reg_bus(shunt_dst,reg_bus);  // also adds shunt to bus->reg_shunts list

    index++;
  }

  // Delete old shunts
  SHUNT_array_del(old_shunt_array,old_num_shunts);
}

void NET_del_shunts(Net* net, Shunt** shunt_ptr_array, int size) {
  /** Removes shunts from the network. The entire shunt array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Shunt* shunt_src;
  Shunt* shunt_dst;
  Shunt* shunt;
  Shunt* old_shunt_array;
  int old_num_shunts;
  char* delete;
  Bus* bus;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !shunt_ptr_array)
    return;

  // Old shunts
  old_shunt_array = net->shunt;
  old_num_shunts = net->num_shunts;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_shunts);
  for (i = 0; i < size; i++) {
    shunt = shunt_ptr_array[i];
    if (shunt) {
      if (shunt == NET_get_shunt(net,SHUNT_get_index(shunt))) { // shunt to delete is present in network
        if (!delete[SHUNT_get_index(shunt)]) {
          delete[SHUNT_get_index(shunt)] = 1;
          num++;
        }
      }
      else
        shunt_ptr_array[i] = NULL;
    }
  }

  // New shunts
  NET_set_shunt_array(net,SHUNT_array_new(net->num_shunts-num,net->num_periods),net->num_shunts-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_shunts; i++) {

    // Delete
    if (delete[i]) {

      shunt = SHUNT_array_get(old_shunt_array,i);
      
      // Clear connections bus - "shunt to be deleted"
      SHUNT_set_bus(shunt,NULL);     // also removes shunt from bus->shunts list
      SHUNT_set_reg_bus(shunt,NULL); // also removes shunt from bus->reg_shunts list
    }

    // Keep
    else {

      shunt_src = SHUNT_array_get(old_shunt_array,i);
      shunt_dst = NET_get_shunt(net,index);

      // Copy data
      SHUNT_copy_from_shunt(shunt_dst,shunt_src);
      
      // Save old connections
      bus = SHUNT_get_bus(shunt_src);        
      reg_bus = SHUNT_get_reg_bus(shunt_src);
      
      // Clear connections bus - old shunt 
      SHUNT_set_bus(shunt_src,NULL);         // also removes shunt from bus->shunts list
      SHUNT_set_reg_bus(shunt_src,NULL);     // also removes shunt from bus->reg_shunts list
      
      // Add connections bus - new shunt
      SHUNT_set_bus(shunt_dst,bus);          // also adds shunt to bus->shunts list
      SHUNT_set_reg_bus(shunt_dst,reg_bus);  // also adds shunt to bus->reg_shunts list
      
      index++;
    }
  }

  // Delete old shunts
  SHUNT_array_del(old_shunt_array,old_num_shunts);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_bats(Net* net, Bat** bat_ptr_array, int size) {
  /** Adds batteries to the network. The entire battery array is
   * relocated, the data is copied (except flags for the new batteries), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Bat* bat_src;
  Bat* bat_dst;
  Bat* bat;
  Bat* old_bat_array;
  int old_num_bats;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bat_ptr_array)
    return;

  // Old batteries
  old_bat_array = net->bat;
  old_num_bats = net->num_bats;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    bat = bat_ptr_array[i];
    if (bat != NET_get_bat(net,BAT_get_index(bat))) // not in the network
      num++;
    else
      bat_ptr_array[i] = NULL;                      // clear to ignore below
  }

  // New batteries
  NET_set_bat_array(net,BAT_array_new(net->num_bats+num,net->num_periods),net->num_bats+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_bats+size; i++) {

    if (i < old_num_bats)
      bat_src = BAT_array_get(old_bat_array,i);
    else
      bat_src = bat_ptr_array[i-old_num_bats];
    bat_dst = NET_get_bat(net,index);

    // Check
    if (!bat_src)
      continue;

    // Copy data
    BAT_copy_from_bat(bat_dst,bat_src);

    // Clear flags
    if (i >= old_num_bats) {
      BAT_clear_flags(bat_dst,FLAG_VARS);
      BAT_clear_flags(bat_dst,FLAG_FIXED);
      BAT_clear_flags(bat_dst,FLAG_BOUNDED);
      BAT_clear_flags(bat_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = BAT_get_bus(bat_src);

    // Clear connections bus - old bat
    BAT_set_bus(bat_src,NULL);         // also removes bat from bus->bats list

    // Add connections bus - new bat
    BAT_set_bus(bat_dst,bus);          // also adds bat to bus->bats list

    index++;
  }

  // Delete old bats
  BAT_array_del(old_bat_array,old_num_bats);
}

void NET_del_bats(Net* net, Bat** bat_ptr_array, int size) {
  /** Removes batteries from the network. The entire battery array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Bat* bat_src;
  Bat* bat_dst;
  Bat* bat;
  Bat* old_bat_array;
  int old_num_bats;
  char* delete;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bat_ptr_array)
    return;

  // Old batteries
  old_bat_array = net->bat;
  old_num_bats = net->num_bats;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_bats);
  for (i = 0; i < size; i++) {
    bat = bat_ptr_array[i];
    if (bat) {
      if (bat == NET_get_bat(net,BAT_get_index(bat))) {
        if (!delete[BAT_get_index(bat)]) {
          delete[BAT_get_index(bat)] = 1;
          num++;
        }
      }
      else
        bat_ptr_array[i] = NULL;
    }
  }

  // New batteries
  NET_set_bat_array(net,BAT_array_new(net->num_bats-num,net->num_periods),net->num_bats-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_bats; i++) {

    // Delete
    if (delete[i]) {

      bat = BAT_array_get(old_bat_array,i);
      
      // Clear connections bus - "bat to be deleted"
      BAT_set_bus(bat,NULL); // also removes bat from bus->bats list
    }

    // Keep
    else {

      bat_src = BAT_array_get(old_bat_array,i);
      bat_dst = NET_get_bat(net,index);

      // Copy data
      BAT_copy_from_bat(bat_dst,bat_src);

      // Save old connections
      bus = BAT_get_bus(bat_src);
      
      // Clear connections bus - old bat 
      BAT_set_bus(bat_src,NULL);       // also removes bat from bus->bats list
      
      // Add connections bus - new bat
      BAT_set_bus(bat_dst,bus);        // also adds bat to bus->bats list
      
      index++;
    }
  }

  // Delete old batteries
  BAT_array_del(old_bat_array,old_num_bats);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_vargens(Net* net, Vargen** gen_ptr_array, int size) {
  /** Adds var generators to the network. The entire var generator array is
   * relocated, the data is copied (except flags for the new var generators), 
   * and the bus connections are stolen.
   */
  
  // Local variables
  Vargen* gen_src;
  Vargen* gen_dst;
  Vargen* gen;
  Vargen* old_gen_array;
  int old_num_vargens;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !gen_ptr_array)
    return;

  // Old var generators
  old_gen_array = net->vargen;
  old_num_vargens = net->num_vargens;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    gen = gen_ptr_array[i];
    if (gen != NET_get_vargen(net,VARGEN_get_index(gen))) // not in the network
      num++;
    else
      gen_ptr_array[i] = NULL;                      // clear to ignore below
  }

  // New var generators
  NET_set_vargen_array(net,VARGEN_array_new(net->num_vargens+num,net->num_periods),net->num_vargens+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_vargens+size; i++) {

    if (i < old_num_vargens)
      gen_src = VARGEN_array_get(old_gen_array,i);
    else
      gen_src = gen_ptr_array[i-old_num_vargens];
    gen_dst = NET_get_vargen(net,index);

    // Check
    if (!gen_src)
      continue;

    // Copy data
    VARGEN_copy_from_vargen(gen_dst,gen_src);

    // Clear flags
    if (i >= old_num_vargens) {
      VARGEN_clear_flags(gen_dst,FLAG_VARS);
      VARGEN_clear_flags(gen_dst,FLAG_FIXED);
      VARGEN_clear_flags(gen_dst,FLAG_BOUNDED);
      VARGEN_clear_flags(gen_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = VARGEN_get_bus(gen_src);

    // Clear connections bus - old vargen
    VARGEN_set_bus(gen_src,NULL);         // also removes vargen from bus->vargens list

    // Add connections bus - new vargens
    VARGEN_set_bus(gen_dst,bus);          // also adds vargen to bus->vargens list

    index++;
  }

  // Delete old vargens
  VARGEN_array_del(old_gen_array,old_num_vargens);
}

void NET_del_vargens(Net* net, Vargen** gen_ptr_array, int size) {
  /** Removes var generators from the network. The entire var generator array is
   * relocated, the data is copied, and the bus connections are set.
   * Network flags are cleared. 
   */
  
  // Local variables
  Vargen* gen_src;
  Vargen* gen_dst;
  Vargen* gen;
  Vargen* old_gen_array;
  int old_num_vargens;
  char* delete;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !gen_ptr_array)
    return;

  // Old var generators
  old_gen_array = net->vargen;
  old_num_vargens = net->num_vargens;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_vargens);
  for (i = 0; i < size; i++) {
    gen = gen_ptr_array[i];
    if (gen) {
      if (gen == NET_get_vargen(net,VARGEN_get_index(gen))) {
        if (!delete[VARGEN_get_index(gen)]) {
          delete[VARGEN_get_index(gen)] = 1;
          num++;
        }
      }
      else
        gen_ptr_array[i] = NULL;
    }
  }

  // New var generators
  NET_set_vargen_array(net,VARGEN_array_new(net->num_vargens-num,net->num_periods),net->num_vargens-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_vargens; i++) {

    // Delete
    if (delete[i]) {

      gen = VARGEN_array_get(old_gen_array,i);
      
      // Clear connections bus - "vargen to be deleted"
      VARGEN_set_bus(gen,NULL); // also removes vargen from bus->vargens list
    }

    // Keep
    else {

      gen_src = VARGEN_array_get(old_gen_array,i);
      gen_dst = NET_get_vargen(net,index);

      // Copy data
      VARGEN_copy_from_vargen(gen_dst,gen_src);

      // Save old connections
      bus = VARGEN_get_bus(gen_src);
      
      // Clear connections bus - old vargen 
      VARGEN_set_bus(gen_src,NULL);       // also removes vargen from bus->vargens list
      
      // Add connections bus - new vargen
      VARGEN_set_bus(gen_dst,bus);        // also adds vargen to bus->vargens list
      
      index++;
    }
  }

  // Delete old var generators
  VARGEN_array_del(old_gen_array,old_num_vargens);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_csc_convs(Net* net, ConvCSC** conv_ptr_array, int size) {
  /** Adds CSC converters to the network. The entire converter array is
   *  relocated, the data is copied (except flags for the new converters), 
   *  and the bus connections are stolen.
   */
  
  // Local variables
  ConvCSC* conv_src;
  ConvCSC* conv_dst;
  ConvCSC* conv;
  ConvCSC* old_conv_array;
  int old_num_convs;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !conv_ptr_array)
    return;

  // Old converters
  old_conv_array = net->csc_conv;
  old_num_convs = net->num_csc_convs;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    conv = conv_ptr_array[i];
    if (conv != NET_get_csc_conv(net,CONVCSC_get_index(conv))) // not in the network
      num++;
    else
      conv_ptr_array[i] = NULL;                         // clear to ignore below
  }

  // New converters
  NET_set_csc_conv_array(net,CONVCSC_array_new(net->num_csc_convs+num,net->num_periods),net->num_csc_convs+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_convs+size; i++) {

    if (i < old_num_convs)
      conv_src = CONVCSC_array_get(old_conv_array,i);
    else
      conv_src = conv_ptr_array[i-old_num_convs];
    conv_dst = NET_get_csc_conv(net,index);

    // Check
    if (!conv_src)
      continue;

    // Copy data
    CONVCSC_copy_from_conv(conv_dst,conv_src);

    // Clear flags
    if (i >= old_num_convs) {
      CONVCSC_clear_flags(conv_dst,FLAG_VARS);
      CONVCSC_clear_flags(conv_dst,FLAG_FIXED);
      CONVCSC_clear_flags(conv_dst,FLAG_BOUNDED);
      CONVCSC_clear_flags(conv_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = CONVCSC_get_ac_bus(conv_src);

    // Clear connections bus - old conv
    CONVCSC_set_ac_bus(conv_src,NULL);     // also removes conv from bus->convs list

    // Add connections bus - new conv
    CONVCSC_set_ac_bus(conv_dst,bus);      // also adds conv to bus->convs list

    index++;
  }

  // Delete old convs
  CONVCSC_array_del(old_conv_array,old_num_convs);
}

void NET_del_csc_convs(Net* net, ConvCSC** conv_ptr_array, int size) {
  /** Removes CSC converters from the network. The entire converter array is
   *  relocated, the data is copied, and the bus connections are set.
   *  Network flags are cleared. 
   */
  
  // Local variables
  ConvCSC* conv_src;
  ConvCSC* conv_dst;
  ConvCSC* conv;
  ConvCSC* old_conv_array;
  int old_num_convs;
  char* delete;
  Bus* bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !conv_ptr_array)
    return;

  // Old converters
  old_conv_array = net->csc_conv;
  old_num_convs = net->num_csc_convs;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_csc_convs);
  for (i = 0; i < size; i++) {
    conv = conv_ptr_array[i];
    if (conv) {
      if (conv == NET_get_csc_conv(net,CONVCSC_get_index(conv))) {
        if (!delete[CONVCSC_get_index(conv)]) {
          delete[CONVCSC_get_index(conv)] = 1;
          num++;
        }
      }
      else
        conv_ptr_array[i] = NULL;
    }
  }

  // New converters
  NET_set_csc_conv_array(net,CONVCSC_array_new(net->num_csc_convs-num,net->num_periods),net->num_csc_convs-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_convs; i++) {

    // Delete
    if (delete[i]) {

      conv = CONVCSC_array_get(old_conv_array,i);
      
      // Clear connections bus - "conv to be deleted"
      CONVCSC_set_ac_bus(conv,NULL); // also removes conv from bus->convs list
    }

    // Keep
    else {

      conv_src = CONVCSC_array_get(old_conv_array,i);
      conv_dst = NET_get_csc_conv(net,index);

      // Copy data
      CONVCSC_copy_from_conv(conv_dst,conv_src);

      // Save old connections
      bus = CONVCSC_get_ac_bus(conv_src);
      
      // Clear connections bus - old conv 
      CONVCSC_set_ac_bus(conv_src,NULL);       // also removes conv from bus->convs list
      
      // Add connections bus - new conv
      CONVCSC_set_ac_bus(conv_dst,bus);        // also adds conv to bus->convs list
      
      index++;
    }
  }

  // Delete old converters
  CONVCSC_array_del(old_conv_array,old_num_convs);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_vsc_convs(Net* net, ConvVSC** conv_ptr_array, int size) {
  /** Adds VSC converters to the network. The entire converter array is
   *  relocated, the data is copied (except flags for the new converters), 
   *  and the bus connections are stolen.
   */
  
  // Local variables
  ConvVSC* conv_src;
  ConvVSC* conv_dst;
  ConvVSC* conv;
  ConvVSC* old_conv_array;
  int old_num_convs;
  Bus* bus;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !conv_ptr_array)
    return;

  // Old converters
  old_conv_array = net->vsc_conv;
  old_num_convs = net->num_vsc_convs;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    conv = conv_ptr_array[i];
    if (conv != NET_get_vsc_conv(net,CONVVSC_get_index(conv))) // not in the network
      num++;
    else
      conv_ptr_array[i] = NULL;                         // clear to ignore below
  }

  // New converters
  NET_set_vsc_conv_array(net,CONVVSC_array_new(net->num_vsc_convs+num,net->num_periods),net->num_vsc_convs+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_convs+size; i++) {

    if (i < old_num_convs)
      conv_src = CONVVSC_array_get(old_conv_array,i);
    else
      conv_src = conv_ptr_array[i-old_num_convs];
    conv_dst = NET_get_vsc_conv(net,index);

    // Check
    if (!conv_src)
      continue;

    // Copy data
    CONVVSC_copy_from_conv(conv_dst,conv_src);

    // Clear flags
    if (i >= old_num_convs) {
      CONVVSC_clear_flags(conv_dst,FLAG_VARS);
      CONVVSC_clear_flags(conv_dst,FLAG_FIXED);
      CONVVSC_clear_flags(conv_dst,FLAG_BOUNDED);
      CONVVSC_clear_flags(conv_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus = CONVVSC_get_ac_bus(conv_src);
    reg_bus = CONVVSC_get_reg_bus(conv_src);

    // Clear connections bus - old conv
    CONVVSC_set_ac_bus(conv_src,NULL);     // also removes conv from bus->convs list
    CONVVSC_set_reg_bus(conv_src,NULL);

    // Add connections bus - new conv
    CONVVSC_set_ac_bus(conv_dst,bus);      // also adds conv to bus->convs list
    CONVVSC_set_reg_bus(conv_dst,reg_bus);

    index++;
  }

  // Delete old convs
  CONVVSC_array_del(old_conv_array,old_num_convs);
}

void NET_del_vsc_convs(Net* net, ConvVSC** conv_ptr_array, int size) {
  /** Removes VSC converters from the network. The entire converter array is
   *  relocated, the data is copied, and the bus connections are set.
   *  Network flags are cleared. 
   */
  
  // Local variables
  ConvVSC* conv_src;
  ConvVSC* conv_dst;
  ConvVSC* conv;
  ConvVSC* old_conv_array;
  int old_num_convs;
  char* delete;
  Bus* bus;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !conv_ptr_array)
    return;

  // Old converters
  old_conv_array = net->vsc_conv;
  old_num_convs = net->num_vsc_convs;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_vsc_convs);
  for (i = 0; i < size; i++) {
    conv = conv_ptr_array[i];
    if (conv) {
      if (conv == NET_get_vsc_conv(net,CONVVSC_get_index(conv))) {
        if (!delete[CONVVSC_get_index(conv)]) {
          delete[CONVVSC_get_index(conv)] = 1;
          num++;
        }
      }
      else
        conv_ptr_array[i] = NULL;
    }
  }

  // New converters
  NET_set_vsc_conv_array(net,CONVVSC_array_new(net->num_vsc_convs-num,net->num_periods),net->num_vsc_convs-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_convs; i++) {

    // Delete
    if (delete[i]) {

      conv = CONVVSC_array_get(old_conv_array,i);
      
      // Clear connections bus - "conv to be deleted"
      CONVVSC_set_ac_bus(conv,NULL);  // also removes conv from bus->convs list
      CONVVSC_set_reg_bus(conv,NULL);
    }

    // Keep
    else {

      conv_src = CONVVSC_array_get(old_conv_array,i);
      conv_dst = NET_get_vsc_conv(net,index);

      // Copy data
      CONVVSC_copy_from_conv(conv_dst,conv_src);

      // Save old connections
      bus = CONVVSC_get_ac_bus(conv_src);
      reg_bus = CONVVSC_get_reg_bus(conv_src);
      
      // Clear connections bus - old conv 
      CONVVSC_set_ac_bus(conv_src,NULL);       // also removes conv from bus->convs list
      CONVVSC_set_reg_bus(conv_src,NULL);
      
      // Add connections bus - new conv
      CONVVSC_set_ac_bus(conv_dst,bus);        // also adds conv to bus->convs list
      CONVVSC_set_reg_bus(conv_dst,reg_bus);
      
      index++;
    }
  }

  // Delete old converters
  CONVVSC_array_del(old_conv_array,old_num_convs);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_dc_buses(Net* net, BusDC** bus_ptr_array, int size) {
  /** Adds dc buses to the network. The entire bus array is
   *  relocated, the data is copied (except flags for the new buses), 
   *  and the bus connections are stolen.
   */
  
  // Local variables
  BusDC* bus_src;
  BusDC* bus_dst;
  BusDC* bus;
  BusDC* old_bus_array; 
  int old_num_buses;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bus_ptr_array)
    return;

  // Old buses
  old_bus_array = net->dc_bus;
  old_num_buses = net->num_dc_buses;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    bus = bus_ptr_array[i];
    if (bus != NET_get_dc_bus(net,BUSDC_get_index(bus))) // not in the network
      num++;
    else
      bus_ptr_array[i] = NULL;                           // clear to ignore below
  }

  // New buses
  NET_set_dc_bus_array(net,BUSDC_array_new(net->num_dc_buses+num,net->num_periods),net->num_dc_buses+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_buses+size; i++) {

    if (i < old_num_buses)
      bus_src = BUSDC_array_get(old_bus_array,i);
    else
      bus_src = bus_ptr_array[i-old_num_buses];
    bus_dst = NET_get_dc_bus(net,index);

    // Check
    if (!bus_src)
      continue;

    // Copy data (except index, hash info, and connections)
    BUSDC_copy_from_dc_bus(bus_dst,bus_src);
    
    // Clear flags
    if (i >= old_num_buses)
      BUSDC_clear_flags(bus_dst,FLAG_ALL);

    // Connection - branch k
    while (BUSDC_get_branch_k(bus_src))
      BRANCHDC_set_bus_k(BUSDC_get_branch_k(bus_src),bus_dst);

    // Connection - branch m
    while (BUSDC_get_branch_m(bus_src))
      BRANCHDC_set_bus_m(BUSDC_get_branch_m(bus_src),bus_dst);

    // Connections - csc conv
    while (BUSDC_get_csc_conv(bus_src))
      CONVCSC_set_dc_bus(BUSDC_get_csc_conv(bus_src),bus_dst);
    
    // Connections - vsc conv
    while (BUSDC_get_vsc_conv(bus_src))
      CONVVSC_set_dc_bus(BUSDC_get_vsc_conv(bus_src),bus_dst);

    index++;
  }

  // Update hash
  NET_update_hash_tables(net);

  // Delete old buses
  BUSDC_array_del(old_bus_array,old_num_buses);
}

void NET_del_dc_buses(Net* net, BusDC** bus_ptr_array, int size) {
  /** Removes dc buses from the network. The entire bus array is
   *  relocated, the data is copied, and the bus connections are set.
   *  Network flags are cleared. 
   */
  
  // Local variables
  BusDC* bus_src;
  BusDC* bus_dst;
  BusDC* bus;
  BusDC* old_bus_array;
  int old_num_buses;
  char* delete;
  int num;
  int index;
  int i;

  // Check
  if (!net || !bus_ptr_array)
    return;

  // Old buses
  old_bus_array = net->dc_bus;
  old_num_buses = net->num_dc_buses;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_dc_buses);
  for (i = 0; i < size; i++) {
    bus = bus_ptr_array[i];
    if (bus) {
      if (bus == NET_get_dc_bus(net,BUSDC_get_index(bus))) { // bus present in the network
        if (!delete[BUSDC_get_index(bus)]) {
          delete[BUSDC_get_index(bus)] = 1;
          num++;
        }
      }
      else
        bus_ptr_array[i] = NULL;                        // bus not present in the network
    }
  }

  // New buses
  NET_set_dc_bus_array(net,BUSDC_array_new(net->num_dc_buses-num,net->num_periods),net->num_dc_buses-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_buses; i++) {

    // Delete
    if (delete[i]) {

      bus = BUSDC_array_get(old_bus_array,i);
      
      // Clear connections
      while (BUSDC_get_branch_k(bus))
        BRANCHDC_set_bus_k(BUSDC_get_branch_k(bus),NULL);
      while (BUSDC_get_branch_m(bus))
        BRANCHDC_set_bus_m(BUSDC_get_branch_m(bus),NULL);
      while (BUSDC_get_csc_conv(bus))
        CONVCSC_set_dc_bus(BUSDC_get_csc_conv(bus),NULL);
      while (BUSDC_get_vsc_conv(bus))
        CONVVSC_set_dc_bus(BUSDC_get_vsc_conv(bus),NULL);
    }
    
    // Keep
    else {
      
      bus_src = BUSDC_array_get(old_bus_array,i);
      bus_dst = NET_get_dc_bus(net,index);

      // Copy data (except index, hash info, and connections)
      BUSDC_copy_from_dc_bus(bus_dst,bus_src);
      
      // Connection - branch k
      while (BUSDC_get_branch_k(bus_src))
        BRANCHDC_set_bus_k(BUSDC_get_branch_k(bus_src),bus_dst);
      
      // Connection - branch m
      while (BUSDC_get_branch_m(bus_src))
        BRANCHDC_set_bus_m(BUSDC_get_branch_m(bus_src),bus_dst);

      // Connections - csc conv
      while (BUSDC_get_csc_conv(bus_src))
        CONVCSC_set_dc_bus(BUSDC_get_csc_conv(bus_src),bus_dst);
            
      // Connections - vsc conv
      while (BUSDC_get_vsc_conv(bus_src))
        CONVVSC_set_dc_bus(BUSDC_get_vsc_conv(bus_src),bus_dst);
      
      index++;
    }
  }

  // Update hash
  NET_update_hash_tables(net);

  // Delete old buses
  BUSDC_array_del(old_bus_array,old_num_buses);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_dc_branches(Net* net, BranchDC** br_ptr_array, int size) {
  /** Adds DC branches to the network. The entire branch array is
   *  relocated, the data is copied (except flags for the new branches), 
   *  and the bus connections are stolen.
   */
  
  // Local variables
  BranchDC* br_src;
  BranchDC* br_dst;
  BranchDC* br;
  BranchDC* old_br_array;
  int old_num_branches;
  BusDC* bus_k;
  BusDC* bus_m;
  int num;
  int index;
  int i;

  // Check
  if (!net || !br_ptr_array)
    return;

  // Old branches
  old_br_array = net->dc_branch;
  old_num_branches = net->num_dc_branches;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    br = br_ptr_array[i];
    if (br != NET_get_dc_branch(net,BRANCHDC_get_index(br))) // not in the network
      num++;
    else
      br_ptr_array[i] = NULL;                      // clear to ignore below
  }
  
  // New branches
  NET_set_dc_branch_array(net,BRANCHDC_array_new(net->num_dc_branches+num,net->num_periods),net->num_dc_branches+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_branches+size; i++) {

    if (i < old_num_branches)
      br_src = BRANCHDC_array_get(old_br_array,i);
    else
      br_src = br_ptr_array[i-old_num_branches];
    br_dst = NET_get_dc_branch(net,index);

    // Check
    if (!br_src)
      continue;

    // Copy data
    BRANCHDC_copy_from_dc_branch(br_dst,br_src);

    // Clear flags
    if (i >= old_num_branches)
      BRANCHDC_clear_flags(br_dst,FLAG_ALL);
    
    // Save old connections
    bus_k = BRANCHDC_get_bus_k(br_src);
    bus_m = BRANCHDC_get_bus_m(br_src);

    // Clear connections bus - old branch 
    BRANCHDC_set_bus_k(br_src,NULL);      // also removes branch from bus->branches_k list
    BRANCHDC_set_bus_m(br_src,NULL);      // also removes branch from bus->branches_m list

    // Add connections bus - new branch
    BRANCHDC_set_bus_k(br_dst,bus_k);      // also adds branch to bus->branches_k list
    BRANCHDC_set_bus_m(br_dst,bus_m);      // also adds branch to bus->branches_m list

    index++;
  }

  // Delete old branches
  BRANCHDC_array_del(old_br_array,old_num_branches);
}

void NET_del_dc_branches(Net* net, BranchDC** br_ptr_array, int size) {
  /** Removes DC branches from the network. The entire branch array is
   *  relocated, the data is copied, and the bus connections are set.
   *  Network flags are cleared. 
   */
  
  // Local variables
  BranchDC* br_src;
  BranchDC* br_dst;
  BranchDC* br;
  BranchDC* old_br_array;
  int old_num_branches;
  char* delete;
  BusDC* bus_k;
  BusDC* bus_m;
  int num;
  int index;
  int i;

  // Check
  if (!net || !br_ptr_array)
    return;

  // Old branches
  old_br_array = net->dc_branch;
  old_num_branches = net->num_dc_branches;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_dc_branches);
  for (i = 0; i < size; i++) {
    br = br_ptr_array[i];
    if (br) {
      if (br == NET_get_dc_branch(net,BRANCHDC_get_index(br))) { // branch to delete is present in network
        if (!delete[BRANCHDC_get_index(br)]) {
          delete[BRANCHDC_get_index(br)] = 1;
          num++;
        }
      }
      else
        br_ptr_array[i] = NULL;
    }
  }

  // New branches
  NET_set_dc_branch_array(net,BRANCHDC_array_new(net->num_dc_branches-num,net->num_periods),net->num_dc_branches-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_branches; i++) {

    // Delete
    if (delete[i]) {

      br = BRANCHDC_array_get(old_br_array,i);

      // Clear connections bus - "branch to be deleted"
      BRANCHDC_set_bus_k(br,NULL);   // also removes branch from bus->branches_k list
      BRANCHDC_set_bus_m(br,NULL);   // also removes branch from bus->branches_m list
    }

    // Keep
    else {

      br_src = BRANCHDC_array_get(old_br_array,i);
      br_dst = NET_get_dc_branch(net,index);

      // Copy data
      BRANCHDC_copy_from_dc_branch(br_dst,br_src);

      // Save old connections
      bus_k = BRANCHDC_get_bus_k(br_src);
      bus_m = BRANCHDC_get_bus_m(br_src);
      
      // Clear connections bus - old branch 
      BRANCHDC_set_bus_k(br_src,NULL);   // also removes branch from bus->branches_k list
      BRANCHDC_set_bus_m(br_src,NULL);   // also removes branch from bus->branches_m list
      
      // Add connections bus - new branch
      BRANCHDC_set_bus_k(br_dst,bus_k);      // also adds branch to bus->branches_k list
      BRANCHDC_set_bus_m(br_dst,bus_m);      // also adds branch to bus->branches_m list
      
      index++;
    }
  }

  // Delete old branches
  BRANCHDC_array_del(old_br_array,old_num_branches);

  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_facts(Net* net, Facts** f_ptr_array, int size) {
  /** Adds facts devices to the network. The entire facts array is
   *  relocated, the data is copied (except flags for the new facts), 
   *  and the bus connections are stolen.
   */
  
  // Local variables
  Facts* f_src;
  Facts* f_dst;
  Facts* f;
  Facts* old_f_array;
  int old_num_facts;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !f_ptr_array)
    return;

  // Old facts
  old_f_array = net->facts;
  old_num_facts = net->num_facts;

  // Count
  num = 0;
  for (i = 0; i < size; i++) {
    f = f_ptr_array[i];
    if (f != NET_get_facts(net,FACTS_get_index(f))) // not in the network
      num++;
    else
      f_ptr_array[i] = NULL;                        // clear to ignore below
  }
  
  // New facts
  NET_set_facts_array(net,FACTS_array_new(net->num_facts+num,net->num_periods),net->num_facts+num);

  // Copy data and steal connections
  index = 0;
  for (i = 0; i < old_num_facts+size; i++) {

    if (i < old_num_facts)
      f_src = FACTS_array_get(old_f_array,i);
    else
      f_src = f_ptr_array[i-old_num_facts];
    f_dst = NET_get_facts(net,index);

    // Check
    if (!f_src)
      continue;

    // Copy data
    FACTS_copy_from_facts(f_dst,f_src);

    // Clear flags
    if (i >= old_num_facts) {
      FACTS_clear_flags(f_dst,FLAG_VARS);
      FACTS_clear_flags(f_dst,FLAG_FIXED);
      FACTS_clear_flags(f_dst,FLAG_BOUNDED);
      FACTS_clear_flags(f_dst,FLAG_SPARSE);
    }
    
    // Save old connections
    bus_k = FACTS_get_bus_k(f_src);
    bus_m = FACTS_get_bus_m(f_src);
    reg_bus = FACTS_get_reg_bus(f_src);

    // Clear connections bus - old facts
    FACTS_set_bus_k(f_src,NULL);      // also removes facts from bus->facts_k list
    FACTS_set_bus_m(f_src,NULL);      // also removes facts from bus->facts_m list
    FACTS_set_reg_bus(f_src,NULL);    // also removes facts from bus->reg_facts list

    // Add connections bus - new facts
    FACTS_set_bus_k(f_dst,bus_k);      // also adds facts to bus->facts_k list
    FACTS_set_bus_m(f_dst,bus_m);      // also adds facts to bus->facts_m list
    FACTS_set_reg_bus(f_dst,reg_bus);  // also adds facts to bus->reg_facts list

    index++;
  }

  // Delete old facts
  FACTS_array_del(old_f_array,old_num_facts);
}

void NET_del_facts(Net* net, Facts** f_ptr_array, int size) {
  /** Removes facts from the network. The entire facts array is
   *  relocated, the data is copied, and the bus connections are set.
   *  Network flags are cleared. 
   */
  
  // Local variables
  Facts* f_src;
  Facts* f_dst;
  Facts* f;
  Facts* old_f_array;
  int old_num_facts;
  char* delete;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;
  int num;
  int index;
  int i;

  // Check
  if (!net || !f_ptr_array)
    return;

  // Old facts
  old_f_array = net->facts;
  old_num_facts = net->num_facts;

  // Count unique and mark for deletion
  num = 0;
  ARRAY_zalloc(delete,char,net->num_facts);
  for (i = 0; i < size; i++) {
    f = f_ptr_array[i];
    if (f) {
      if (f == NET_get_facts(net,FACTS_get_index(f))) { // facts to delete is present in network
        if (!delete[FACTS_get_index(f)]) {
          delete[FACTS_get_index(f)] = 1;
          num++;
        }
      }
      else
        f_ptr_array[i] = NULL;
    }
  }
  
  // New facts
  NET_set_facts_array(net,FACTS_array_new(net->num_facts-num,net->num_periods),net->num_facts-num);

  // Copy data and set connections
  index = 0;
  for (i = 0; i < old_num_facts; i++) {

    // Delete
    if (delete[i]) {

      f = FACTS_array_get(old_f_array,i);

      // Clear connections bus - "facts to be deleted"
      FACTS_set_bus_k(f,NULL);   // also removes facts from bus->facts_k list
      FACTS_set_bus_m(f,NULL);   // also removes facts from bus->facts_m list
      FACTS_set_reg_bus(f,NULL); // also removes facts from bus->reg_facts list
    }

    // Keep
    else {

      f_src = FACTS_array_get(old_f_array,i);
      f_dst = NET_get_facts(net,index);

      // Copy data
      FACTS_copy_from_facts(f_dst,f_src);

      // Save old connections
      bus_k = FACTS_get_bus_k(f_src);
      bus_m = FACTS_get_bus_m(f_src);
      reg_bus = FACTS_get_reg_bus(f_src);
      
      // Clear connections bus - old facts
      FACTS_set_bus_k(f_src,NULL);   // also removes facts from bus->facts_k list
      FACTS_set_bus_m(f_src,NULL);   // also removes facts from bus->facts_m list
      FACTS_set_reg_bus(f_src,NULL); // also removes facts from bus->reg_facts list
      
      // Add connections bus - new facts
      FACTS_set_bus_k(f_dst,bus_k);      // also adds facts to bus->facts_k list
      FACTS_set_bus_m(f_dst,bus_m);      // also adds facts to bus->facts_m list
      FACTS_set_reg_bus(f_dst,reg_bus);  // also adds facts to bus->reg_facts list
      
      index++;
    }
  }

  // Delete old facts
  FACTS_array_del(old_f_array,old_num_facts);
  
  // Delete delete flags
  free(delete);

  // Clear flags
  NET_clear_flags(net);
}

void NET_add_vargens_from_params(Net* net, Bus* bus_list, REAL power_capacity, REAL power_base, REAL power_std, REAL corr_radius, REAL corr_value) {
  
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
  VARGEN_array_del(net->vargen,net->num_vargens);
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
  VARGEN_array_del(net->vargen,net->num_vargens);
  NET_set_vargen_array(net,VARGEN_array_new(num,net->num_periods),num);

  // Set buses
  NET_set_vargen_buses(net,bus_list);

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
}

void NET_add_batteries_from_params(Net* net, Bus* bus_list, REAL power_capacity,  REAL energy_capacity, REAL eta_c, REAL eta_d) {
  
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

void NET_dc_bus_hash_number_add(Net* net, BusDC* bus) {
  if (net)
    net->dc_bus_hash_number = BUSDC_hash_number_add(net->dc_bus_hash_number,bus);
}

BusDC* NET_dc_bus_hash_number_find(Net* net, int number) {
  if (net)
    return BUSDC_hash_number_find(net->dc_bus_hash_number,number);
  else
    return NULL;
}

void NET_dc_bus_hash_name_add(Net* net, BusDC* bus) {
  if (net)
    net->dc_bus_hash_name = BUSDC_hash_name_add(net->dc_bus_hash_name,bus);
}

BusDC* NET_dc_bus_hash_name_find(Net* net, char* name) {
  if (net)
    return BUSDC_hash_name_find(net->dc_bus_hash_name,name);
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

  // CSC convs

  // VSC convs

  // DC buses

  // DC branches

  // Facts

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
  BUSDC_hash_number_del(net->dc_bus_hash_number);
  BUSDC_hash_name_del(net->dc_bus_hash_name);
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;
  net->dc_bus_hash_number = NULL;
  net->dc_bus_hash_name = NULL;

  // Free components
  BUS_array_del(net->bus,net->num_buses);
  BRANCH_array_del(net->branch,net->num_branches);
  GEN_array_del(net->gen,net->num_gens);
  SHUNT_array_del(net->shunt,net->num_shunts);
  LOAD_array_del(net->load,net->num_loads);
  VARGEN_array_del(net->vargen,net->num_vargens);
  BAT_array_del(net->bat,net->num_bats);
  CONVCSC_array_del(net->csc_conv,net->num_csc_convs);
  CONVVSC_array_del(net->vsc_conv,net->num_vsc_convs);
  BUSDC_array_del(net->dc_bus,net->num_dc_buses);
  BRANCHDC_array_del(net->dc_branch,net->num_dc_branches);
  FACTS_array_del(net->facts,net->num_facts);

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
  ConvCSC* csc_conv;
  ConvVSC* vsc_conv;
  BusDC* dc_bus;
  BranchDC* dc_branch;
  Facts* facts;
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

  // CSC convs
  for (i = 0; i < net->num_csc_convs; i++) {
    csc_conv = CONVCSC_array_get(net->csc_conv,i);
    CONVCSC_clear_flags(csc_conv,FLAG_ALL);
  }

  // VSC convs
  for (i = 0; i < net->num_vsc_convs; i++) {
    vsc_conv = CONVVSC_array_get(net->vsc_conv,i);
    CONVVSC_clear_flags(vsc_conv,FLAG_ALL);
  }

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++) {
    dc_bus = BUSDC_array_get(net->dc_bus,i);
    BUSDC_clear_flags(dc_bus,FLAG_ALL);
  }

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++) {
    dc_branch = BRANCHDC_array_get(net->dc_branch,i);
    BRANCHDC_clear_flags(dc_branch,FLAG_ALL);
  }

  // Facts
  for (i = 0; i < net->num_facts; i++) {
    facts = FACTS_array_get(net->facts,i);
    FACTS_clear_flags(facts,FLAG_ALL);
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

    // Vargen

    // Battery

    // CSC converter

    // VSC converter

    // DC bus

    // DC branch

    // Facts

    // Actions
    net->num_actions[t] = 0;
  }

  // Mismatches
  for (i = 0; i < net->num_buses; i++)
    BUS_clear_mismatches(BUS_array_get(net->bus,i));
  for (i = 0; i < net->num_dc_buses; i++)
    BUSDC_clear_mismatches(BUSDC_array_get(net->dc_bus,i));
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
  for (i = 0; i < net->num_vargens; i++)
    VARGEN_clear_sensitivities(VARGEN_array_get(net->vargen,i));

  // Shunts
  for (i = 0; i < net->num_shunts; i++)
    SHUNT_clear_sensitivities(SHUNT_array_get(net->shunt,i));

  // Batteries
  for (i = 0; i < net->num_bats; i++)
    BAT_clear_sensitivities(BAT_array_get(net->bat,i));

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++)
    CONVCSC_clear_sensitivities(CONVCSC_array_get(net->csc_conv,i));

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++)
    CONVVSC_clear_sensitivities(CONVVSC_array_get(net->vsc_conv,i));

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++)
    BUSDC_clear_sensitivities(BUSDC_array_get(net->dc_bus,i));

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++)
    BRANCHDC_clear_sensitivities(BRANCHDC_array_get(net->dc_branch,i));

  // FACTS
  for (i = 0; i < net->num_facts; i++)
    FACTS_clear_sensitivities(FACTS_array_get(net->facts,i));
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

void NET_copy_from_net(Net* net, Net* other_net) {
  /** Copies data from another network except
   *  topological information.
   */
  
  // Local variables
  Bus* bus = NULL;
  Bus* other_bus = NULL;
  Branch* branch = NULL;
  Branch* other_branch = NULL;
  Gen* gen = NULL;
  Gen* other_gen = NULL;
  Vargen* vargen = NULL;
  Vargen* other_vargen = NULL;
  Load* load = NULL;
  Load* other_load = NULL;
  Shunt* shunt = NULL;
  Shunt* other_shunt = NULL;
  Bat* bat = NULL;
  Bat* other_bat = NULL;
  ConvCSC* csc_conv = NULL;
  ConvCSC* other_csc_conv = NULL;
  ConvVSC* vsc_conv = NULL;
  ConvVSC* other_vsc_conv = NULL;
  BusDC* dc_bus = NULL;
  BusDC* other_dc_bus = NULL;
  BranchDC* dc_branch = NULL;
  BranchDC* other_dc_branch = NULL;
  Facts* facts = NULL;
  Facts* other_facts = NULL;
  int i;

  // Check
  if (!net || !other_net)
    return;

  // Free hash tables
  BUS_hash_number_del(net->bus_hash_number);
  BUS_hash_name_del(net->bus_hash_name);
  BUSDC_hash_number_del(net->dc_bus_hash_number);
  BUSDC_hash_name_del(net->dc_bus_hash_name);
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;
  net->dc_bus_hash_number = NULL;
  net->dc_bus_hash_name = NULL;
  
  // Error
  net->error_flag = other_net->error_flag;
  strcpy(net->error_string, other_net->error_string);

  // Output
  strcpy(net->output_string, other_net->output_string);

  // Base power
  net->base_power = other_net->base_power;

  // Num flags
  net->num_vars = other_net->num_vars;
  net->num_fixed = other_net->num_fixed;
  net->num_bounded = other_net->num_bounded;
  net->num_sparse = other_net->num_sparse;  

  // Buses
  for (i = 0; i < net->num_buses; i++) {
    bus = NET_get_bus(net,i);
    other_bus = NET_get_bus(other_net,i);
    BUS_copy_from_bus(bus,other_bus);
    NET_bus_hash_number_add(net,bus);
    NET_bus_hash_name_add(net,bus);
  }

  // Branches
  for (i = 0; i < net->num_branches; i++) {  
    branch = NET_get_branch(net,i);
    other_branch = NET_get_branch(other_net,i);
    BRANCH_copy_from_branch(branch,other_branch);
  }

  // Generators
  for (i = 0; i < net->num_gens; i++) {  
    gen = NET_get_gen(net,i);
    other_gen = NET_get_gen(other_net,i);
    GEN_copy_from_gen(gen,other_gen);
  }

  // Var generators
  for (i = 0; i < net->num_vargens; i++) {  
    vargen = NET_get_vargen(net,i);
    other_vargen = NET_get_vargen(other_net,i);
    VARGEN_copy_from_vargen(vargen,other_vargen);
  }
  
  // Shunts
  for (i = 0; i < net->num_shunts; i++) { 
    shunt = NET_get_shunt(net,i);
    other_shunt = NET_get_shunt(other_net,i);
    SHUNT_copy_from_shunt(shunt,other_shunt);
  }

  // Loads
  for (i = 0; i < net->num_loads; i++) {
    load = NET_get_load(net,i);
    other_load = NET_get_load(other_net,i);
    LOAD_copy_from_load(load,other_load);
  }

  // Batteries
  for (i = 0; i < net->num_bats; i++) {
    bat = NET_get_bat(net,i);
    other_bat = NET_get_bat(other_net,i);
    BAT_copy_from_bat(bat,other_bat);
  }

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++) {
    csc_conv = NET_get_csc_conv(net,i);
    other_csc_conv = NET_get_csc_conv(other_net,i);
    CONVCSC_copy_from_conv(csc_conv,other_csc_conv);
  }

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++) {
    vsc_conv = NET_get_vsc_conv(net,i);
    other_vsc_conv = NET_get_vsc_conv(other_net,i);
    CONVVSC_copy_from_conv(vsc_conv,other_vsc_conv);
  }

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++) {
    dc_bus = NET_get_dc_bus(net,i);
    other_dc_bus = NET_get_dc_bus(other_net,i);
    BUSDC_copy_from_dc_bus(dc_bus,other_dc_bus);
    NET_dc_bus_hash_number_add(net,dc_bus);
    NET_dc_bus_hash_name_add(net,dc_bus);
  }

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++) {  
    dc_branch = NET_get_dc_branch(net,i);
    other_dc_branch = NET_get_dc_branch(other_net,i);
    BRANCHDC_copy_from_dc_branch(dc_branch,other_dc_branch);
  }

  // Facts
  for (i = 0; i < net->num_facts; i++) {  
    facts = NET_get_facts(net,i);
    other_facts = NET_get_facts(other_net,i);
    FACTS_copy_from_facts(facts,other_facts);
  }
}

Net* NET_get_copy(Net* net) {
  /** Gets deep copy of network.
   */

  // Local variables
  Net* new_net = NULL;
  Bus* bus = NULL;
  Bus* new_bus = NULL;
  Branch* branch = NULL;
  Branch* new_branch = NULL;
  Gen* gen = NULL;
  Gen* new_gen = NULL;
  Vargen* vargen = NULL;
  Vargen* new_vargen = NULL;
  Load* load = NULL;
  Load* new_load = NULL;
  Shunt* shunt = NULL;
  Shunt* new_shunt = NULL;
  Bat* bat = NULL;
  Bat* new_bat = NULL;
  ConvCSC* csc_conv = NULL;
  ConvCSC* new_csc_conv = NULL;
  ConvVSC* vsc_conv = NULL;
  ConvVSC* new_vsc_conv = NULL;
  BusDC* dc_bus = NULL;
  BusDC* new_dc_bus = NULL;
  BranchDC* dc_branch = NULL;
  BranchDC* new_dc_branch = NULL;
  Facts* facts = NULL;
  Facts* new_facts = NULL;
  int i;

  // Check
  if (!net)
    return new_net;

  // Allocate
  new_net = NET_new(net->num_periods);
  NET_set_bus_array(new_net,BUS_array_new(net->num_buses,net->num_periods), net->num_buses);
  NET_set_branch_array(new_net,BRANCH_array_new(net->num_branches,net->num_periods),net->num_branches);
  NET_set_gen_array(new_net,GEN_array_new(net->num_gens,net->num_periods),net->num_gens);
  NET_set_vargen_array(new_net,VARGEN_array_new(net->num_vargens,net->num_periods),net->num_vargens);
  NET_set_shunt_array(new_net,SHUNT_array_new(net->num_shunts,net->num_periods),net->num_shunts);
  NET_set_load_array(new_net,LOAD_array_new(net->num_loads,net->num_periods),net->num_loads);
  NET_set_bat_array(new_net,BAT_array_new(net->num_bats,net->num_periods),net->num_bats);
  NET_set_csc_conv_array(new_net,CONVCSC_array_new(net->num_csc_convs,net->num_periods),net->num_csc_convs);
  NET_set_vsc_conv_array(new_net,CONVVSC_array_new(net->num_vsc_convs,net->num_periods),net->num_vsc_convs);
  NET_set_dc_bus_array(new_net,BUSDC_array_new(net->num_dc_buses,net->num_periods), net->num_dc_buses);
  NET_set_dc_branch_array(new_net,BRANCHDC_array_new(net->num_dc_branches,net->num_periods), net->num_dc_branches);
  NET_set_facts_array(new_net,FACTS_array_new(net->num_facts,net->num_periods), net->num_facts);
  
  // Buses
  for (i = 0; i < net->num_buses; i++) {

    bus = NET_get_bus(net,i);
    new_bus = NET_get_bus(new_net,i);

    // Connections gen
    for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen))
      BUS_add_gen(new_bus,NET_get_gen(new_net,GEN_get_index(gen)));

    // Connections reg gen
    for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen))
      BUS_add_reg_gen(new_bus,NET_get_gen(new_net,GEN_get_index(gen)));

    // Connections load
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load))
      BUS_add_load(new_bus,NET_get_load(new_net,LOAD_get_index(load)));

    // Connections shunt
    for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt))
      BUS_add_shunt(new_bus,NET_get_shunt(new_net,SHUNT_get_index(shunt)));

    // Connections reg_shunt
    for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt))
      BUS_add_reg_shunt(new_bus,NET_get_shunt(new_net,SHUNT_get_index(shunt)));

    // Connections branch_k
    for (branch = BUS_get_branch_k(bus); branch != NULL; branch = BRANCH_get_next_k(branch))
      BUS_add_branch_k(new_bus,NET_get_branch(new_net,BRANCH_get_index(branch)));

    // Connections branch_m
    for (branch = BUS_get_branch_m(bus); branch != NULL; branch = BRANCH_get_next_m(branch))
      BUS_add_branch_m(new_bus,NET_get_branch(new_net,BRANCH_get_index(branch)));

    // Connections reg_tran
    for (branch = BUS_get_reg_tran(bus); branch != NULL; branch = BRANCH_get_reg_next(branch))
      BUS_add_reg_tran(new_bus,NET_get_branch(new_net,BRANCH_get_index(branch)));

    // Connections vargen
    for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen))
      BUS_add_vargen(new_bus, NET_get_vargen(new_net,VARGEN_get_index(vargen)));

    // Connections bat
    for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat))
      BUS_add_bat(new_bus,NET_get_bat(new_net,BAT_get_index(bat)));

    // Connections csc conv
    for (csc_conv = BUS_get_csc_conv(bus); csc_conv != NULL; csc_conv = CONVCSC_get_next_ac(csc_conv))
      BUS_add_csc_conv(new_bus,NET_get_csc_conv(new_net,CONVCSC_get_index(csc_conv)));

    // Connections vsc conv
    for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv))
      BUS_add_vsc_conv(new_bus,NET_get_vsc_conv(new_net,CONVVSC_get_index(vsc_conv)));

    // Connections reg vsc conv
    for (vsc_conv = BUS_get_reg_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_reg_next(vsc_conv))
      BUS_add_reg_vsc_conv(new_bus,NET_get_vsc_conv(new_net,CONVVSC_get_index(vsc_conv)));

    // Connections facts_k
    for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts))
      BUS_add_facts_k(new_bus,NET_get_facts(new_net,FACTS_get_index(facts)));

    // Connections facts_m
    for (facts = BUS_get_facts_m(bus); facts != NULL; facts = FACTS_get_next_m(facts))
      BUS_add_facts_m(new_bus,NET_get_facts(new_net,FACTS_get_index(facts)));

    // Connections reg_facts
    for (facts = BUS_get_reg_facts(bus); facts != NULL; facts = FACTS_get_reg_next(facts))
      BUS_add_reg_facts(new_bus,NET_get_facts(new_net,FACTS_get_index(facts)));
  }

  // Branches
  for (i = 0; i < net->num_branches; i++) {

    branch = NET_get_branch(net,i);
    new_branch = NET_get_branch(new_net,i);
    
    // Connections bus
    BRANCH_set_bus_k(new_branch,NET_get_bus(new_net,BUS_get_index(BRANCH_get_bus_k(branch))));
    BRANCH_set_bus_m(new_branch,NET_get_bus(new_net,BUS_get_index(BRANCH_get_bus_m(branch))));
    BRANCH_set_reg_bus(new_branch,NET_get_bus(new_net,BUS_get_index(BRANCH_get_reg_bus(branch))));
  }

  // Generators
  for (i = 0; i < net->num_gens; i++) {

    gen = NET_get_gen(net,i);
    new_gen = NET_get_gen(new_net,i);
    
    // Connections bus
    GEN_set_bus(new_gen,NET_get_bus(new_net,BUS_get_index(GEN_get_bus(gen))));
    GEN_set_reg_bus(new_gen,NET_get_bus(new_net,BUS_get_index(GEN_get_reg_bus(gen))));
  }

  // Var generators
  for (i = 0; i < net->num_vargens; i++) {

    vargen = NET_get_vargen(net,i);
    new_vargen = NET_get_vargen(new_net,i);

    // Connections bus
    VARGEN_set_bus(new_vargen,NET_get_bus(new_net,BUS_get_index(VARGEN_get_bus(vargen))));
  }

  // Shunts
  for (i = 0; i < net->num_shunts; i++) {

    shunt = NET_get_shunt(net,i);
    new_shunt = NET_get_shunt(new_net,i);

    // Connections bus
    SHUNT_set_bus(new_shunt,NET_get_bus(new_net,BUS_get_index(SHUNT_get_bus(shunt))));
    SHUNT_set_reg_bus(new_shunt,NET_get_bus(new_net,BUS_get_index(SHUNT_get_reg_bus(shunt))));
  }

  // Loads
  for (i = 0; i < net->num_loads; i++) {

    load = NET_get_load(net,i);
    new_load = NET_get_load(new_net,i);
    
    // Connections bus
    LOAD_set_bus(new_load,NET_get_bus(new_net,BUS_get_index(LOAD_get_bus(load))));
  }

  // Batteries
  for (i = 0; i < net->num_bats; i++) {

    bat = NET_get_bat(net,i);
    new_bat = NET_get_bat(new_net,i);
    
    // Connections bus
    BAT_set_bus(new_bat,NET_get_bus(new_net,BUS_get_index(BAT_get_bus(bat))));
  }

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++) {

    csc_conv = NET_get_csc_conv(net,i);
    new_csc_conv = NET_get_csc_conv(new_net,i);
    
    // Connections ac bus
    CONVCSC_set_ac_bus(new_csc_conv,NET_get_bus(new_net,BUS_get_index(CONVCSC_get_ac_bus(csc_conv))));

    // Connections dc bus
    CONVCSC_set_dc_bus(new_csc_conv,NET_get_dc_bus(new_net,BUSDC_get_index(CONVCSC_get_dc_bus(csc_conv))));
  }

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++) {

    vsc_conv = NET_get_vsc_conv(net,i);
    new_vsc_conv = NET_get_vsc_conv(new_net,i);
    
    // Connections ac bus
    CONVVSC_set_ac_bus(new_vsc_conv,NET_get_bus(new_net,BUS_get_index(CONVVSC_get_ac_bus(vsc_conv))));
    CONVVSC_set_reg_bus(new_vsc_conv,NET_get_bus(new_net,BUS_get_index(CONVVSC_get_reg_bus(vsc_conv))));
    
    // Connections dc bus
    CONVVSC_set_dc_bus(new_vsc_conv,NET_get_dc_bus(new_net,BUSDC_get_index(CONVVSC_get_dc_bus(vsc_conv))));
  }

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++) {

    dc_bus = NET_get_dc_bus(net,i);
    new_dc_bus = NET_get_dc_bus(new_net,i);

    // Connections branch_k
    for (dc_branch = BUSDC_get_branch_k(dc_bus); dc_branch != NULL; dc_branch = BRANCHDC_get_next_k(dc_branch))
      BUSDC_add_branch_k(new_dc_bus,NET_get_dc_branch(new_net,BRANCHDC_get_index(dc_branch)));

    // Connections branch_m
    for (dc_branch = BUSDC_get_branch_m(dc_bus); dc_branch != NULL; dc_branch = BRANCHDC_get_next_m(dc_branch))
      BUSDC_add_branch_m(new_dc_bus,NET_get_dc_branch(new_net,BRANCHDC_get_index(dc_branch)));

    // Connections csc conv
    for (csc_conv = BUSDC_get_csc_conv(dc_bus); csc_conv != NULL; csc_conv = CONVCSC_get_next_dc(csc_conv))
      BUSDC_add_csc_conv(new_dc_bus,NET_get_csc_conv(new_net,CONVCSC_get_index(csc_conv)));

    // Connections vsc conv
    for (vsc_conv = BUSDC_get_vsc_conv(dc_bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_dc(vsc_conv))
      BUSDC_add_vsc_conv(new_dc_bus,NET_get_vsc_conv(new_net,CONVVSC_get_index(vsc_conv)));
  }

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++) {

    dc_branch = NET_get_dc_branch(net,i);
    new_dc_branch = NET_get_dc_branch(new_net,i);
    
    // Connections bus
    BRANCHDC_set_bus_k(new_dc_branch,NET_get_dc_bus(new_net,BUSDC_get_index(BRANCHDC_get_bus_k(dc_branch))));
    BRANCHDC_set_bus_m(new_dc_branch,NET_get_dc_bus(new_net,BUSDC_get_index(BRANCHDC_get_bus_m(dc_branch))));
  }

  // Facts
  for (i = 0; i < net->num_facts; i++) {

    facts = NET_get_facts(net,i);
    new_facts = NET_get_facts(new_net,i);
    
    // Connections bus
    FACTS_set_bus_k(new_facts,NET_get_bus(new_net,BUS_get_index(FACTS_get_bus_k(facts))));
    FACTS_set_bus_m(new_facts,NET_get_bus(new_net,BUS_get_index(FACTS_get_bus_m(facts))));
    FACTS_set_reg_bus(new_facts,NET_get_bus(new_net,BUS_get_index(FACTS_get_reg_bus(facts))));
  }
  
  // Copy rest of data
  NET_copy_from_net(new_net,net);
    
  // Return
  return new_net;
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

Net* NET_extract_subnet(Net* net, Bus** bus_ptr_array, int size) {
  /** Extracts subnetwork containing the specified buses.
   *  It does not include any HVDC components.
   */

  // Local variables
  Net* new_net;
  Bus* bus;
  Gen* gen;
  Branch* br;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  Vargen* vargen;
  Facts* facts;
  
  char* keep_bus;
  char* delete_br;
  char* delete_facts;

  Bus** del_bus_ptr_array;
  Gen** del_gen_ptr_array;
  Branch** del_br_ptr_array;
  Shunt** del_shunt_ptr_array;
  Load** del_load_ptr_array;
  Bat** del_bat_ptr_array;
  Vargen** del_vargen_ptr_array;
  ConvCSC** del_csc_conv_ptr_array;
  ConvVSC** del_vsc_conv_ptr_array;
  BusDC** del_dc_bus_ptr_array;
  BranchDC** del_dc_br_ptr_array;
  Facts** del_facts_ptr_array;
  
  int num_del_bus;
  int num_del_br;
  int num_del_gen;
  int num_del_shunt;
  int num_del_load;
  int num_del_bat;
  int num_del_vargen;
  int num_del_facts;

  int i;

  // Check
  if (!net || !bus_ptr_array)
    return NULL;

  // Copy
  new_net = NET_get_copy(net);

  // Flags
  ARRAY_zalloc(keep_bus,char,net->num_buses);
  ARRAY_zalloc(delete_br,char,net->num_branches);
  ARRAY_zalloc(delete_facts,char,net->num_facts);

  // Del ptr arrays
  del_bus_ptr_array = (Bus**)malloc(sizeof(Bus*)*net->num_buses);
  del_br_ptr_array = (Branch**)malloc(sizeof(Branch*)*net->num_branches);
  del_gen_ptr_array = (Gen**)malloc(sizeof(Gen*)*net->num_gens);
  del_load_ptr_array = (Load**)malloc(sizeof(Load*)*net->num_loads);
  del_shunt_ptr_array = (Shunt**)malloc(sizeof(Shunt*)*net->num_shunts);
  del_bat_ptr_array = (Bat**)malloc(sizeof(Bat*)*net->num_bats);
  del_vargen_ptr_array = (Vargen**)malloc(sizeof(Vargen*)*net->num_vargens);
  del_csc_conv_ptr_array = (ConvCSC**)malloc(sizeof(ConvCSC*)*net->num_csc_convs);
  del_vsc_conv_ptr_array = (ConvVSC**)malloc(sizeof(ConvVSC*)*net->num_vsc_convs);
  del_dc_bus_ptr_array = (BusDC**)malloc(sizeof(BusDC*)*net->num_dc_buses);
  del_dc_br_ptr_array = (BranchDC**)malloc(sizeof(BranchDC*)*net->num_dc_branches);
  del_facts_ptr_array = (Facts**)malloc(sizeof(Facts*)*net->num_facts);

  // Mark buses to keep
  for (i = 0; i < size; i++) {
    
    bus = bus_ptr_array[i];

    // No bus
    if (!bus)
      continue;

    // Not on network
    if (bus != NET_get_bus(net, BUS_get_index(bus)))
      continue;
    
    keep_bus[BUS_get_index(bus)] = 1;	
  }

  // Fill del ptr arrays
  num_del_bus = 0;
  num_del_br = 0;
  num_del_gen = 0;
  num_del_shunt =0;
  num_del_load = 0;
  num_del_bat = 0;
  num_del_vargen = 0;
  num_del_facts = 0;
  for (i = 0; i < new_net->num_buses; i++) {
    
    bus = NET_get_bus(new_net,i);

    if (keep_bus[BUS_get_index(bus)])
      continue;

    // Bus
    del_bus_ptr_array[num_del_bus] = bus;
    num_del_bus++;

    // Generators
    for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
      del_gen_ptr_array[num_del_gen] = gen;
      num_del_gen++;
    }

    // Loads
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
      del_load_ptr_array[num_del_load] = load;
      num_del_load++;
    }

    // Shunts
    for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
      del_shunt_ptr_array[num_del_shunt] = shunt;
      num_del_shunt++;
    }

    // Batteries
    for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
      del_bat_ptr_array[num_del_bat] = bat;
      num_del_bat++;
    }

    // Var generators
    for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
      del_vargen_ptr_array[num_del_vargen] = vargen;
      num_del_vargen++;
    }
   
    // Branches k
    for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {
      if (!delete_br[BRANCH_get_index(br)]) {
        delete_br[BRANCH_get_index(br)] = 1;
        del_br_ptr_array[num_del_br] = br;
        num_del_br++;
      }
    }

    // Branches m
    for (br = BUS_get_branch_m(bus); br != NULL; br = BRANCH_get_next_m(br)) {
      if (!delete_br[BRANCH_get_index(br)]) {
        delete_br[BRANCH_get_index(br)] = 1;
        del_br_ptr_array[num_del_br] = br;
        num_del_br++;
      }
    }

    // Facts k
    for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {
      if (!delete_facts[FACTS_get_index(facts)]) {
        delete_facts[FACTS_get_index(facts)] = 1;
        del_facts_ptr_array[num_del_facts] = facts;
        num_del_facts++;
      }
    }

    // Facts m
    for (facts = BUS_get_facts_m(bus); facts != NULL; facts = FACTS_get_next_m(facts)) {
      if (!delete_facts[FACTS_get_index(facts)]) {
        delete_facts[FACTS_get_index(facts)] = 1;
        del_facts_ptr_array[num_del_facts] = facts;
        num_del_facts++;
      }
    }
  }

  // Delete all CSC converters
  for (i = 0; i < new_net->num_csc_convs; i++)
    del_csc_conv_ptr_array[i] = NET_get_csc_conv(new_net, i);

  // Delete all VSC converters
  for (i = 0; i < new_net->num_vsc_convs; i++)
    del_vsc_conv_ptr_array[i] = NET_get_vsc_conv(new_net, i);

  // Delete all DC buses
  for (i = 0; i < new_net->num_dc_buses; i++)
    del_dc_bus_ptr_array[i] = NET_get_dc_bus(new_net, i);

  // Delete all DC branches
  for (i = 0; i < new_net->num_dc_branches; i++)
    del_dc_br_ptr_array[i] = NET_get_dc_branch(new_net, i);

  // Delete components
  NET_del_buses(new_net, del_bus_ptr_array, num_del_bus);
  NET_del_branches(new_net, del_br_ptr_array, num_del_br);
  NET_del_gens(new_net, del_gen_ptr_array, num_del_gen);
  NET_del_loads(new_net, del_load_ptr_array, num_del_load);
  NET_del_shunts(new_net, del_shunt_ptr_array, num_del_shunt);
  NET_del_bats(new_net, del_bat_ptr_array, num_del_bat);
  NET_del_vargens(new_net, del_vargen_ptr_array, num_del_vargen);
  NET_del_csc_convs(new_net,del_csc_conv_ptr_array,new_net->num_csc_convs);
  NET_del_vsc_convs(new_net,del_vsc_conv_ptr_array,new_net->num_vsc_convs);
  NET_del_dc_buses(new_net,del_dc_bus_ptr_array,new_net->num_dc_buses);
  NET_del_dc_branches(new_net,del_dc_br_ptr_array,new_net->num_dc_branches);
  NET_del_facts(new_net,del_facts_ptr_array,num_del_facts);

  // Clean flag arrays
  free(keep_bus);
  free(delete_br);
  free(delete_facts);

  // Clean del ptr arrays
  free(del_bus_ptr_array);
  free(del_br_ptr_array);
  free(del_gen_ptr_array);
  free(del_load_ptr_array);
  free(del_shunt_ptr_array);
  free(del_bat_ptr_array);
  free(del_vargen_ptr_array);
  free(del_csc_conv_ptr_array);
  free(del_vsc_conv_ptr_array);
  free(del_dc_bus_ptr_array);
  free(del_dc_br_ptr_array);
  free(del_facts_ptr_array);
  
  // Clear flags
  NET_clear_flags(net);
  
  // Return
  return new_net;
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
  net->csc_conv = NULL;
  net->vsc_conv = NULL;
  net->dc_bus = NULL;
  net->dc_branch = NULL;
  net->facts = NULL;

  // Hash tables
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;
  net->dc_bus_hash_number = NULL;
  net->dc_bus_hash_name = NULL;

  // Number of components
  net->num_buses = 0;
  net->num_branches = 0;
  net->num_gens = 0;
  net->num_loads = 0;
  net->num_shunts = 0;
  net->num_vargens = 0;
  net->num_bats = 0;
  net->num_csc_convs = 0;
  net->num_vsc_convs = 0;
  net->num_dc_buses = 0;
  net->num_dc_branches = 0;
  net->num_facts = 0;

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

  // State tag
  net->state_tag = 0;
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

BusDC* NET_get_dc_bus_hash_number(Net* net) {
  if (!net)
    return NULL;
  else
    return net->dc_bus_hash_number;
}

BusDC* NET_get_dc_bus_hash_name(Net* net) {
  if (!net)
    return NULL;
  else
    return net->dc_bus_hash_name;
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

ConvCSC* NET_get_csc_conv(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_csc_convs)
    return NULL;
  else
    return CONVCSC_array_get(net->csc_conv,index);
}

ConvVSC* NET_get_vsc_conv(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_vsc_convs)
    return NULL;
  else
    return CONVVSC_array_get(net->vsc_conv,index);
}

BusDC* NET_get_dc_bus(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_dc_buses)
    return NULL;
  else
    return BUSDC_array_get(net->dc_bus,index);
}

BranchDC* NET_get_dc_branch(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_dc_branches)
    return NULL;
  else
    return BRANCHDC_array_get(net->dc_branch,index);
}

Facts* NET_get_facts(Net* net, int index) {
  if (!net || index < 0 || index >= net->num_facts)
    return NULL;
  else
    return FACTS_array_get(net->facts,index);
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

Gen* NET_get_gen_from_name_and_bus_number(Net* net, char* name, int number) {
  Gen* gen;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    if (strcmp(GEN_get_name(gen), name) == 0)
      return gen;
  }
  return NULL;
}

Branch* NET_get_branch_from_name_and_bus_numbers(Net* net, char* name, int number1, int number2) {
  Branch* br;
  Bus* bus1 = NET_bus_hash_number_find(net, number1);
  Bus* bus2 = NET_bus_hash_number_find(net, number2);
  for (br = BUS_get_branch_k(bus1); br != NULL; br = BRANCH_get_next_k(br)) {
    if (bus2 == BRANCH_get_bus_m(br) &&
        strcmp(BRANCH_get_name(br), name) == 0)
      return br;
  }
  for (br = BUS_get_branch_m(bus1); br != NULL; br = BRANCH_get_next_m(br)) {
    if (bus2 == BRANCH_get_bus_k(br) &&
        strcmp(BRANCH_get_name(br), name) == 0)
      return br;
  }
  return NULL;
}

Shunt* NET_get_shunt_from_name_and_bus_number(Net* net, char* name, int number) {
  Shunt* shunt;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    if (strcmp(SHUNT_get_name(shunt), name) == 0)
      return shunt;
  }
  return NULL;
}

Shunt* NET_get_fixed_shunt_from_name_and_bus_number(Net* net, char* name, int number) {
  Shunt* shunt;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    if (SHUNT_is_fixed(shunt) && strcmp(SHUNT_get_name(shunt), name) == 0)
      return shunt;
  }
  return NULL;
}

Shunt* NET_get_switched_shunt_from_name_and_bus_number(Net* net, char* name, int number) {
  Shunt* shunt;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    if (SHUNT_is_switched(shunt) && strcmp(SHUNT_get_name(shunt), name) == 0)
      return shunt;
  }
  return NULL;
}

Load* NET_get_load_from_name_and_bus_number(Net* net, char* name, int number) {
  Load* load;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    if (strcmp(LOAD_get_name(load), name) == 0)
      return load;
  }
  return NULL;
}

Vargen* NET_get_vargen_from_name_and_bus_number(Net* net, char* name, int number) {
  Vargen* vargen;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    if (strcmp(VARGEN_get_name(vargen), name) == 0)
      return vargen;
  }
  return NULL;
}

Bat* NET_get_bat_from_name_and_bus_number(Net* net, char* name, int number) {
  Bat* bat;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    if (strcmp(BAT_get_name(bat), name) == 0)
      return bat;
  }
  return NULL;
}

ConvCSC* NET_get_csc_conv_from_name_and_ac_bus_number(Net* net, char* name, int number) {
  ConvCSC* conv;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (conv = BUS_get_csc_conv(bus); conv != NULL; conv = CONVCSC_get_next_ac(conv)) {
    if (strcmp(CONVCSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

ConvCSC* NET_get_csc_conv_from_name_and_dc_bus_number(Net* net, char* name, int number) {
  ConvCSC* conv;
  BusDC* bus = NET_dc_bus_hash_number_find(net, number);
  for (conv = BUSDC_get_csc_conv(bus); conv != NULL; conv = CONVCSC_get_next_dc(conv)) {
    if (strcmp(CONVCSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

ConvCSC* NET_get_csc_conv_from_name_and_dc_bus_name(Net* net, char* name, char* bus_name) {
  ConvCSC* conv;
  BusDC* bus = NET_dc_bus_hash_name_find(net, bus_name);
  for (conv = BUSDC_get_csc_conv(bus); conv != NULL; conv = CONVCSC_get_next_dc(conv)) {
    if (strcmp(CONVCSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

ConvVSC* NET_get_vsc_conv_from_name_and_ac_bus_number(Net* net, char* name, int number) {
  ConvVSC* conv;
  Bus* bus = NET_bus_hash_number_find(net, number);
  for (conv = BUS_get_vsc_conv(bus); conv != NULL; conv = CONVVSC_get_next_ac(conv)) {
    if (strcmp(CONVVSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

ConvVSC* NET_get_vsc_conv_from_name_and_dc_bus_number(Net* net, char* name, int number) {
  ConvVSC* conv;
  BusDC* bus = NET_dc_bus_hash_number_find(net, number);
  for (conv = BUSDC_get_vsc_conv(bus); conv != NULL; conv = CONVVSC_get_next_dc(conv)) {
    if (strcmp(CONVVSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

ConvVSC* NET_get_vsc_conv_from_name_and_dc_bus_name(Net* net, char* name, char* bus_name) {
  ConvVSC* conv;
  BusDC* bus = NET_dc_bus_hash_name_find(net, bus_name);
  for (conv = BUSDC_get_vsc_conv(bus); conv != NULL; conv = CONVVSC_get_next_dc(conv)) {
    if (strcmp(CONVVSC_get_name(conv), name) == 0)
      return conv;
  }
  return NULL;
}

BranchDC* NET_get_dc_branch_from_name_and_dc_bus_numbers(Net* net, char* name, int number1, int number2) {
  BranchDC* br;
  BusDC* bus1 = NET_dc_bus_hash_number_find(net, number1);
  BusDC* bus2 = NET_dc_bus_hash_number_find(net, number2);
  for (br = BUSDC_get_branch_k(bus1); br != NULL; br = BRANCHDC_get_next_k(br)) {
    if (bus2 == BRANCHDC_get_bus_m(br) && strcmp(BRANCHDC_get_name(br), name) == 0)
      return br;
  }
  for (br = BUSDC_get_branch_m(bus1); br != NULL; br = BRANCHDC_get_next_m(br)) {
    if (bus2 == BRANCHDC_get_bus_k(br) && strcmp(BRANCHDC_get_name(br), name) == 0)
      return br;
  }
  return NULL;
}

BranchDC* NET_get_dc_branch_from_name_and_dc_bus_names(Net* net, char* name, char* bus1_name, char* bus2_name) {
  BranchDC* br;
  BusDC* bus1 = NET_dc_bus_hash_name_find(net, bus1_name);
  BusDC* bus2 = NET_dc_bus_hash_name_find(net, bus2_name);
  for (br = BUSDC_get_branch_k(bus1); br != NULL; br = BRANCHDC_get_next_k(br)) {
    if (bus2 == BRANCHDC_get_bus_m(br) && strcmp(BRANCHDC_get_name(br), name) == 0)
      return br;
  }
  for (br = BUSDC_get_branch_m(bus1); br != NULL; br = BRANCHDC_get_next_m(br)) {
    if (bus2 == BRANCHDC_get_bus_k(br) && strcmp(BRANCHDC_get_name(br), name) == 0)
      return br;
  }
  return NULL;
}

Facts* NET_get_facts_from_name_and_bus_numbers(Net* net, char* name, int number1, int number2) {
  Facts* f;
  Bus* bus1 = NET_bus_hash_number_find(net,number1);
  Bus* bus2 = NET_bus_hash_number_find(net,number2);
  for (f = BUS_get_facts_k(bus1); f != NULL; f = FACTS_get_next_k(f)) {
    if (bus2 == FACTS_get_bus_m(f) &&
        strcmp(FACTS_get_name(f), name) == 0)
      return f;
  }
  for (f = BUS_get_facts_m(bus1); f != NULL; f = FACTS_get_next_m(f)) {
    if (bus2 == FACTS_get_bus_k(f) &&
        strcmp(FACTS_get_name(f), name) == 0)
      return f;
  }
  return NULL;
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

int NET_get_num_star_buses(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_star(BUS_array_get(net->bus,i)))
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

int NET_get_num_buses_reg_by_vsc_conv(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_regulated_by_vsc_conv(BUS_array_get(net->bus,i)))
      n++;
  }
  return n;
}

int NET_get_num_buses_reg_by_facts(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_buses; i++) {
    if (BUS_is_regulated_by_facts(BUS_array_get(net->bus,i)))
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

int NET_get_num_branches_on_outage(Net* net) {
  if (net)
    return net->num_branches-NET_get_num_branches_not_on_outage(net);
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

int NET_get_num_gens_on_outage(Net* net) {
  if (net)
    return net->num_gens-NET_get_num_gens_not_on_outage(net);
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

int NET_get_num_vdep_loads(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_loads; i++) {
    if (LOAD_is_vdep(LOAD_array_get(net->load,i)))
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

int NET_get_num_switched_v_shunts(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for (i = 0; i < net->num_shunts; i++) {
    if (SHUNT_is_switched_v(SHUNT_array_get(net->shunt,i)))
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

int NET_get_num_csc_convs(Net* net) {
  if (net)
    return net->num_csc_convs;
  else
    return 0;
}

int NET_get_num_vsc_convs(Net* net) {
  if (net)
    return net->num_vsc_convs;
  else
    return 0;
}

int NET_get_num_vsc_convs_in_P_dc_mode(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_vsc_convs; i++) {
    if (CONVVSC_is_in_P_dc_mode(NET_get_vsc_conv(net,i)))
      n++;
  }
  return n;
}

int NET_get_num_vsc_convs_in_v_dc_mode(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_vsc_convs; i++) {
    if (CONVVSC_is_in_v_dc_mode(NET_get_vsc_conv(net,i)))
      n++;
  }
  return n;
}

int NET_get_num_vsc_convs_in_v_ac_mode(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_vsc_convs; i++) {
    if (CONVVSC_is_in_v_ac_mode(NET_get_vsc_conv(net,i)))
      n++;
  }
  return n;
}

int NET_get_num_vsc_convs_in_f_ac_mode(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_vsc_convs; i++) {
    if (CONVVSC_is_in_f_ac_mode(NET_get_vsc_conv(net,i)))
      n++;
  }
  return n;
}

int NET_get_num_dc_buses(Net* net) {
  if (net)
    return net->num_dc_buses;
  else
    return 0;
}

int NET_get_num_dc_branches(Net* net) {
  if (net)
    return net->num_dc_branches;
  else
    return 0;
}

int NET_get_num_facts(Net* net) {
  if (net)
    return net->num_facts;
  else
    return 0;
}

int NET_get_num_facts_in_normal_series_mode(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_facts; i++) {
    if (FACTS_is_in_normal_series_mode(NET_get_facts(net,i)))
      n++;
  }
  return n;
}

int NET_get_num_reg_facts(Net* net) {
  int i;
  int n = 0;
  if (!net)
    return 0;
  for(i = 0; i < net->num_facts; i++) {
    if (FACTS_is_regulator(FACTS_array_get(net->facts,i)))
      n++;
  }
  return n;
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

  // No net
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

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++)
    CONVCSC_get_var_values(CONVCSC_array_get(net->csc_conv,i),values,code);

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++)
    CONVVSC_get_var_values(CONVVSC_array_get(net->vsc_conv,i),values,code);

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++)
    BUSDC_get_var_values(BUSDC_array_get(net->dc_bus,i),values,code);

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++)
    BRANCHDC_get_var_values(BRANCHDC_array_get(net->dc_branch,i),values,code);

  // Facts
  for (i = 0; i < net->num_facts; i++)
    FACTS_get_var_values(FACTS_array_get(net->facts,i),values,code);
  
  // Return
  return values;
}

char* NET_get_var_info_string(Net* net, int index) {
  
  // Local variables
  char* info;
  int i;

  // No net
  if (!net)
    return NULL;

  // Buses
  for (i = 0; i < net->num_buses; i++) {
    info = BUS_get_var_info_string(BUS_array_get(net->bus,i),index);
    if (info)
      return info;
  }

  // Generators
  for (i = 0; i < net->num_gens; i++) {
    info = GEN_get_var_info_string(GEN_array_get(net->gen,i),index);
    if (info)
      return info;
  }

  // Branches
  for (i = 0; i < net->num_branches; i++) {
    info = BRANCH_get_var_info_string(BRANCH_array_get(net->branch,i),index);
    if (info)
      return info;
  }

  // Shunts
  for (i = 0; i < net->num_shunts; i++) {
    info = SHUNT_get_var_info_string(SHUNT_array_get(net->shunt,i),index);
    if (info)
      return info;
  }

  // Loads
  for (i = 0; i < net->num_loads; i++) {
    info = LOAD_get_var_info_string(LOAD_array_get(net->load,i),index);
    if (info)
      return info;
  }
  
  // Variable generators
  for (i = 0; i < net->num_vargens; i++) {
    info = VARGEN_get_var_info_string(VARGEN_array_get(net->vargen,i),index);
    if (info)
      return info;
  }

  // Batteries
  for (i = 0; i < net->num_bats; i++) {
    info = BAT_get_var_info_string(BAT_array_get(net->bat,i),index);
    if (info)
      return info;
  }

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++) {
    info = CONVCSC_get_var_info_string(CONVCSC_array_get(net->csc_conv,i),index);
    if (info)
      return info;
  }

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++) {
    info = CONVVSC_get_var_info_string(CONVVSC_array_get(net->vsc_conv,i),index);
    if (info)
      return info;
  }

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++) {
    info = BUSDC_get_var_info_string(BUSDC_array_get(net->dc_bus,i),index);
    if (info)
      return info;
  }

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++) {
    info = BRANCHDC_get_var_info_string(BRANCHDC_array_get(net->dc_branch,i),index);
    if (info)
      return info;
  }

  // Facts
  for (i = 0; i < net->num_facts; i++) {
    info = FACTS_get_var_info_string(FACTS_array_get(net->facts,i),index);
    if (info)
      return info;
  }

  // Return
  return NULL;  
}

Mat* NET_get_var_projection(Net* net, char obj_type, char prop_mask, unsigned char var, int t_start, int t_end) {

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
  if ((obj_type == OBJ_ALL) && ((var != ALL_VARS) || (prop_mask != ANY_PROP))) {
    sprintf(net->error_string,"component-specific flag or properties cannot be used on all components");
    net->error_flag = TRUE;
    return NULL;
  }

  // Count
  num_subvars = 0;
  if ((obj_type == OBJ_BUS) || (obj_type == OBJ_ALL)) { // Bus or all
    for (i = 0; i < net->num_buses; i++) {
      if (BUS_has_properties(NET_get_bus(net,i),prop_mask))
        num_subvars += BUS_get_num_vars(NET_get_bus(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_GEN) || (obj_type == OBJ_ALL)) { // Gen or all
    for (i = 0; i < net->num_gens; i++) {
      if (GEN_has_properties(NET_get_gen(net,i),prop_mask))
        num_subvars += GEN_get_num_vars(NET_get_gen(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_LOAD) || (obj_type == OBJ_ALL)) { // Load or all
    for (i = 0; i < net->num_loads; i++) {
      if (LOAD_has_properties(NET_get_load(net,i),prop_mask))
        num_subvars += LOAD_get_num_vars(NET_get_load(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_BRANCH) || (obj_type == OBJ_ALL)) { // Branch or all
    for (i = 0; i < net->num_branches; i++) {
      if (BRANCH_has_properties(NET_get_branch(net,i),prop_mask))
        num_subvars += BRANCH_get_num_vars(NET_get_branch(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_SHUNT) || (obj_type == OBJ_ALL)) { // Shunt or all
    for (i = 0; i < net->num_shunts; i++) {
      if (SHUNT_has_properties(NET_get_shunt(net,i),prop_mask))
        num_subvars += SHUNT_get_num_vars(NET_get_shunt(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_VARGEN) || (obj_type == OBJ_ALL)) { // Vargen or all
    for (i = 0; i < net->num_vargens; i++) {
      if (VARGEN_has_properties(NET_get_vargen(net,i),prop_mask))
        num_subvars += VARGEN_get_num_vars(NET_get_vargen(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_BAT) || (obj_type == OBJ_ALL)) { // Battery or all
    for (i = 0; i < net->num_bats; i++) {
      if (BAT_has_properties(NET_get_bat(net,i),prop_mask))
        num_subvars += BAT_get_num_vars(NET_get_bat(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_CONVCSC) || (obj_type == OBJ_ALL)) { // CSC converter or all
    for (i = 0; i < net->num_csc_convs; i++) {
      if (CONVCSC_has_properties(NET_get_csc_conv(net,i),prop_mask))
        num_subvars += CONVCSC_get_num_vars(NET_get_csc_conv(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_CONVVSC) || (obj_type == OBJ_ALL)) { // VSC converter or all
    for (i = 0; i < net->num_vsc_convs; i++) {
      if (CONVVSC_has_properties(NET_get_vsc_conv(net,i),prop_mask))
        num_subvars += CONVVSC_get_num_vars(NET_get_vsc_conv(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_BUSDC) || (obj_type == OBJ_ALL)) { // DC bus or all
    for (i = 0; i < net->num_dc_buses; i++) {
      if (BUSDC_has_properties(NET_get_dc_bus(net,i),prop_mask))
        num_subvars += BUSDC_get_num_vars(NET_get_dc_bus(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_BRANCHDC) || (obj_type == OBJ_ALL)) { // DC branch or all
    for (i = 0; i < net->num_dc_branches; i++) {
      if (BRANCHDC_has_properties(NET_get_dc_branch(net,i),prop_mask))
        num_subvars += BRANCHDC_get_num_vars(NET_get_dc_branch(net,i),var,t_start,t_end);
    }
  }
  if ((obj_type == OBJ_FACTS) || (obj_type == OBJ_ALL)) { // Facts or all
    for (i = 0; i < net->num_facts; i++) {
      if (FACTS_has_properties(NET_get_facts(net,i),prop_mask))
        num_subvars += FACTS_get_num_vars(NET_get_facts(net,i),var,t_start,t_end);
    }
  }

  // Allocate
  proj = MAT_new(num_subvars,
                 net->num_vars,
                 num_subvars);

  // Fill
  num_subvars = 0;
  if ((obj_type == OBJ_BUS) || (obj_type == OBJ_ALL)) { // Bus or all
    for (i = 0; i < net->num_buses; i++) {
      if (BUS_has_properties(NET_get_bus(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_GEN) || (obj_type == OBJ_ALL)) { // Gen or all
    for (i = 0; i < net->num_gens; i++) {
      if (GEN_has_properties(NET_get_gen(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_LOAD) || (obj_type == OBJ_ALL)) { // Load or all
    for (i = 0; i < net->num_loads; i++) {
      if (LOAD_has_properties(NET_get_load(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_BRANCH) || (obj_type == OBJ_ALL)) { // Branch or all
    for (i = 0; i < net->num_branches; i++) {
      if (BRANCH_has_properties(NET_get_branch(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_SHUNT) || (obj_type == OBJ_ALL)) { // Shunt or all
    for (i = 0; i < net->num_shunts; i++) {
      if (SHUNT_has_properties(NET_get_shunt(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_VARGEN) || (obj_type == OBJ_ALL)) { // Vargen or all
    for (i = 0; i < net->num_vargens; i++) {
      if (VARGEN_has_properties(NET_get_vargen(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_BAT) || (obj_type == OBJ_ALL)) { // Battery or all
    for (i = 0; i < net->num_bats; i++) {
      if (BAT_has_properties(NET_get_bat(net,i),prop_mask)) {
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
  }
  if ((obj_type == OBJ_CONVCSC) || (obj_type == OBJ_ALL)) { // Converter CSC or all
    for (i = 0; i < net->num_csc_convs; i++) {
      if (CONVCSC_has_properties(NET_get_csc_conv(net,i),prop_mask)) {
        indices = CONVCSC_get_var_indices(NET_get_csc_conv(net,i),var,t_start,t_end);
        for (j = 0; j < VEC_get_size(indices); j++) {
          MAT_set_i(proj,num_subvars,num_subvars);
          MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
          MAT_set_d(proj,num_subvars,1.);
          num_subvars++;
        }
        VEC_del(indices);
      }
    }
  }
  if ((obj_type == OBJ_CONVVSC) || (obj_type == OBJ_ALL)) { // Converter VSC or all
    for (i = 0; i < net->num_vsc_convs; i++) {
      if (CONVVSC_has_properties(NET_get_vsc_conv(net,i),prop_mask)) {
        indices = CONVVSC_get_var_indices(NET_get_vsc_conv(net,i),var,t_start,t_end);
        for (j = 0; j < VEC_get_size(indices); j++) {
          MAT_set_i(proj,num_subvars,num_subvars);
          MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
          MAT_set_d(proj,num_subvars,1.);
          num_subvars++;
        }
        VEC_del(indices);
      }
    }
  }
  if ((obj_type == OBJ_BUSDC) || (obj_type == OBJ_ALL)) { // DC bus or all
    for (i = 0; i < net->num_dc_buses; i++) {
      if (BUSDC_has_properties(NET_get_dc_bus(net,i),prop_mask)) {
        indices = BUSDC_get_var_indices(NET_get_dc_bus(net,i),var,t_start,t_end);
        for (j = 0; j < VEC_get_size(indices); j++) {
          MAT_set_i(proj,num_subvars,num_subvars);
          MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
          MAT_set_d(proj,num_subvars,1.);
          num_subvars++;
        }
        VEC_del(indices);
      }
    }
  }
  if ((obj_type == OBJ_BRANCHDC) || (obj_type == OBJ_ALL)) { // DC branch or all
    for (i = 0; i < net->num_dc_branches; i++) {
      if (BRANCHDC_has_properties(NET_get_dc_branch(net,i),prop_mask)) {
        indices = BRANCHDC_get_var_indices(NET_get_dc_branch(net,i),var,t_start,t_end);
        for (j = 0; j < VEC_get_size(indices); j++) {
          MAT_set_i(proj,num_subvars,num_subvars);
          MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
          MAT_set_d(proj,num_subvars,1.);
          num_subvars++;
        }
        VEC_del(indices);
      }
    }
  }
  if ((obj_type == OBJ_FACTS) || (obj_type == OBJ_ALL)) { // Facts or all
    for (i = 0; i < net->num_facts; i++) {
      if (FACTS_has_properties(NET_get_facts(net,i),prop_mask)) {
        indices = FACTS_get_var_indices(NET_get_facts(net,i),var,t_start,t_end);
        for (j = 0; j < VEC_get_size(indices); j++) {
          MAT_set_i(proj,num_subvars,num_subvars);
          MAT_set_j(proj,num_subvars,(int)VEC_get(indices,j));
          MAT_set_d(proj,num_subvars,1.);
          num_subvars++;
        }
        VEC_del(indices);
      }
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
  char* output_start;
  int max_size;

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
              BAT_BUFFER_SIZE*BAT_NUM_JSON_FIELDS*net->num_bats +
              CONVCSC_BUFFER_SIZE*CONVCSC_NUM_JSON_FIELDS*net->num_csc_convs +
              CONVVSC_BUFFER_SIZE*CONVVSC_NUM_JSON_FIELDS*net->num_vsc_convs +
              BUSDC_BUFFER_SIZE*BUSDC_NUM_JSON_FIELDS*net->num_dc_buses +
              BRANCHDC_BUFFER_SIZE*BRANCHDC_NUM_JSON_FIELDS*net->num_dc_branches +
              FACTS_BUFFER_SIZE*FACTS_NUM_JSON_FIELDS*net->num_facts)*net->num_periods;

  // Output
  output = (char*)malloc(sizeof(char)*max_size);
  output_start = output;
  
  // Write
  JSON_start(output);
  JSON_int(temp,output,"num_periods",net->num_periods,FALSE);
  JSON_float(temp,output,"base_power",net->base_power,FALSE);
  JSON_str(temp,output,"version",VERSION,FALSE);
  JSON_array_json(temp,output,"buses",net->bus,BUS_array_get,net->num_buses,BUS_get_json_string,FALSE);
  JSON_array_json(temp,output,"branches",net->branch,BRANCH_array_get,net->num_branches,BRANCH_get_json_string,FALSE);
  JSON_array_json(temp,output,"generators",net->gen,GEN_array_get,net->num_gens,GEN_get_json_string,FALSE);
  JSON_array_json(temp,output,"loads",net->load,LOAD_array_get,net->num_loads,LOAD_get_json_string,FALSE);
  JSON_array_json(temp,output,"shunts",net->shunt,SHUNT_array_get,net->num_shunts,SHUNT_get_json_string,FALSE);
  JSON_array_json(temp,output,"var_generators",net->vargen,VARGEN_array_get,net->num_vargens,VARGEN_get_json_string,FALSE);
  JSON_array_json(temp,output,"batteries",net->bat,BAT_array_get,net->num_bats,BAT_get_json_string,FALSE);
  JSON_array_json(temp,output,"csc_converters",net->csc_conv,CONVCSC_array_get,net->num_csc_convs,CONVCSC_get_json_string,FALSE);
  JSON_array_json(temp,output,"vsc_converters",net->vsc_conv,CONVVSC_array_get,net->num_vsc_convs,CONVVSC_get_json_string,FALSE);
  JSON_array_json(temp,output,"dc_buses",net->dc_bus,BUSDC_array_get,net->num_dc_buses,BUSDC_get_json_string,FALSE);
  JSON_array_json(temp,output,"dc_branches",net->dc_branch,BRANCHDC_array_get,net->num_dc_branches,BRANCHDC_get_json_string,FALSE);
  JSON_array_json(temp,output,"facts",net->facts,FACTS_array_get,net->num_facts,FACTS_get_json_string,TRUE);

  JSON_end(output);
  
  // Resize
  output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

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
  int i;
  if (net) {

    // Pointers
    net->branch = branch;
    net->num_branches = num;

    // Network
    for (i = 0; i < net->num_branches; i++)
      BRANCH_set_network(BRANCH_array_get(net->branch,i), net);
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
  int i;
  if (net) {

    // Pointers
    net->bus = bus;
    net->num_buses = num;

    // Network
    for (i = 0; i < net->num_buses; i++)
      BUS_set_network(BUS_array_get(net->bus,i), net);
  }
}

void NET_set_gen_array(Net* net, Gen* gen, int num) {
  int i;
  if (net) {

    // Pointers
    net->gen = gen;
    net->num_gens = num;

    // Network
    for (i = 0; i < net->num_gens; i++)
      GEN_set_network(GEN_array_get(net->gen,i), net);
  }
}

void NET_set_vargen_array(Net* net, Vargen* gen, int num) {
  if (net) {
    net->vargen = gen;
    net->num_vargens = num;
  }
}

void NET_set_bat_array(Net* net, Bat* bat, int num) {
  if (net) {
    net->bat = bat;
    net->num_bats = num;
  }
}

void NET_set_csc_conv_array(Net* net, ConvCSC* conv, int num) {
  if (net) {
    net->csc_conv = conv;
    net->num_csc_convs = num;
  }
}

void NET_set_vsc_conv_array(Net* net, ConvVSC* conv, int num) {
  if (net) {
    net->vsc_conv = conv;
    net->num_vsc_convs = num;
  }
}

void NET_set_dc_bus_array(Net* net, BusDC* bus, int num) {
  int i;
  if (net) {

    // Pointers
    net->dc_bus = bus;
    net->num_dc_buses = num;

    // Network
    for (i = 0; i < net->num_dc_buses; i++)
      BUSDC_set_network(BUSDC_array_get(net->dc_bus,i), net);
  }
}

void NET_set_dc_branch_array(Net* net, BranchDC* branch, int num) {
  if (net) {
    net->dc_branch = branch;
    net->num_dc_branches = num;
  }
}

void NET_set_facts_array(Net* net, Facts* facts, int num) {
  if (net) {
    net->facts = facts;
    net->num_facts = num;
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

  // Connect
  i = 0;
  bus = bus_list;
  while (i < net->num_vargens && bus) {
    gen = VARGEN_array_get(net->vargen,i);
    VARGEN_set_bus(gen,bus); // also removes old connection
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

  // Connect
  i = 0;
  bus = bus_list;
  while (i < net->num_bats && bus) {
    bat = BAT_array_get(net->bat,i);
    BAT_set_bus(bat,bus);    // also removes old connection
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
  case OBJ_CONVCSC:
    num = net->num_csc_convs;
    array = net->csc_conv;
    get_element = &CONVCSC_array_get;
    set_flags = &CONVCSC_set_flags;
    has_properties = &CONVCSC_has_properties;
    break;
  case OBJ_CONVVSC:
    num = net->num_vsc_convs;
    array = net->vsc_conv;
    get_element = &CONVVSC_array_get;
    set_flags = &CONVVSC_set_flags;
    has_properties = &CONVVSC_has_properties;
    break;
  case OBJ_BUSDC:
    num = net->num_dc_buses;
    array = net->dc_bus;
    get_element = &BUSDC_array_get;
    set_flags = &BUSDC_set_flags;
    has_properties = &BUSDC_has_properties;
    break;
  case OBJ_BRANCHDC:
    num = net->num_dc_branches;
    array = net->dc_branch;
    get_element = &BRANCHDC_array_get;
    set_flags = &BRANCHDC_set_flags;
    has_properties = &BRANCHDC_has_properties;
    break;
  case OBJ_FACTS:
    num = net->num_facts;
    array = net->facts;
    get_element = &FACTS_array_get;
    set_flags = &FACTS_set_flags;
    has_properties = &FACTS_has_properties;
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
  case OBJ_CONVCSC:
    set_flags = &CONVCSC_set_flags;
    get_obj_type = &CONVCSC_get_obj_type;
    break;
  case OBJ_CONVVSC:
    set_flags = &CONVVSC_set_flags;
    get_obj_type = &CONVVSC_get_obj_type;
    break;
  case OBJ_BUSDC:
    set_flags = &BUSDC_set_flags;
    get_obj_type = &BUSDC_get_obj_type;
    break;
  case OBJ_BRANCHDC:
    set_flags = &BRANCHDC_set_flags;
    get_obj_type = &BRANCHDC_get_obj_type;
    break;
  case OBJ_FACTS:
    set_flags = &FACTS_set_flags;
    get_obj_type = &FACTS_get_obj_type;
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

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++)
    CONVCSC_set_var_values(CONVCSC_array_get(net->csc_conv,i),values);

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++)
    CONVVSC_set_var_values(CONVVSC_array_get(net->vsc_conv,i),values);

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++)
    BUSDC_set_var_values(BUSDC_array_get(net->dc_bus,i),values);

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++)
    BRANCHDC_set_var_values(BRANCHDC_array_get(net->dc_branch,i),values);

  // Facts
  for (i = 0; i < net->num_facts; i++)
    FACTS_set_var_values(FACTS_array_get(net->facts,i),values);
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
  sprintf(out+strlen(out),"  star           : %d\n",NET_get_num_star_buses(net));
  sprintf(out+strlen(out),"shunts           : %d\n",NET_get_num_shunts(net));
  sprintf(out+strlen(out),"  fixed          : %d\n",NET_get_num_fixed_shunts(net));
  sprintf(out+strlen(out),"  switched       : %d\n",NET_get_num_switched_shunts(net));
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
  sprintf(out+strlen(out),"facts            : %d\n",NET_get_num_facts(net));
  sprintf(out+strlen(out),"csc converters   : %d\n",NET_get_num_csc_convs(net));
  sprintf(out+strlen(out),"vsc converters   : %d\n",NET_get_num_vsc_convs(net));
  sprintf(out+strlen(out),"  P dc mode      : %d\n",NET_get_num_vsc_convs_in_P_dc_mode(net));
  sprintf(out+strlen(out),"  v dc mode      : %d\n",NET_get_num_vsc_convs_in_v_dc_mode(net));
  sprintf(out+strlen(out),"  v ac mode      : %d\n",NET_get_num_vsc_convs_in_v_ac_mode(net));
  sprintf(out+strlen(out),"  f ac mode      : %d\n",NET_get_num_vsc_convs_in_f_ac_mode(net));
  sprintf(out+strlen(out),"dc buses         : %d\n",NET_get_num_dc_buses(net));
  sprintf(out+strlen(out),"dc branches      : %d\n",NET_get_num_dc_branches(net));

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
      case BUS_SENS_V_SET_REG:
        printf("%15s ","  v_set_reg ");
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
    for (i = 0; i < NET_get_num_buses(net); i++)
      NET_update_properties_step(net,NET_get_bus(net,i),NULL,t,values);
    for (i = 0; i < NET_get_num_dc_buses(net); i++)
      NET_update_properties_step(net,NULL,NET_get_dc_bus(net,i),t,values);
  }
}

void NET_update_properties_step(Net* net, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Branch* br;
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Shunt* shunt;
  Bat* bat;
  ConvCSC* csc_conv;
  ConvVSC* vsc_conv;
  Facts* facts;

  REAL P;
  REAL Q;
  REAL dQ;
  REAL dP;

  REAL v;
  REAL dv;

  REAL a;
  REAL da;
  REAL phi;
  REAL dphi;

  REAL shunt_b;
  REAL shunt_db;
  REAL shunt_g;

  // Check pointers
  if (!net || !bus)
    return;

  // Voltage magnitude
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && var_values)
    v = VEC_get(var_values,BUS_get_index_v_mag(bus,t));
  else
    v = BUS_get_v_mag(bus,t);

  // Maximum and minimum voltage magnitudes
  //***************************************
  if (net->bus_v_max[t] == 0 && net->bus_v_min[t] == 0) {
    net->bus_v_max[t] = v;
    net->bus_v_min[t] = v;
  }
  else {
    if (v > net->bus_v_max[t])
      net->bus_v_max[t] = v;
    if (v < net->bus_v_min[t])
      net->bus_v_min[t] = v;
  }

  // Normal voltage magnitude limit violations
  //******************************************
  dv = 0;
  if (v > BUS_get_v_max_norm(bus))
    dv = v-BUS_get_v_max_norm(bus);
  if (v < BUS_get_v_min_norm(bus))
    dv = BUS_get_v_min_norm(bus)-v;
  if (dv > net->bus_v_vio[t])
    net->bus_v_vio[t] = dv;

  // Regulation voltage magntiude limit violations
  //**********************************************
  dv = 0;
  if (v > BUS_get_v_max_reg(bus))
    dv = v-BUS_get_v_max_reg(bus);
  if (v < BUS_get_v_min_reg(bus))
    dv = BUS_get_v_min_reg(bus)-v;
  if (BUS_is_regulated_by_tran(bus)) { // false if all reg trans are on outage
    if (dv > net->tran_v_vio[t])
      net->tran_v_vio[t] = dv;
  }
  if (BUS_is_regulated_by_shunt(bus)) {
    if (dv > net->shunt_v_vio[t])
      net->shunt_v_vio[t] = dv;
  }

  // Bus regulated by gen
  if (BUS_is_regulated_by_gen(bus)) { // false if all reg gens are on outage

    // Voltage set point deviation
    //****************************
    if (fabs(v-BUS_get_v_set(bus,t)) > net->gen_v_dev[t])
      net->gen_v_dev[t] = fabs(v-BUS_get_v_set(bus,t));
    
    // Voltage set point action
    //*************************
    dv = BUS_get_v_max_reg(bus)-BUS_get_v_min_reg(bus);
    if (dv < NET_CONTROL_EPS)
      dv = NET_CONTROL_EPS;
    if (100.*fabs(v-BUS_get_v_set(bus,t))/dv > NET_CONTROL_ACTION_PCT)
      net->num_actions[t]++;
  }

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

    // Skip if gen on outage
    if (GEN_is_on_outage(gen))
      continue;
    
    // P Q
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
    
    // Reacive power of regulator
    if (GEN_is_regulator(gen) && !GEN_is_slack(gen)) {
      
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
    BUS_inject_P(bus,-shunt_g*v*v,t);
    BUS_inject_Q(bus,shunt_b*v*v,t);
    
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

  // VSC converters
  for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv)) {
    
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P) && var_values)
      P = VEC_get(var_values,CONVVSC_get_index_P(vsc_conv,t));
    else
      P = CONVVSC_get_P(vsc_conv,t);
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q) && var_values)
      Q = VEC_get(var_values,CONVVSC_get_index_Q(vsc_conv,t));
    else
      Q = CONVVSC_get_Q(vsc_conv,t);
    
    // Injections
    BUS_inject_P(bus,P,t);
    BUS_inject_Q(bus,Q,t);
  }

  // CSC converters
  for (csc_conv = BUS_get_csc_conv(bus); csc_conv != NULL; csc_conv = CONVCSC_get_next_ac(csc_conv)) {
    
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_P) && var_values)
      P = VEC_get(var_values,CONVCSC_get_index_P(csc_conv,t));
    else
      P = CONVCSC_get_P(csc_conv,t);
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_Q) && var_values)
      Q = VEC_get(var_values,CONVCSC_get_index_Q(csc_conv,t));
    else
      Q = CONVCSC_get_Q(csc_conv,t);
    
    // Injections
    BUS_inject_P(bus,P,t);
    BUS_inject_Q(bus,Q,t);
  }

  //FACTS
  for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {
    
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P) && var_values)
      P = VEC_get(var_values,FACTS_get_index_P_k(facts,t));
    else
      P = FACTS_get_P_k(facts,t);
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_Q) && var_values)
      Q = VEC_get(var_values,FACTS_get_index_Q_k(facts,t));
    else
      Q = FACTS_get_Q_k(facts,t);
    
    // Injections
    BUS_inject_P(bus,P,t);
    BUS_inject_Q(bus,Q,t);
  }
  for (facts = BUS_get_facts_m(bus); facts != NULL; facts = FACTS_get_next_m(facts)) {
    
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P) && var_values)
      P = VEC_get(var_values,FACTS_get_index_P_m(facts,t));
    else
      P = FACTS_get_P_m(facts,t);
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_Q) && var_values)
      Q = VEC_get(var_values,FACTS_get_index_Q_m(facts,t));
    else
      Q = FACTS_get_Q_m(facts,t);

    // Injections
    BUS_inject_P(bus,P,t);
    BUS_inject_Q(bus,Q,t);
  }
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {
  
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
    if (BRANCH_is_tap_changer(br) && !BRANCH_is_on_outage(br)) {

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
    if (BRANCH_is_phase_shifter(br) && !BRANCH_is_on_outage(br)) {
      
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
    if (!BRANCH_is_on_outage(br)) {
      BUS_inject_P(BRANCH_get_bus_k(br),-BRANCH_get_P_km(br,var_values,t),t);
      BUS_inject_Q(BRANCH_get_bus_k(br),-BRANCH_get_Q_km(br,var_values,t),t);
      BUS_inject_P(BRANCH_get_bus_m(br),-BRANCH_get_P_mk(br,var_values,t),t);
      BUS_inject_Q(BRANCH_get_bus_m(br),-BRANCH_get_Q_mk(br,var_values,t),t);
    }
  }
    
  // Power mismatches
  if (BUS_get_index(bus) == net->num_buses-1) {
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

void NET_update_hash_tables(Net* net) {

  // Local variables
  Bus* bus;
  BusDC* dc_bus;
  int i;

  // Check
  if (!net)
    return;

  // Clear
  BUS_hash_number_del(net->bus_hash_number);
  BUS_hash_name_del(net->bus_hash_name);
  BUSDC_hash_number_del(net->dc_bus_hash_number);
  BUSDC_hash_name_del(net->dc_bus_hash_name);
  net->bus_hash_number = NULL;
  net->bus_hash_name = NULL;
  net->dc_bus_hash_number = NULL;
  net->dc_bus_hash_name = NULL;

  // Update
  for (i = 0; i < net->num_buses; i++) {
    bus = NET_get_bus(net,i);
    NET_bus_hash_number_add(net,bus);
    NET_bus_hash_name_add(net,bus);
  }
  for (i = 0; i < net->num_dc_buses; i++) {
    dc_bus = NET_get_dc_bus(net,i);
    NET_dc_bus_hash_number_add(net,dc_bus);
    NET_dc_bus_hash_name_add(net,dc_bus);
  }
}

void NET_propagate_data_in_time(Net* net, int start, int end) {

  // Local variables
  int i;

  // No net
  if (!net)
    return;

  // Buses
  for (i = 0; i < net->num_buses; i++)
    BUS_propagate_data_in_time(BUS_array_get(net->bus,i),start,end);

  // Branches
  for (i = 0; i < net->num_branches; i++)
    BRANCH_propagate_data_in_time(BRANCH_array_get(net->branch,i),start,end);
  
  // Generators
  for (i = 0; i < net->num_gens; i++)
    GEN_propagate_data_in_time(GEN_array_get(net->gen,i),start,end);

  // Loads
  for (i = 0; i < net->num_loads; i++)
    LOAD_propagate_data_in_time(LOAD_array_get(net->load,i),start,end);

  // Vargens
  for (i = 0; i < net->num_vargens; i++)
    VARGEN_propagate_data_in_time(VARGEN_array_get(net->vargen,i),start,end);

  // Shunts
  for (i = 0; i < net->num_shunts; i++)
    SHUNT_propagate_data_in_time(SHUNT_array_get(net->shunt,i),start,end);

  // Batteries
  for (i = 0; i < net->num_bats; i++)
    BAT_propagate_data_in_time(BAT_array_get(net->bat,i),start,end);

  // CSC converters
  for (i = 0; i < net->num_csc_convs; i++)
    CONVCSC_propagate_data_in_time(CONVCSC_array_get(net->csc_conv,i),start,end);

  // VSC converters
  for (i = 0; i < net->num_vsc_convs; i++)
    CONVVSC_propagate_data_in_time(CONVVSC_array_get(net->vsc_conv,i),start,end);

  // DC buses
  for (i = 0; i < net->num_dc_buses; i++)
    BUSDC_propagate_data_in_time(BUSDC_array_get(net->dc_bus,i),start,end);

  // DC branches
  for (i = 0; i < net->num_dc_branches; i++)
    BRANCHDC_propagate_data_in_time(BRANCHDC_array_get(net->dc_branch,i),start,end);

  // Facts
  for (i = 0; i < net->num_facts; i++)
    FACTS_propagate_data_in_time(FACTS_array_get(net->facts,i),start,end);
}
