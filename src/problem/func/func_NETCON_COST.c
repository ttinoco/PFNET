/** @file func_NETCON_COST.c
 *  @brief This file defines the data structure and routines associated with the function of type NETCON_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_NETCON_COST.h>

void FUNC_NETCON_COST_init(Func* f) {
  // Nothing
}

void FUNC_NETCON_COST_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);
  
  // gphi
  // Constant
  
  // Hphi
  // Zero
    
  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_NETCON_COST_count_branch(Func* f, Branch* br) {
  // nothing
}

void FUNC_NETCON_COST_allocate(Func* f) {
  
  // Local variables
  int num_vars;
  
  num_vars = NET_get_num_vars(FUNC_get_network(f));

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));
  VEC_set_zero(FUNC_get_gphi(f));
  
  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  0));
}

void FUNC_NETCON_COST_analyze_branch(Func* f, Branch* br) {

  /*
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  int bus_index[2];
  char* bus_counted;
  Vec* g;
  int k;

  // Constr data
  g = FUNC_get_gphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  MAT_set_i(H,*Hcounter,LOAD_get_index_P(load));
	  MAT_set_j(H,*Hcounter,LOAD_get_index_P(load));
	  MAT_set_d(H,*Hcounter,2.*LOAD_get_util_coeff_Q2(load));
	  (*Hcounter)++;
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }

  */
}

void FUNC_NETCON_COST_eval_branch(Func* f, Branch* br, Vec* var_values) {
  // need to compute phi (nothing else)
}

void FUNC_NETCON_COST_free(Func* f) {
  // Nothing
}
