/** @file func_GEN_COST.c
 *  @brief This file defines the data structure and routines associated with the function of type GEN_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_GEN_COST.h>

void FUNC_GEN_COST_init(Func* f) {
  // Nothing
}

void FUNC_GEN_COST_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);
  
  // gphi
  VEC_set_zero(FUNC_get_gphi(f));
  
  // Hphi
  // Constant so not clear it
  
  // Counter
  FUNC_set_Hcounter(f,0);
  
  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_GEN_COST_count_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  int k;
  
  // Constr data
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

    if (!bus_counted[bus_index[k]*BRANCH_get_num_periods(br)]) {
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  (*Hcounter)++;
      }     
    }
    
    // Update counted flag
    bus_counted[bus_index[k]*BRANCH_get_num_periods(br)] = TRUE;
  }
}

void FUNC_GEN_COST_allocate(Func* f) {
  
  // Local variables
  int num_vars;
  int Hcounter;
  Net* net;
  
  net = FUNC_get_network(f);
  num_vars = NET_get_num_vars(net);
  Hcounter = FUNC_get_Hcounter(f);

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));

  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  Hcounter*NET_get_num_periods(net)));
}

void FUNC_GEN_COST_analyze_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  int t;
  int T;

  // Constr data
  H = FUNC_get_Hphi(f);
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

  // Time loop
  T = BRANCH_get_num_periods(br);
  for (t = 0; t < T; t++) {

    // Buses
    for (k = 0; k < 2; k++) {
    
      bus = buses[k];
      
      if (!bus_counted[bus_index[k]*T+t]) {
	for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	  if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	    MAT_set_i(H,*Hcounter,GEN_get_index_P(gen,t));
	    MAT_set_j(H,*Hcounter,GEN_get_index_P(gen,t));
	    MAT_set_d(H,*Hcounter,2.*GEN_get_cost_coeff_Q2(gen));
	    (*Hcounter)++;
	  }
	}
      }
      
      // Update counted flag
      bus_counted[bus_index[k]*T+t] = TRUE;
    }
  }
}

void FUNC_GEN_COST_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_P;
  REAL P;
  REAL Q0;
  REAL Q1;
  REAL Q2;
  int k;
  int t;
  int T;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !gphi || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Time loop
  T = BRANCH_get_num_periods(br);
  for (t = 0; t < T; t++) {

    // Buses
    for (k = 0; k < 2; k++) {
      
      bus = buses[k];
      
      if (!bus_counted[bus_index[k]*T+t]) {
	
	for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	  
	  Q0 = GEN_get_cost_coeff_Q0(gen);
	  Q1 = GEN_get_cost_coeff_Q1(gen);
	  Q2 = GEN_get_cost_coeff_Q2(gen);
	  
	  // Variable
	  if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	    
	    // Index
	    index_P = GEN_get_index_P(gen,t);
	    
	    // P
	    P = VEC_get(var_values,index_P);
	    
	    // phi
	    (*phi) += Q0 + Q1*P + Q2*pow(P,2.);
	    
	    // gphi
	    gphi[index_P] = Q1 + 2.*Q2*P;
	  }

	  // Constant
	  else {
	    
	    // P
	    P = GEN_get_P(gen,t);
	    
	    // phi
	    (*phi) += Q0 + Q1*P + Q2*pow(P,2.);
	  }
	}
      }
      
      // Update counted flag
      bus_counted[bus_index[k]*T+t] = TRUE;
    }
  }
}

void FUNC_GEN_COST_free(Func* f) {
  // Nothing
}
