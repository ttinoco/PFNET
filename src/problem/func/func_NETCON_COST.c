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

void FUNC_NETCON_COST_count_step(Func* f, Branch* br, int t) {
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

void FUNC_NETCON_COST_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  int bus_index[2];
  char* bus_counted;
  REAL price;
  Vec* gphi;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  gphi = FUNC_get_gphi(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!gphi || !bus_counted)
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
    
    price = BUS_get_price(bus,t);
    
    if (!bus_counted[bus_index[k]*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  VEC_set(gphi,GEN_get_index_P(gen,t),-price);
      }
      
      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
	  VEC_set(gphi,VARGEN_get_index_P(vargen,t),-price);
      }
      
      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  VEC_set(gphi,LOAD_get_index_P(load,t),price);
      }
      
      // Battery charging
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
	  VEC_set(gphi,BAT_get_index_Pc(bat,t),price);
	  VEC_set(gphi,BAT_get_index_Pd(bat,t),-price);
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]*T+t] = TRUE;
  }
}

void FUNC_NETCON_COST_eval_step(Func* f, Branch* br, int t, Vec* var_values) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  int bus_index[2];
  char* bus_counted;
  REAL price;
  REAL* phi;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !bus_counted)
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
    
    price = BUS_get_price(bus,t);
    
    if (!bus_counted[bus_index[k]*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  (*phi) -= price*VEC_get(var_values,GEN_get_index_P(gen,t));
	else
	  (*phi) -= price*GEN_get_P(gen,t);
      }
      
      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
	  (*phi) -= price*VEC_get(var_values,VARGEN_get_index_P(vargen,t));
	else
	  (*phi) -= price*VARGEN_get_P(vargen,t);
      }
      
      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  (*phi) += price*VEC_get(var_values,LOAD_get_index_P(load,t));
	else
	  (*phi) += price*LOAD_get_P(load,t);
      }
      
      // Battery charging
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
	  (*phi) += price*VEC_get(var_values,BAT_get_index_Pc(bat,t));
	  (*phi) -= price*VEC_get(var_values,BAT_get_index_Pd(bat,t));
	}
	else {
	  (*phi) += price*BAT_get_P(bat,t);
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]*T+t] = TRUE;
  }  
}

void FUNC_NETCON_COST_free(Func* f) {
  // Nothing
}
