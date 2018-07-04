/** @file func_NETCON_COST.c
 *  @brief This file defines the data structure and routines associated with the function of type NETCON_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_NETCON_COST.h>

Func* FUNC_NETCON_COST_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f, &FUNC_NETCON_COST_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_NETCON_COST_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_NETCON_COST_eval_step);
  FUNC_set_name(f,"net consumption cost");
  return f;
}

void FUNC_NETCON_COST_count_step(Func* f, Bus* bus, int t) {
  // nothing
}

void FUNC_NETCON_COST_analyze_step(Func* f, Bus* bus, int t) {
  // nothing
}

void FUNC_NETCON_COST_eval_step(Func* f, Bus* bus, int t, Vec* var_values) {

  // Local variables
  Load* load;
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  REAL price;
  REAL* phi;
  REAL* gphi;
  int k;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));

  // Check pointers
  if (!phi || !gphi)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
      (*phi) -= price*VEC_get(var_values,GEN_get_index_P(gen,t));
      gphi[GEN_get_index_P(gen,t)] = -price;
    }
    else
      (*phi) -= price*GEN_get_P(gen,t);
  }

  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
      (*phi) -= price*VEC_get(var_values,VARGEN_get_index_P(vargen,t));
      gphi[VARGEN_get_index_P(vargen,t)] = -price;
    }
    else
      (*phi) -= price*VARGEN_get_P(vargen,t);
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
      (*phi) += price*VEC_get(var_values,LOAD_get_index_P(load,t));
      gphi[LOAD_get_index_P(load,t)] = price;
    }
    else
      (*phi) += price*LOAD_get_P(load,t);
  }
  
  // Battery charging
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
      (*phi) += price*VEC_get(var_values,BAT_get_index_Pc(bat,t));
      (*phi) -= price*VEC_get(var_values,BAT_get_index_Pd(bat,t));
      gphi[BAT_get_index_Pc(bat,t)] = price;
      gphi[BAT_get_index_Pd(bat,t)] = -price;
    }
    else {
      (*phi) += price*BAT_get_P(bat,t);
    }
  }
}
