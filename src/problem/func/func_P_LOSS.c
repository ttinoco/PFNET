/** @file func_P_LOSS.c
 *  @brief This file defines the data structure and routines associated with the function of type P_LOSS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_P_LOSS.h>

Func* FUNC_P_LOSS_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_P_LOSS_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_P_LOSS_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_P_LOSS_eval_step);
  FUNC_set_name(f,"active power loss");
  return f;
}

void FUNC_P_LOSS_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

    // pass No Hessian
}

void FUNC_P_LOSS_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {
  
    // pass No Hessian
}

void FUNC_P_LOSS_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  Load* load;
  REAL* phi;
  REAL* gphi;

  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));

  // Check pointers
  if (!phi || !gphi || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

    // Out of service
    if (!GEN_is_in_service(gen))
      continue;

    if (GEN_has_flags(gen, FLAG_VARS, GEN_VAR_P)) { // Pg var

      // phi
      (*phi) += VEC_get(var_values, GEN_get_index_P(gen, t));

      // gphi
      gphi[GEN_get_index_P(gen, t)] = 1.0;
    }
    else {

      // phi
      (*phi) += GEN_get_P(gen, t);
    }
  }

  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)){

    // Out of service
    if (!VARGEN_is_in_service(vargen))
      continue;
    
    // Pg var
    if (VARGEN_has_flags(vargen, FLAG_VARS, VARGEN_VAR_P)){ 

      // phi
      (*phi) += VEC_get(var_values, VARGEN_get_index_P(vargen, t));

      // gphi
      gphi[VARGEN_get_index_P(vargen, t)] = 1.0;
    }
    else {

      // phi
      (*phi) += VARGEN_get_P(vargen, t);
    }
  }

  // Battery 
  for(bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)){

    // Out of service 
    if (!BAT_is_in_service(bat))
      continue;
    
    // P_bat is var
    if(BAT_has_flags(bat, FLAG_VARS, BAT_VAR_P)){

      // phi 
      (*phi) -= VEC_get(var_values, BAT_get_index_Pc(bat, t));   // charging
      (*phi) += VEC_get(var_values, BAT_get_index_Pd(bat, t));   // discharging

      // gphi
      gphi[BAT_get_index_Pc(bat, t)] = -1;
      gphi[BAT_get_index_Pd(bat, t)] =  1;

    }
    else{

      //phi
      (*phi) -= BAT_get_P(bat, t);
    }
  }

  // Loads
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;

    // Pl Variable
    if (LOAD_has_flags(load, FLAG_VARS, LOAD_VAR_P)) {

      // phi
      (*phi) -= VEC_get(var_values, LOAD_get_index_P(load, t));

      // gphi
      gphi[LOAD_get_index_P(load, t)] = -1.0;

    }

    // Pl Constant
    else {

      // phi
      (*phi) -= LOAD_get_P(load, t);
    }
  }
}
