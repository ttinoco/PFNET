/** @file func_GEN_RED.c
 *  @brief This file defines the data structure and routines associated with the function of type GEN_RED.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_GEN_RED.h>

Func* FUNC_GEN_RED_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_GEN_RED_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_GEN_RED_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_GEN_RED_eval_step);
  FUNC_set_name(f,"generation redispatch penalty");
  return f;
}

void FUNC_GEN_RED_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Gen* gen;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Out of service
    if (!GEN_is_in_service(gen))
      continue;
        
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) // P var
      (*Hphi_nnz)++;
  }
}

void FUNC_GEN_RED_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Gen* gen;
  int* Hphi_nnz;
  Mat* Hphi;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Out of service
    if (!GEN_is_in_service(gen))
      continue;
        
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
      
      MAT_set_i(Hphi,*Hphi_nnz,GEN_get_index_P(gen,t));
      MAT_set_j(Hphi,*Hphi_nnz,GEN_get_index_P(gen,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_GEN_RED_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Gen* gen;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL P;
  REAL P0;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Out of service
    if (!GEN_is_in_service(gen))
      continue;
    
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
      
      // Value
      P0 = GEN_get_P(gen,t);
      P = VEC_get(var_values,GEN_get_index_P(gen,t));
      
      // phi
      (*phi) += 0.5*pow(P-P0,2.);
      
      // gphi
      gphi[GEN_get_index_P(gen,t)] = P-P0;
      
      // Hphi
      Hphi[*Hphi_nnz] = 1.;
      (*Hphi_nnz)++;
    }
    else {
            
      // phi
      (*phi) += 0.;
    }
  }
}
