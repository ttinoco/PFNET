/** @file func_GEN_COST.c
 *  @brief This file defines the data structure and routines associated with the function of type GEN_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_GEN_COST.h>

Func* FUNC_GEN_COST_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_GEN_COST_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_GEN_COST_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_GEN_COST_eval_step);
  FUNC_set_name(f,"generation cost");
  return f;
}

void FUNC_GEN_COST_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

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
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) && !GEN_is_on_outage(gen))
      (*Hphi_nnz)++;
  }
}

void FUNC_GEN_COST_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

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
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) && !GEN_is_on_outage(gen)) {
      MAT_set_i(Hphi,*Hphi_nnz,GEN_get_index_P(gen,t));
      MAT_set_j(Hphi,*Hphi_nnz,GEN_get_index_P(gen,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_GEN_COST_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Gen* gen;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_P;
  REAL P;
  REAL Q0;
  REAL Q1;
  REAL Q2;

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

    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    // Cost coefficients
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
      
      // Hphi
      Hphi[*Hphi_nnz] = 2.*Q2;
      (*Hphi_nnz)++;
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
