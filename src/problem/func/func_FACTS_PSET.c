/** @file func_FACTS_PSET.c
 *  @brief This file defines the data structure and routines associated with the function of type FACTS_PSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_FACTS_PSET.h>

Func* FUNC_FACTS_PSET_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f, &FUNC_FACTS_PSET_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_FACTS_PSET_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_FACTS_PSET_eval_step);
  FUNC_set_name(f,"FACTS active power control");
  return f;
}

void FUNC_FACTS_PSET_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Facts* facts;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;
  
  // FACTS
  for (facts = BUS_get_facts_k(bus); facts !=NULL; facts = FACTS_get_next_k(facts)) {    
    if (FACTS_is_in_normal_series_mode(facts) &&
        FACTS_get_P_max_dc(facts) > 0. &&
        FACTS_is_in_service(facts) &&
        FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P))
      (*Hphi_nnz)++;
  }
}

void FUNC_FACTS_PSET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Facts* facts;
  int* Hphi_nnz;
  Mat* H;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;
      
  // FACTS
  for (facts = BUS_get_facts_k(bus); facts !=NULL; facts = FACTS_get_next_k(facts)) {
    if (FACTS_is_in_normal_series_mode(facts) &&
        FACTS_get_P_max_dc(facts) > 0. &&
        FACTS_is_in_service(facts) &&
        FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P)) {
      MAT_set_i(H,*Hphi_nnz,FACTS_get_index_P_m(facts,t));
      MAT_set_j(H,*Hphi_nnz,FACTS_get_index_P_m(facts,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_FACTS_PSET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Facts* facts;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL P;
  REAL P_set;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts !=NULL; facts = FACTS_get_next_k(facts)) {
    
    if (FACTS_is_in_normal_series_mode(facts) && FACTS_get_P_max_dc(facts) > 0. && FACTS_is_in_service(facts)) {
      
      // Set point
      P_set = FACTS_get_P_set(facts,t);
      
      if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P)) { // P var
        
        // Value
        P = VEC_get(var_values,FACTS_get_index_P_m(facts,t));
        
        // phi
        (*phi) += 0.5*pow(P-P_set,2.);
        
        // gphi
        gphi[FACTS_get_index_P_m(facts,t)] = P-P_set;

        // Hphi
        Hphi[*Hphi_nnz] = 1.;
        (*Hphi_nnz)++;
      }
      else {
        
        // Value
        P = FACTS_get_P_m(facts,t);
        
        // phi
        (*phi) += 0.5*pow(P-P_set,2.);
      }
    }
  }
}
