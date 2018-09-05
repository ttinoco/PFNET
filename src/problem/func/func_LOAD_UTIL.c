/** @file func_LOAD_UTIL.c
 *  @brief This file defines the data structure and routines associated with the function of type LOAD_UTIL.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_LOAD_UTIL.h>

Func* FUNC_LOAD_UTIL_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_LOAD_UTIL_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_LOAD_UTIL_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_LOAD_UTIL_eval_step);
  FUNC_set_name(f,"consumption utility");
  return f;
}

void FUNC_LOAD_UTIL_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Load* load;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
      (*Hphi_nnz)++;
  }
}

void FUNC_LOAD_UTIL_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Load* load;
  int* Hphi_nnz;
  Mat* Hphi;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
      MAT_set_i(Hphi,*Hphi_nnz,LOAD_get_index_P(load,t));
      MAT_set_j(Hphi,*Hphi_nnz,LOAD_get_index_P(load,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_LOAD_UTIL_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Load* load;
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

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    Q0 = LOAD_get_util_coeff_Q0(load);
    Q1 = LOAD_get_util_coeff_Q1(load);
    Q2 = LOAD_get_util_coeff_Q2(load);
    
    // Variable
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
      
      // Index
      index_P = LOAD_get_index_P(load,t);
      
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
      P = LOAD_get_P(load,t);
      
      // phi
      (*phi) += Q0 + Q1*P + Q2*pow(P,2.);
    }
  }
}
