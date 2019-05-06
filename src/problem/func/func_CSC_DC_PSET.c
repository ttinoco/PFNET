/** @file func_CSC_DC_PSET.c
 *  @brief This file defines the data structure and routines associated with the function of type CSC_DC_PSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_CSC_DC_PSET.h>

Func* FUNC_CSC_DC_PSET_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_CSC_DC_PSET_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_CSC_DC_PSET_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_CSC_DC_PSET_eval_step);
  FUNC_set_name(f,"CSC DC power control");
  return f;
}

void FUNC_CSC_DC_PSET_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvCSC* csc;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;
  
  // CSC converters
  for (csc = BUS_get_csc_conv(bus); csc != NULL; csc = CONVCSC_get_next_ac(csc)) {   
    if (CONVCSC_is_in_P_dc_mode(csc) &&
        CONVCSC_has_flags(csc,FLAG_VARS,CONVCSC_VAR_PDC) &&
        CONVCSC_is_in_service(csc))
      (*Hphi_nnz)++;
  }
}

void FUNC_CSC_DC_PSET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvCSC* csc;
  int* Hphi_nnz;
  Mat* H;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // CSC converters
  for (csc = BUS_get_csc_conv(bus); csc != NULL; csc = CONVCSC_get_next_ac(csc)) {
    if (CONVCSC_is_in_P_dc_mode(csc) &&
        CONVCSC_has_flags(csc,FLAG_VARS,CONVCSC_VAR_PDC) &&
        CONVCSC_is_in_service(csc)) {
      MAT_set_i(H,*Hphi_nnz,CONVCSC_get_index_P_dc(csc,t));
      MAT_set_j(H,*Hphi_nnz,CONVCSC_get_index_P_dc(csc,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_CSC_DC_PSET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  ConvCSC* csc;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL P_dc;
  REAL P_dc_set;
  int index_P_dc;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;
  
  // CSC converters
  for (csc = BUS_get_csc_conv(bus); csc != NULL; csc = CONVCSC_get_next_ac(csc)) {

    // No P_dc_mode
    if (!CONVCSC_is_in_P_dc_mode(csc))
      continue;

    // Out of service
    if (!CONVCSC_is_in_service(csc))
      continue;
    
    // Set point
    P_dc_set = CONVCSC_get_P_dc_set(csc,t);
    
    if (CONVCSC_has_flags(csc,FLAG_VARS,CONVCSC_VAR_PDC)) {

      // Index
      index_P_dc = CONVCSC_get_index_P_dc(csc,t);
      
      // Value
      P_dc = VEC_get(var_values,index_P_dc);
      
      // phi
      (*phi) += 0.5*pow(P_dc-P_dc_set,2.);
      
      // gphi
      gphi[index_P_dc] = P_dc-P_dc_set;

      // Hphi
      Hphi[*Hphi_nnz] = 1.;
      (*Hphi_nnz)++;
    }
    else {
      
      // Value
      P_dc = CONVCSC_get_P_dc(csc,t);
      
      // phi
      (*phi) += 0.5*pow(P_dc-P_dc_set,2.);
    }
  }
}
