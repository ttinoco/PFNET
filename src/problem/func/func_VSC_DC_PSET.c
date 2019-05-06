/** @file func_VSC_DC_PSET.c
 *  @brief This file defines the data structure and routines associated with the function of type VSC_DC_PSET.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_VSC_DC_PSET.h>

Func* FUNC_VSC_DC_PSET_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_VSC_DC_PSET_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_VSC_DC_PSET_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_VSC_DC_PSET_eval_step);
  FUNC_set_name(f,"VSC DC power control");
  return f;
}

void FUNC_VSC_DC_PSET_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;
  
  // VSC converters
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {   
    if (CONVVSC_is_in_P_dc_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_is_in_service(vsc))
      (*Hphi_nnz)++;
  }
}

void FUNC_VSC_DC_PSET_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
  int* Hphi_nnz;
  Mat* H;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // VSC converters
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {
    if (CONVVSC_is_in_P_dc_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_is_in_service(vsc)) {
      MAT_set_i(H,*Hphi_nnz,CONVVSC_get_index_P(vsc,t));
      MAT_set_j(H,*Hphi_nnz,CONVVSC_get_index_P(vsc,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_VSC_DC_PSET_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  ConvVSC* vsc;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL P_ac;
  REAL P_dc_set;
  int index_P;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;
  
  // VSC converters
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {

    // No P_dc_mode
    if (!CONVVSC_is_in_P_dc_mode(vsc))
      continue;

    // Out of service
    if (!CONVVSC_is_in_service(vsc))
      continue;
    
    // Set point
    P_dc_set = CONVVSC_get_P_dc_set(vsc,t);
    
    if (CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P)) {

      // Index
      index_P = CONVVSC_get_index_P(vsc,t);
      
      // Value
      P_ac = VEC_get(var_values,index_P);
      
      // phi
      (*phi) += 0.5*pow(P_ac+P_dc_set,2.);
      
      // gphi
      gphi[index_P] = P_ac+P_dc_set;

      // Hphi
      Hphi[*Hphi_nnz] = 1.;
      (*Hphi_nnz)++;
    }
    else {
      
      // Value
      P_ac = CONVVSC_get_P(vsc,t);
      
      // phi
      (*phi) += 0.5*pow(P_ac+P_dc_set,2.);
    }
  }
}
