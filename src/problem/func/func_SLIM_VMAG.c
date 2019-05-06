/** @file func_SLIM_VMAG.c
 *  @brief This file defines the data structure and routines associated with the function of type SLIM_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_SLIM_VMAG.h>

Func* FUNC_SLIM_VMAG_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_SLIM_VMAG_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_SLIM_VMAG_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_SLIM_VMAG_eval_step);
  FUNC_set_name(f,"soft voltage magnitude limits");
  return f;
}

void FUNC_SLIM_VMAG_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) // v var
    (*Hphi_nnz)++;
}

void FUNC_SLIM_VMAG_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Mat* Hphi;
  REAL dv;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  dv = BUS_get_v_max_norm(bus)-BUS_get_v_min_norm(bus);
  if (dv < FUNC_SLIM_VMAG_PARAM)
    dv = FUNC_SLIM_VMAG_PARAM;

  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
    MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
    MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
    (*Hphi_nnz)++;
  }
}

void FUNC_SLIM_VMAG_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_v_mag;
  REAL v;
  REAL vmid;
  REAL dv;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  dv = BUS_get_v_max_norm(bus)-BUS_get_v_min_norm(bus);
  if (dv < FUNC_SLIM_VMAG_PARAM)
    dv = FUNC_SLIM_VMAG_PARAM;
  
  vmid = 0.5*(BUS_get_v_max_norm(bus)+BUS_get_v_min_norm(bus));
  
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
    
    // Index
    index_v_mag = BUS_get_index_v_mag(bus,t);
    
    // v
    v = VEC_get(var_values,index_v_mag);
    
    // phi
    (*phi) += 0.5*pow((v-vmid)/dv,2.);
    
    // gphi
    gphi[index_v_mag] = (v-vmid)/(dv*dv);
    
    // Hphi
    Hphi[*Hphi_nnz] = 1./(dv*dv);
    (*Hphi_nnz)++;	
  }
  else{
    
    // v
    v = BUS_get_v_mag(bus,t);
    
    // phi
    (*phi) += 0.5*pow((v-vmid)/dv,2.);
  }
}
