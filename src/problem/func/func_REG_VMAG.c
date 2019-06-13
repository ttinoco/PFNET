/** @file func_REG_VMAG.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_VMAG.h>

struct Func_REG_VMAG_Data {
  BOOL v_set_ref;
};

Func* FUNC_REG_VMAG_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_init(f,&FUNC_REG_VMAG_init);
  FUNC_set_func_count_step(f,&FUNC_REG_VMAG_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_REG_VMAG_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_VMAG_eval_step);
  FUNC_set_func_free(f,&FUNC_REG_VMAG_free);
  FUNC_set_func_set_parameter(f,&FUNC_REG_VMAG_set_parameter);
  FUNC_set_name(f,"voltage magnitude regularization");
  FUNC_init(f);
  return f;
}

void FUNC_REG_VMAG_init(Func* f) {

  // Local variables
  Func_REG_VMAG_Data* data;
  
  // Init
  data = (Func_REG_VMAG_Data*)malloc(sizeof(Func_REG_VMAG_Data));
  data->v_set_ref = TRUE;
  FUNC_set_data(f,data);
}

void FUNC_REG_VMAG_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;

  // Func data
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

void FUNC_REG_VMAG_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Mat* Hphi;

  // Func data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
    MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
    MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
    (*Hphi_nnz)++;
  }
}

void FUNC_REG_VMAG_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_v_mag;
  REAL v;
  REAL vref;
  REAL dv = FUNC_REG_VMAG_PARAM;
  Func_REG_VMAG_Data* data;

  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  data = (Func_REG_VMAG_Data*)FUNC_get_data(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus || !data)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // Reference
  if (data->v_set_ref || BUS_is_v_set_regulated(bus,TRUE))
    vref = BUS_get_v_set(bus,t);
  else
    vref = BUS_get_v_mag(bus,t);
  
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
    
    // Index
    index_v_mag = BUS_get_index_v_mag(bus,t);
    
    // v
    v = VEC_get(var_values,index_v_mag);

    // phi
    (*phi) += 0.5*pow((v-vref)/dv,2.);
    
    // gphi
    gphi[index_v_mag] = (v-vref)/(dv*dv);
    
    // Hphi
    Hphi[*Hphi_nnz] = 1./(dv*dv);
    (*Hphi_nnz)++;
  }
  else {
    
    // v
    v = BUS_get_v_mag(bus,t);
    
    // phi
    (*phi) += 0.5*pow((v-vref)/dv,2.);
  }
}

void FUNC_REG_VMAG_free(Func* f) {

  // Local variables
  Func_REG_VMAG_Data* data;

  // Get data
  data = (Func_REG_VMAG_Data*)FUNC_get_data(f);

  // Free
  if (data)
    free(data);

  // Set data
  FUNC_set_data(f,NULL);
}

void FUNC_REG_VMAG_set_parameter(Func* f, char* key, void* value) {

  // Local variables
  Func_REG_VMAG_Data* data = (Func_REG_VMAG_Data*)FUNC_get_data(f);

  // Check
  if (!data)
    return;

  // Set 
  if (strcmp(key,"v_set_reference") == 0) 
    data->v_set_ref = *(BOOL*)value;
  else // unknown
    FUNC_set_error(f,"invalid parameter");
}
