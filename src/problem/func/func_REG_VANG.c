/** @file func_REG_VANG.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VANG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_VANG.h>

Func* FUNC_REG_VANG_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_REG_VANG_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_REG_VANG_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_VANG_eval_step);
  FUNC_set_name(f,"voltage angle regularization");
  return f;
}

void FUNC_REG_VANG_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;

  // Func data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Bus
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG))
    (*Hphi_nnz)++; // w and w
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
  
    if (BUS_has_flags(BRANCH_get_bus_k(br),FLAG_VARS,BUS_VAR_VANG)) {
      (*Hphi_nnz)++; // wk and wk

      if (BUS_has_flags(BRANCH_get_bus_m(br),FLAG_VARS,BUS_VAR_VANG))
        (*Hphi_nnz)++; // wk and wm
    }
    
    if (BUS_has_flags(BRANCH_get_bus_m(br),FLAG_VARS,BUS_VAR_VANG))
      (*Hphi_nnz)++; // wm and wm
  }
}

void FUNC_REG_VANG_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;
  Bus* bus_k;
  Bus* bus_m;
  Mat* Hphi;
  
  // Func data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Bus
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
    MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus,t));
    MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus,t));
    (*Hphi_nnz)++; // w and w
  }

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);
    
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) {
      MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_k,t));
      MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_k,t));
      (*Hphi_nnz)++; // wk and wk
      
      if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) {
        MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_k,t));
        MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_m,t));
        (*Hphi_nnz)++; // wk and wm
      }
    }
    
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) {
      MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_m,t));
      MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_ang(bus_m,t));
      (*Hphi_nnz)++; // wm and wm
    }
  }
}

void FUNC_REG_VANG_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Branch* br;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  Bus* bus_k;
  Bus* bus_m;
  REAL wdiff;
  REAL dw = FUNC_REG_VANG_PARAM;
  REAL wk;
  REAL wm;

  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  
  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;

  // Bus
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
    wk = VEC_get(var_values,BUS_get_index_v_ang(bus,t));
    gphi[BUS_get_index_v_ang(bus,t)] += wk/(dw*dw);
    Hphi[*Hphi_nnz] = 1./(dw*dw);
    (*Hphi_nnz)++; // w and w
  }    
  else
    wk = BUS_get_v_ang(bus,t);

  // phi
  (*phi) += 0.5*pow(wk/dw,2.);
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;

    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);

    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG))
      wm = VEC_get(var_values,BUS_get_index_v_ang(bus_m,t));
    else
      wm = BUS_get_v_ang(bus_m,t);
    
    // Difference
    wdiff = wk-wm-BRANCH_get_phase(br,t);
    
    // gphi
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) {
      gphi[BUS_get_index_v_ang(bus_k,t)] += wdiff/(dw*dw);
      Hphi[*Hphi_nnz] = 1./(dw*dw);
      (*Hphi_nnz)++; // wk and wk
      if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) {
        Hphi[*Hphi_nnz] = -1./(dw*dw);
        (*Hphi_nnz)++; // wk and wm
      }
    }
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) {
      gphi[BUS_get_index_v_ang(bus_m,t)] -= wdiff/(dw*dw);
      Hphi[*Hphi_nnz] = 1./(dw*dw);
      (*Hphi_nnz)++; // wm and wm
    }

    // phi
    (*phi) += 0.5*pow(wdiff/dw,2.);
  }
}
