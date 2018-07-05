/** @file func_REG_RATIO.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_RATIO.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_RATIO.h>

Func* FUNC_REG_RATIO_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_REG_RATIO_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_REG_RATIO_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_RATIO_eval_step);
  FUNC_set_name(f,"tap ratio regularization");
  return f;
}

void FUNC_REG_RATIO_count_step(Func* f, Bus* bus, int t) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      return;
  
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) // ratio var
      (*Hphi_nnz)++;
  }
}

void FUNC_REG_RATIO_analyze_step(Func* f, Bus* bus, int t) {

  // Local variables
  Branch* br;
  int* Hphi_nnz;
  Mat* Hphi;
  REAL da;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointer
  if (!Hphi_nnz || !Hphi)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      return;
    
    // Normalization factor
    da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
    if (da < FUNC_REG_RATIO_PARAM)
      da = FUNC_REG_RATIO_PARAM;
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var
      MAT_set_i(Hphi,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
      MAT_set_j(Hphi,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_REG_RATIO_eval_step(Func* f, Bus* bus, int t, Vec* var_values) {

  // Local variables
  Branch* br;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL a;
  REAL da;
  REAL a0;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      return;
    
    // Normalizatin factor
    da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
    if (da < FUNC_REG_RATIO_PARAM)
      da = FUNC_REG_RATIO_PARAM;
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var    
      
      a0 = BRANCH_get_ratio(br,t);
      a = VEC_get(var_values,BRANCH_get_index_ratio(br,t));
      (*phi) += 0.5*pow((a-a0)/da,2.);
      gphi[BRANCH_get_index_ratio(br,t)] = (a-a0)/(da*da);
      Hphi[*Hphi_nnz] = 1./(da*da);
      (*Hphi_nnz)++;
    }
    else {
      // nothing because a0-a0 = 0
    }
  }
}
