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
  FUNC_set_func_init(f,&FUNC_REG_RATIO_init);
  FUNC_set_func_count_step(f,&FUNC_REG_RATIO_count_step);
  FUNC_set_func_allocate(f,&FUNC_REG_RATIO_allocate);
  FUNC_set_func_clear(f,&FUNC_REG_RATIO_clear);
  FUNC_set_func_analyze_step(f,&FUNC_REG_RATIO_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_RATIO_eval_step);
  FUNC_set_func_free(f,&FUNC_REG_RATIO_free);
  FUNC_init(f);
  return f;
}

void FUNC_REG_RATIO_init(Func* f) {
  
  FUNC_set_name(f,"tap ratio regularization");
}

void FUNC_REG_RATIO_clear(Func* f) {
  
  // phi
  FUNC_set_phi(f,0);
  
  // gphi
  VEC_set_zero(FUNC_get_gphi(f));
  
  // Hphi
  // Constant so not clear it
  
  // Counter
  FUNC_set_Hphi_nnz(f,0);
}

void FUNC_REG_RATIO_count_step(Func* f, Branch* br, int t) {

  // Local variables
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) // ratio var
    (*Hphi_nnz)++;
}

void FUNC_REG_RATIO_allocate(Func* f) {
  
  // Local variables
  int num_vars;
  int Hphi_nnz;
  
  num_vars = NET_get_num_vars(FUNC_get_network(f));
  Hphi_nnz = FUNC_get_Hphi_nnz(f);

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));

  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  Hphi_nnz));
}

void FUNC_REG_RATIO_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  int* Hphi_nnz;
  Mat* H;
  REAL da;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointer
  if (!Hphi_nnz)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Normalization factor
  da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
  if (da < FUNC_REG_RATIO_PARAM)
    da = FUNC_REG_RATIO_PARAM;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var
    MAT_set_i(H,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
    MAT_set_j(H,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
    MAT_set_d(H,*Hphi_nnz,1./(da*da));
    (*Hphi_nnz)++;
  }
}

void FUNC_REG_RATIO_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  REAL* phi;
  REAL* gphi;
  REAL a;
  REAL da;
  REAL a0;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));

  // Check pointers
  if (!phi || !gphi)
    return;

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
  }
  else {
    // nothing because a0-a0 = 0
  }
}
    
void FUNC_REG_RATIO_free(Func* f) {
  // Nothing
}
