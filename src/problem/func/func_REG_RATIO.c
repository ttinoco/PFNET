/** @file func_REG_RATIO.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_RATIO.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_RATIO.h>

Func* FUNC_REG_RATIO_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_name(f,"tap ratio regularization");
  FUNC_set_func_init(f, &FUNC_REG_RATIO_init);
  FUNC_set_func_count_step(f, &FUNC_REG_RATIO_count_step);
  FUNC_set_func_allocate(f, &FUNC_REG_RATIO_allocate);
  FUNC_set_func_clear(f, &FUNC_REG_RATIO_clear);
  FUNC_set_func_analyze_step(f, &FUNC_REG_RATIO_analyze_step);
  FUNC_set_func_eval_setp(f, &FUNC_REG_RATIO_eval_step);
  FUNC_set_func_free(f, &FUNC_REG_RATIO_free);
  return f;
}

void FUNC_REG_RATIO_init(Func* f) {
  // Nothing
}

void FUNC_REG_RATIO_clear(Func* f) {
  
  // phi
  FUNC_set_phi(f,0);
  
  // gphi
  VEC_set_zero(FUNC_get_gphi(f));
  
  // Hphi
  // Constant so not clear it
  
  // Counter
  FUNC_set_Hcounter(f,0);
}

void FUNC_REG_RATIO_count_step(Func* f, Branch* br, int t) {

  // Local variables
  int* Hcounter;

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);

  // Check pointers
  if (!Hcounter)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) // ratio var
    (*Hcounter)++;

  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // ratio dev var
    (*Hcounter)++;
    (*Hcounter)++;
  }
}

void FUNC_REG_RATIO_allocate(Func* f) {
  
  // Local variables
  int num_vars;
  int Hcounter;
  
  num_vars = NET_get_num_vars(FUNC_get_network(f));
  Hcounter = FUNC_get_Hcounter(f);

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));

  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  Hcounter));
}

void FUNC_REG_RATIO_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  int* Hcounter;
  Mat* H;
  REAL da;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);

  // Check pointer
  if (!Hcounter)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Normalization factor
  da = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
  if (da < FUNC_REG_RATIO_PARAM)
    da = FUNC_REG_RATIO_PARAM;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var
    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio(br,t));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio(br,t));
    MAT_set_d(H,*Hcounter,1./(da*da));
    (*Hcounter)++;
  }
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // ratio dev var
   
    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio_y(br,t));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio_y(br,t));
    MAT_set_d(H,*Hcounter,1./(da*da));
    (*Hcounter)++;

    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio_z(br,t));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio_z(br,t));
    MAT_set_d(H,*Hcounter,1./(da*da));
    (*Hcounter)++;
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
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // ratio dev var
    
    a = VEC_get(var_values,BRANCH_get_index_ratio_y(br,t));
    (*phi) += 0.5*pow(a/da,2.);
    gphi[BRANCH_get_index_ratio_y(br,t)] = a/(da*da);

    a = VEC_get(var_values,BRANCH_get_index_ratio_z(br,t));
    (*phi) += 0.5*pow(a/da,2.);
    gphi[BRANCH_get_index_ratio_z(br,t)] = a/(da*da);
  }
  else {
    // nothing becuase a0-a0 = 0
  }
}
    
void FUNC_REG_RATIO_free(Func* f) {
  // Nothing
}
