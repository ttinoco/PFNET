/** @file func_REG_PHASE.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_PHASE.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_PHASE.h>

void FUNC_REG_PHASE_init(Func* f) {
  // Nothing
}

void FUNC_REG_PHASE_clear(Func* f) {
  
  // phi
  FUNC_set_phi(f,0);
  
  // gphi
  VEC_set_zero(FUNC_get_gphi(f));
  
  // Hphi
  // Constant so not clear it
  
  // Counter
  FUNC_set_Hcounter(f,0);
}

void FUNC_REG_PHASE_count_branch(Func* f, Branch *br) {

  // Local variables
  int* Hcounter;

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);
  if (!Hcounter)
    return;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) // phase var
    (*Hcounter)++;
}

void FUNC_REG_PHASE_allocate(Func* f) {
  
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

void FUNC_REG_PHASE_analyze_branch(Func* f, Branch *br) {

  // Local variables
  int* Hcounter;
  Mat* H;
  REAL dp;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);
  if (!Hcounter)
    return;

  // Normalization factor
  dp = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
  if (dp < FUNC_REG_PHASE_PARAM)
    dp = FUNC_REG_PHASE_PARAM;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phase var
    MAT_set_i(H,*Hcounter,BRANCH_get_index_phase(br));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_phase(br));
    MAT_set_d(H,*Hcounter,1./(dp*dp));
    (*Hcounter)++;
  }
}

void FUNC_REG_PHASE_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  REAL* phi;
  REAL* gphi;
  REAL p;
  REAL dp;
  REAL p0;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  if (!phi || !gphi)
    return;
  
  // Normalizatin factor
  dp = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
  if (dp < FUNC_REG_PHASE_PARAM)
    dp = FUNC_REG_PHASE_PARAM;

  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phase var    
    
    p0 = BRANCH_get_phase(br);
    p = VEC_get(var_values,BRANCH_get_index_phase(br));
    (*phi) += 0.5*pow((p-p0)/dp,2.);
    gphi[BRANCH_get_index_phase(br)] = (p-p0)/(dp*dp);
  }
}
    
void FUNC_REG_PHASE_free(Func* f) {
  // Nothing
}
