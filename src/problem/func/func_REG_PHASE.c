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

Func* FUNC_REG_PHASE_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f, &FUNC_REG_PHASE_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_REG_PHASE_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_REG_PHASE_eval_step);
  FUNC_set_name(f,"phase shift regularization");
  return f;
}

void FUNC_REG_PHASE_count_step(Func* f, Bus* bus, int t) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointer
  if (!Hphi_nnz)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      return;
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) // phase var
      (*Hphi_nnz)++;
  }
}

void FUNC_REG_PHASE_analyze_step(Func* f, Bus* bus, int t) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;
  Mat* Hphi;
  REAL dp;

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
    dp = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
    if (dp < FUNC_REG_PHASE_PARAM)
      dp = FUNC_REG_PHASE_PARAM;
  
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phase var
      MAT_set_i(Hphi,*Hphi_nnz,BRANCH_get_index_phase(br,t));
      MAT_set_j(Hphi,*Hphi_nnz,BRANCH_get_index_phase(br,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_REG_PHASE_eval_step(Func* f, Bus* bus, int t, Vec* var_values) {

  // Local variables
  int* Hphi_nnz;
  Branch* br;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  REAL p;
  REAL dp;
  REAL p0;

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
    dp = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
    if (dp < FUNC_REG_PHASE_PARAM)
      dp = FUNC_REG_PHASE_PARAM;
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phase var    
      
      p0 = BRANCH_get_phase(br,t);
      p = VEC_get(var_values,BRANCH_get_index_phase(br,t));
      (*phi) += 0.5*pow((p-p0)/dp,2.);
      gphi[BRANCH_get_index_phase(br,t)] = (p-p0)/(dp*dp);
      Hphi[*Hphi_nnz] = 1./(dp*dp);
      (*Hphi_nnz)++;
    }
    else {
      // nothing because p0-p0 = 0
    }
  }
}
