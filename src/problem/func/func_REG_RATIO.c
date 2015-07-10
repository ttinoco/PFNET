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

void FUNC_REG_RATIO_count_branch(Func* f, Branch *br) {

  // Local variables
  int* Hcounter;

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);
  if (!Hcounter)
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

void FUNC_REG_RATIO_analyze_branch(Func* f, Branch *br) {

  // Local variables
  int* Hcounter;
  Mat* H;
  REAL dt;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);
  if (!Hcounter)
    return;

  // Normalization factor
  dt = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
  if (dt < FUNC_REG_RATIO_PARAM)
    dt = FUNC_REG_RATIO_PARAM;
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var
    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio(br));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio(br));
    MAT_set_d(H,*Hcounter,1./(dt*dt));
    (*Hcounter)++;
  }
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // ratio dev var
    
    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio_y(br));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio_y(br));
    MAT_set_d(H,*Hcounter,1./(dt*dt));
    (*Hcounter)++;

    MAT_set_i(H,*Hcounter,BRANCH_get_index_ratio_z(br));
    MAT_set_j(H,*Hcounter,BRANCH_get_index_ratio_z(br));
    MAT_set_d(H,*Hcounter,1./(dt*dt));
    (*Hcounter)++;
  }
}

void FUNC_REG_RATIO_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  REAL* phi;
  REAL* gphi;
  REAL t;
  REAL dt;
  REAL t0;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  if (!phi || !gphi)
    return;
  
  // Normalizatin factor
  dt = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
  if (dt < FUNC_REG_RATIO_PARAM)
    dt = FUNC_REG_RATIO_PARAM;

  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // ratio var    
    
    t0 = BRANCH_get_ratio(br);
    t = VEC_get(var_values,BRANCH_get_index_ratio(br));
    (*phi) += 0.5*pow((t-t0)/dt,2.);
    gphi[BRANCH_get_index_ratio(br)] = (t-t0)/(dt*dt);
  }

  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // ratio dev var
    
    t = VEC_get(var_values,BRANCH_get_index_ratio_y(br));
    (*phi) += 0.5*pow(t/dt,2.);
    gphi[BRANCH_get_index_ratio_y(br)] = t/(dt*dt);

    t = VEC_get(var_values,BRANCH_get_index_ratio_z(br));
    (*phi) += 0.5*pow(t/dt,2.);
    gphi[BRANCH_get_index_ratio_z(br)] = t/(dt*dt);
  }
}
    
void FUNC_REG_RATIO_free(Func* f) {
  // Nothing
}
