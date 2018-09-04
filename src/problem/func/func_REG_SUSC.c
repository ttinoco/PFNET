/** @file func_REG_SUSC.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_SUSC.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_SUSC.h>

Func* FUNC_REG_SUSC_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_REG_SUSC_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_REG_SUSC_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_SUSC_eval_step);
  FUNC_set_name(f,"susceptance regularization");
  return f;
}

void FUNC_REG_SUSC_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Shunt* shunt;
  int* Hphi_nnz;

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) // b var
      (*Hphi_nnz)++;
  }
}

void FUNC_REG_SUSC_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Shunt* shunt;
  int* Hphi_nnz;
  Mat* Hphi;
  REAL db;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

    db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt); // p.u.
    if (db < FUNC_REG_SUSC_PARAM)
      db = FUNC_REG_SUSC_PARAM;
    
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
      MAT_set_i(Hphi,*Hphi_nnz,SHUNT_get_index_b(shunt,t));
      MAT_set_j(Hphi,*Hphi_nnz,SHUNT_get_index_b(shunt,t));
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_REG_SUSC_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Shunt* shunt;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL b0;
  REAL b;
  REAL db;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    // Normalization factor
    db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt); // p.u.
    if (db < FUNC_REG_SUSC_PARAM)
      db = FUNC_REG_SUSC_PARAM;
    
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
      b0 = SHUNT_get_b(shunt,t);
      b = VEC_get(var_values,SHUNT_get_index_b(shunt,t));
      (*phi) += 0.5*pow((b-b0)/db,2.);
      gphi[SHUNT_get_index_b(shunt,t)] = (b-b0)/(db*db);
      Hphi[*Hphi_nnz] = 1./(db*db);
      (*Hphi_nnz)++;
    }
    else {
      // nothing because b0 - b0 = 0
    }
  }
}
