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
  FUNC_set_func_count_step(f, &FUNC_REG_VANG_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_REG_VANG_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_REG_VANG_eval_step);
  FUNC_set_name(f,"voltage angle regularization");
  return f;
}

void FUNC_REG_VANG_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  int bus_index_t[2];
  Bus* bus;
  int* Hphi_nnz;
  char* bus_counted;
  int k;

  // Func data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(buses[k],t);

  // Branch
  if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br)) { // wk var
    (*Hphi_nnz)++; // wk and wk

    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br))
      (*Hphi_nnz)++; // wk and wm
  }
  if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br)) // wm var
    (*Hphi_nnz)++; // wm and wm

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];
    
    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) // w var
	(*Hphi_nnz)++; // w and w
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VANG_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  int bus_index_t[2];
  int index_v_ang[2];
  Bus* bus;
  int* Hphi_nnz;
  char* bus_counted;
  Mat* Hphi;
  int k;

  // Func data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted || !Hphi)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index_t(buses[k],t);
    index_v_ang[k] = BUS_get_index_v_ang(buses[k],t);
  }

  // Branch
  if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br)) { // wk var
    MAT_set_i(Hphi,*Hphi_nnz,index_v_ang[0]);
    MAT_set_j(Hphi,*Hphi_nnz,index_v_ang[0]);
    (*Hphi_nnz)++; // wk and wk

    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br)) {
      MAT_set_i(Hphi,*Hphi_nnz,index_v_ang[0]);
      MAT_set_j(Hphi,*Hphi_nnz,index_v_ang[1]);
      (*Hphi_nnz)++; // wk and wm
    }
  }
  if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG) && !BRANCH_is_on_outage(br)) { // wm var
    MAT_set_i(Hphi,*Hphi_nnz,index_v_ang[1]);
    MAT_set_j(Hphi,*Hphi_nnz,index_v_ang[1]);
    (*Hphi_nnz)++; // wm and wm
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) { // v var
	MAT_set_i(Hphi,*Hphi_nnz,index_v_ang[k]);
	MAT_set_j(Hphi,*Hphi_nnz,index_v_ang[k]);
	(*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VANG_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  int bus_index_t[2];
  int index_v_ang[2];
  REAL w[2];
  BOOL var_w[2];
  Bus* bus;
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  REAL shift;
  int k;
  REAL wdiff;
  REAL dw = FUNC_REG_VANG_PARAM;

  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !gphi || !bus_counted || !Hphi || !Hphi_nnz)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index_t(buses[k],t);
    index_v_ang[k] = BUS_get_index_v_ang(buses[k],t);
    var_w[k] = BUS_has_flags(buses[k],FLAG_VARS,BUS_VAR_VANG);
    if (var_w[k])
      w[k] = VEC_get(var_values,index_v_ang[k]);
    else
      w[k] = BUS_get_v_ang(buses[k],t);
  }

  // No branch outage
  if (!BRANCH_is_on_outage(br)) {

    // Branch data
    shift = BRANCH_get_phase(br,t); // radians

    // Difference
    wdiff = w[0]-w[1]-shift;

    // phi
    (*phi) += 0.5*pow(wdiff/dw,2.);

    // gphi
    if (var_w[0]) { // wk var
      gphi[index_v_ang[0]] += wdiff/(dw*dw);
      Hphi[*Hphi_nnz] = 1./(dw*dw);
      (*Hphi_nnz)++;
      if (var_w[1]) {
	Hphi[*Hphi_nnz] = -1./(dw*dw);
	(*Hphi_nnz)++;
      }
    }
    if (var_w[1]) { // wm var
      gphi[index_v_ang[1]] -= wdiff/(dw*dw);
      Hphi[*Hphi_nnz] = 1./(dw*dw);
      (*Hphi_nnz)++;
    }
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // phi
      (*phi) += 0.5*pow(w[k]/dw,2.);

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) { // v var

	// gphi
	gphi[index_v_ang[k]] += w[k]/(dw*dw);

	// Hphi
	Hphi[*Hphi_nnz] = 1./(dw*dw);
	(*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}
