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

Func* FUNC_REG_VMAG_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_REG_VMAG_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_REG_VMAG_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_VMAG_eval_step);
  FUNC_set_name(f,"voltage magnitude regularization");
  return f;
}

void FUNC_REG_VMAG_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
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

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) // v var
	(*Hphi_nnz)++;
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
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
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(buses[k],t);

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	MAT_set_i(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	(*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_v_mag;
  REAL v;
  REAL vt;
  REAL dv = FUNC_REG_VMAG_PARAM;
  int k;

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
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(buses[k],t);

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // Set point
      vt = BUS_get_v_set(bus,t);

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	// Index
	index_v_mag = BUS_get_index_v_mag(bus,t);

	// v
	v = VEC_get(var_values,index_v_mag);

	// phi
	(*phi) += 0.5*pow((v-vt)/dv,2.);

	// gphi
	gphi[index_v_mag] = (v-vt)/(dv*dv);

	// Hphi
	Hphi[*Hphi_nnz] = 1./(dv*dv);
	(*Hphi_nnz)++;
      }
      else {

	// v
	v = BUS_get_v_mag(bus,t);

	// phi
	(*phi) += 0.5*pow((v-vt)/dv,2.);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}
