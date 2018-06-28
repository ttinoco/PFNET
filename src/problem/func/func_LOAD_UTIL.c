/** @file func_LOAD_UTIL.c
 *  @brief This file defines the data structure and routines associated with the function of type LOAD_UTIL.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_LOAD_UTIL.h>

Func* FUNC_LOAD_UTIL_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f, &FUNC_LOAD_UTIL_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_LOAD_UTIL_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_LOAD_UTIL_eval_step);
  FUNC_set_name(f,"consumption utility");
  return f;
}

void FUNC_LOAD_UTIL_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  int k;

  // Constr data
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
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  (*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_LOAD_UTIL_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;

  // Constr data
  H = FUNC_get_Hphi(f);
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
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  MAT_set_i(H,*Hphi_nnz,LOAD_get_index_P(load,t));
	  MAT_set_j(H,*Hphi_nnz,LOAD_get_index_P(load,t));
	  (*Hphi_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_LOAD_UTIL_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_P;
  REAL P;
  REAL Q0;
  REAL Q1;
  REAL Q2;
  int k;

  // Constr data
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

      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	Q0 = LOAD_get_util_coeff_Q0(load);
	Q1 = LOAD_get_util_coeff_Q1(load);
	Q2 = LOAD_get_util_coeff_Q2(load);

	// Variable
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {

	  // Index
	  index_P = LOAD_get_index_P(load,t);

	  // P
	  P = VEC_get(var_values,index_P);

	  // phi
	  (*phi) += Q0 + Q1*P + Q2*pow(P,2.);

	  // gphi
	  gphi[index_P] = Q1 + 2.*Q2*P;

	  // Hphi
	  Hphi[*Hphi_nnz] = 2.*Q2;
	  (*Hphi_nnz)++;
	}

	// Constant
	else {

	  // P
	  P = LOAD_get_P(load,t);

	  // phi
	  (*phi) += Q0 + Q1*P + Q2*pow(P,2.);
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}
