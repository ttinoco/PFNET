/** @file func_LOAD_UTIL.c
 *  @brief This file defines the data structure and routines associated with the function of type LOAD_UTIL.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_LOAD_UTIL.h>

void FUNC_LOAD_UTIL_init(Func* f) {
  // Nothing
}

void FUNC_LOAD_UTIL_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);

  // gphi
  VEC_set_zero(FUNC_get_gphi(f));

  // Hphi
  // Constant so not clear it

  // Counter
  FUNC_set_Hcounter(f,0);

  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_LOAD_UTIL_count_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  int k;

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  (*Hcounter)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_LOAD_UTIL_allocate(Func* f) {

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

void FUNC_LOAD_UTIL_analyze_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;

  // Constr data
  H = FUNC_get_Hphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  MAT_set_i(H,*Hcounter,LOAD_get_index_P(load));
	  MAT_set_j(H,*Hcounter,LOAD_get_index_P(load));
	  MAT_set_d(H,*Hcounter,2.*LOAD_get_util_coeff_Q2(load));
	  (*Hcounter)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_LOAD_UTIL_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int bus_index[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_P;
  REAL P;
  REAL Q0;
  REAL Q1;
  REAL Q2;
  int k;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !gphi || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	Q0 = LOAD_get_util_coeff_Q0(load);
	Q1 = LOAD_get_util_coeff_Q1(load);
	Q2 = LOAD_get_util_coeff_Q2(load);

	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {

	  // Index
	  index_P = LOAD_get_index_P(load);

	  // P
	  P = VEC_get(var_values,index_P);

	  // phi
	  (*phi) += Q0 + Q1*P + Q2*pow(P,2.);

	  // gphi
	  gphi[index_P] = Q1 + 2.*Q2*P;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_LOAD_UTIL_free(Func* f) {
  // Nothing
}
