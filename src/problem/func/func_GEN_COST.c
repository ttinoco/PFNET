/** @file func_GEN_COST.c
 *  @brief This file defines the data structure and routines associated with the function of type GEN_COST.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_GEN_COST.h>

Func* FUNC_GEN_COST_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_init(f, &FUNC_GEN_COST_init);
  FUNC_set_func_count_step(f, &FUNC_GEN_COST_count_step);
  FUNC_set_func_allocate(f, &FUNC_GEN_COST_allocate);
  FUNC_set_func_clear(f, &FUNC_GEN_COST_clear);
  FUNC_set_func_analyze_step(f, &FUNC_GEN_COST_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_GEN_COST_eval_step);
  FUNC_set_func_free(f, &FUNC_GEN_COST_free);
  FUNC_init(f);
  return f;
}

void FUNC_GEN_COST_init(Func* f) {
  
  FUNC_set_name(f,"generation cost");
}

void FUNC_GEN_COST_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);

  // gphi
  VEC_set_zero(FUNC_get_gphi(f));

  // Hphi
  // Constant so not clear it

  // Counter
  FUNC_set_Hphi_nnz(f,0);

  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_GEN_COST_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  (*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_GEN_COST_allocate(Func* f) {

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

void FUNC_GEN_COST_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  MAT_set_i(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  MAT_set_j(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  MAT_set_d(H,*Hphi_nnz,2.*GEN_get_cost_coeff_Q2(gen));
	  (*Hphi_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_GEN_COST_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_P;
  REAL P;
  REAL Q0;
  REAL Q1;
  REAL Q2;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

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
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	Q0 = GEN_get_cost_coeff_Q0(gen);
	Q1 = GEN_get_cost_coeff_Q1(gen);
	Q2 = GEN_get_cost_coeff_Q2(gen);

	// Variable
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {

	  // Index
	  index_P = GEN_get_index_P(gen,t);

	  // P
	  P = VEC_get(var_values,index_P);

	  // phi
	  (*phi) += Q0 + Q1*P + Q2*pow(P,2.);

	  // gphi
	  gphi[index_P] = Q1 + 2.*Q2*P;
	}

	// Constant
	else {

	  // P
	  P = GEN_get_P(gen,t);

	  // phi
	  (*phi) += Q0 + Q1*P + Q2*pow(P,2.);
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_GEN_COST_free(Func* f) {
  // Nothing
}
