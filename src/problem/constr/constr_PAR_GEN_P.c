/** @file constr_PAR_GEN_P.c
 *  @brief This file defines the data structure and routines associated with the constraint of type PAR_GEN_P.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_PAR_GEN_P.h>

Constr* CONSTR_PAR_GEN_P_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_PAR_GEN_P_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_PAR_GEN_P_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_PAR_GEN_P_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_PAR_GEN_P_store_sens_step);
  CONSTR_set_name(c,"generator active power participation");
  return c;
}

void CONSTR_PAR_GEN_P_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen1;
  Gen* gen2;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Active power of slack generators
      if (BUS_is_slack(bus)) {
	for (gen1 = BUS_get_gen(bus); gen1 != NULL; gen1 = GEN_get_next(gen1)) {
	  if (!GEN_is_on_outage(gen1) && GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P))
	    break;
	}
	for (gen2 = GEN_get_next(gen1); gen2 != NULL; gen2 = GEN_get_next(gen2)) {
	  if (GEN_is_on_outage(gen2))
	    continue;
	  if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_P)) {
	    (*A_nnz)++;
	    (*A_nnz)++;
	    (*A_row)++;
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_P_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen1;
  Gen* gen2;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  Vec* b;
  Mat* A;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Active power of slack generators
      if (BUS_is_slack(bus)) {
	for (gen1 = BUS_get_gen(bus); gen1 != NULL; gen1 = GEN_get_next(gen1)) {
	  if (!GEN_is_on_outage(gen1) && GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P))
	    break;
	}
	for (gen2 = GEN_get_next(gen1); gen2 != NULL; gen2 = GEN_get_next(gen2)) {
	  if (GEN_is_on_outage(gen2))
	    continue;
	  if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_P)) {
	    VEC_set(b,*A_row,0.);
	    MAT_set_i(A,*A_nnz,*A_row);
	    MAT_set_j(A,*A_nnz,GEN_get_index_P(gen1,t));
	    MAT_set_d(A,*A_nnz,1.);
	    (*A_nnz)++;
	    MAT_set_i(A,*A_nnz,*A_row);
	    MAT_set_j(A,*A_nnz,GEN_get_index_P(gen2,t));
	    MAT_set_d(A,*A_nnz,-1.);
	    (*A_nnz)++;
	    (*A_row)++;
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_P_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_PAR_GEN_P_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
