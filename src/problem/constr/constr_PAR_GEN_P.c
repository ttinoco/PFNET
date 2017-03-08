/** @file constr_PAR_GEN_P.c
 *  @brief This file defines the data structure and routines associated with the constraint of type PAR_GEN_P.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_PAR_GEN_P.h>
#include <assert.h>

void CONSTR_PAR_GEN_P_init(Constr* c) {

  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_PAR_GEN_P_clear(Constr* c) {

  // Counters
  CONSTR_set_A_nnz(c,0);
  CONSTR_set_A_row(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
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

  // Check outage
  if (BRANCH_is_on_outage(br))
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
	gen1 = BUS_get_gen(bus);
	for (gen2 = GEN_get_next(gen1); gen2 != NULL; gen2 = GEN_get_next(gen2)) {
	  if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P))
	    (*A_nnz)++;
	  if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_P))
	    (*A_nnz)++;
	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_P_allocate(Constr* c) {

  // Local variables
  int num_constr;
  int num_vars;
  int A_nnz;

  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_constr = CONSTR_get_A_row(c);
  A_nnz = CONSTR_get_A_nnz(c);

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // b
  CONSTR_set_b(c,VEC_new(num_constr));

  // A
  CONSTR_set_A(c,MAT_new(num_constr, // size1 (rows)
			 num_vars,   // size2 (rows)
			 A_nnz)); // nnz
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

  // Check outage
  if (BRANCH_is_on_outage(br))
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
	gen1 = BUS_get_gen(bus);
	for (gen2 = GEN_get_next(gen1); gen2 != NULL; gen2 = GEN_get_next(gen2)) {
	  VEC_set(b,*A_row,0.);
	  if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P)) {
	    MAT_set_i(A,*A_nnz,*A_row);
	    MAT_set_j(A,*A_nnz,GEN_get_index_P(gen1,t));
	    MAT_set_d(A,*A_nnz,1.);
	    (*A_nnz)++;
	  }
	  else
	    VEC_add_to_entry(b,*A_row,-GEN_get_P(gen1,t));
	  if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_P)) {
	    MAT_set_i(A,*A_nnz,*A_row);
	    MAT_set_j(A,*A_nnz,GEN_get_index_P(gen2,t));
	    MAT_set_d(A,*A_nnz,-1.);
	    (*A_nnz)++;
	  }
	  else
	    VEC_add_to_entry(b,*A_row,GEN_get_P(gen2,t));
	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_P_eval_step(Constr* c, Branch* br, int t, Vec* values) {
  // Nothing to do
}

void CONSTR_PAR_GEN_P_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_PAR_GEN_P_free(Constr* c) {
  // Nothing to do
}
