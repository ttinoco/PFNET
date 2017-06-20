/** @file constr_LOAD_PF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LOAD_PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LOAD_PF.h>

Constr* CONSTR_LOAD_PF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_LOAD_PF_init);
  CONSTR_set_func_count_step(c,&CONSTR_LOAD_PF_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_LOAD_PF_allocate);
  CONSTR_set_func_clear(c,&CONSTR_LOAD_PF_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_LOAD_PF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_LOAD_PF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_LOAD_PF_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_LOAD_PF_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_LOAD_PF_init(Constr* c) {

  // Init
  CONSTR_set_name(c,"load constant power factor");
  CONSTR_set_data(c,NULL);
}

void CONSTR_LOAD_PF_clear(Constr* c) {

  // Counters
  CONSTR_set_A_nnz(c,0);
  CONSTR_set_A_row(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_LOAD_PF_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
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

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Variables
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P) && LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {
	  (*A_nnz)++; // P
	  (*A_nnz)++; // Q
	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_LOAD_PF_allocate(Constr* c) {

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

  // A b
  CONSTR_set_A(c,MAT_new(num_constr, // rows
			 num_vars,   // cols
			 A_nnz));    // nnz
  CONSTR_set_b(c,VEC_new(num_constr));

  // G l u
  CONSTR_set_l(c,VEC_new(0));
  CONSTR_set_u(c,VEC_new(0));
  CONSTR_set_G(c,MAT_new(0,        // rows
			 num_vars, // cols
			 0));      // nnz
}

void CONSTR_LOAD_PF_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Load* load;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  REAL gamma;
  REAL factor;
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

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Variables
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P) && LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {

	  gamma = LOAD_get_target_power_factor(load);
	  factor = sqrt(1.-gamma*gamma)/gamma;

	  // A
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t));
	  MAT_set_d(A,*A_nnz,-factor);
	  (*A_nnz)++;

	  // A
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_Q(load,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;

	  // b
	  VEC_set(b,*A_row,0.);

	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_LOAD_PF_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_LOAD_PF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing for now
}

void CONSTR_LOAD_PF_free(Constr* c) {
  // Nothing to do
}
