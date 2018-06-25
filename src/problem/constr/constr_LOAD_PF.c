/** @file constr_LOAD_PF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LOAD_PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LOAD_PF.h>

Constr* CONSTR_LOAD_PF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_LOAD_PF_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_LOAD_PF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_LOAD_PF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_LOAD_PF_store_sens_step);
  CONSTR_set_name(c,"load constant power factor");
  return c;
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
	  factor = sqrt((1.-gamma*gamma)/(gamma*gamma));

	  // A
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t));
	  if (LOAD_get_P(load,t)*LOAD_get_Q(load,t) >= 0)
	    MAT_set_d(A,*A_nnz,-factor);
	  else
	    MAT_set_d(A,*A_nnz,factor);
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
