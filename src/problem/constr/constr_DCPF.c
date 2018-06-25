/** @file constr_DCPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_DCPF.h>

Constr* CONSTR_DCPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_DCPF_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_DCPF_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_DCPF_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_DCPF_store_sens_step);
  CONSTR_set_name(c,"DC power balance");
  return c;
}

void CONSTR_DCPF_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Load* load;
  Vargen* vargen;
  Bat* bat;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int m;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(bus[k],t);

  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    // Outage
    if (BRANCH_is_on_outage(br))
      break;

    if (k == 0)
      m = 1;
    else
      m = 0;

    //***********
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // wk var

      // A
      (*A_nnz)++; // Pk
    }

    //***********
    if (BUS_has_flags(bus[m],FLAG_VARS,BUS_VAR_VANG)) { // wm var

      // A
      (*A_nnz)++; // Pk
    }

    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var

      // A
      (*A_nnz)++; // Pk
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var

	  // A
	  (*A_nnz)++; // Pk
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // P var

	  // A
	  (*A_nnz)++; // Pk
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // P var

	  // A
	  (*A_nnz)++; // Pk
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var

	  // A
	  (*A_nnz)++; // Pc
	  (*A_nnz)++; // Pd
	}
      }

      // Counter constraints
      (*A_row)++;
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_DCPF_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Mat* A;
  Vec* rhs;
  int* A_nnz;
  char* bus_counted;
  int bus_index_t[2];
  REAL b;
  REAL sign_phi;
  int k;
  int m;

  // Constr data
  A = CONSTR_get_A(c);
  rhs = CONSTR_get_b(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(bus[k],t);

  // Branch data
  b = BRANCH_get_b(br);

  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    // Outage
    if (BRANCH_is_on_outage(br))
      break;

    if (k == 0) {
      m = 1;
      sign_phi = 1;
    }
    else {
      m = 0;
      sign_phi = -1;
    }

    //***********
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // wk var

      // A
      MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus[k],t)); // wk
      MAT_set_d(A,*A_nnz,b);
      (*A_nnz)++;
    }
    else {

      // b
      VEC_add_to_entry(rhs,bus_index_t[k],-b*BUS_get_v_ang(bus[k],t));
    }

    //***********
    if (BUS_has_flags(bus[m],FLAG_VARS,BUS_VAR_VANG)) { // wm var

      // A
      MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus[m],t)); // wk
      MAT_set_d(A,*A_nnz,-b);
      (*A_nnz)++;
    }
    else {

      // b
      VEC_add_to_entry(rhs,bus_index_t[k],b*BUS_get_v_ang(bus[m],t));
    }

    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var

      // A
      MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
      MAT_set_j(A,*A_nnz,BRANCH_get_index_phase(br,t)); // phi
      MAT_set_d(A,*A_nnz,-b*sign_phi);
      (*A_nnz)++;
    }
    else {

      // b
      VEC_add_to_entry(rhs,bus_index_t[k],b*BRANCH_get_phase(br,t)*sign_phi);
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var

	  // A
	  MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
	  MAT_set_j(A,*A_nnz,GEN_get_index_P(gen,t)); // Pg
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	}
	else {

	  // b
	  VEC_add_to_entry(rhs,bus_index_t[k],-GEN_get_P(gen,t));
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var

	  // A
	  MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
	  MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t)); // Pl
	  MAT_set_d(A,*A_nnz,-1.);
	  (*A_nnz)++;
	}
	else {

	  // b
	  VEC_add_to_entry(rhs,bus_index_t[k],LOAD_get_P(load,t));
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // A
	  MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
	  MAT_set_j(A,*A_nnz,VARGEN_get_index_P(vargen,t)); // Pg
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	}
	else {

	  // b
	  VEC_add_to_entry(rhs,bus_index_t[k],-VARGEN_get_P(vargen,t));
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var

	  // A
	  MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pc(bat,t)); // Pc
	  MAT_set_d(A,*A_nnz,-1.);
	  (*A_nnz)++;

	  // A
	  MAT_set_i(A,*A_nnz,bus_index_t[k]); // Pk
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pd(bat,t)); // Pd
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	}
	else {

	  // b
	  VEC_add_to_entry(rhs,bus_index_t[k],BAT_get_P(bat,t));
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_DCPF_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_DCPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* bus[2];
  int bus_index_t[2];
  char* bus_counted;
  int k;

  // Constr data
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index_t(bus[k],t);

  // Buses
  for (k = 0; k < 2; k++) {

    // Store P balance sensitivity
    if (!bus_counted[bus_index_t[k]])
      BUS_set_sens_P_balance(bus[k],VEC_get(sA,bus_index_t[k]),t);

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}
