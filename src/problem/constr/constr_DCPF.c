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

void CONSTR_DCPF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Branch* br;
  Bus* bus_k;
  Bus* bus_m;
  Gen* gen;
  Load* load;
  Vargen* vargen;
  Bat* bat;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // dP index
  BUS_set_dP_index(bus,*A_row,t);

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Out of service
    if (!GEN_is_in_service(gen))
      continue;
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
      
      // A
      (*A_nnz)++; // Pk
    }
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // P var
      
      // A
      (*A_nnz)++; // Pk
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

    // Out of service
    if (!VARGEN_is_in_service(vargen))
      continue;
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // P var
      
      // A
      (*A_nnz)++; // Pk
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

    // Out of service
    if (!BAT_is_in_service(bat))
      continue;
    
    //*****************************
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var
      
      // A
      (*A_nnz)++; // Pc
      (*A_nnz)++; // Pd
    }
  }

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    // Bus data
    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);
  
    //***********
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) { // wk var
	
      // A
      (*A_nnz)++; // Pk
    }
    
    //***********
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A
      (*A_nnz)++; // Pk
    }
    
    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // A
      (*A_nnz)++; // Pk
    }
  }
  for (br = BUS_get_branch_m(bus); br != NULL; br = BRANCH_get_next_m(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    // Bus data
    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);
  
    //***********
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A
      (*A_nnz)++; // Pm
    }
    
    //***********
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A
      (*A_nnz)++; // Pm
    }
    
    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // A
      (*A_nnz)++; // Pm
    }
  }
  
  // Rows
  (*A_row)++;
}

void CONSTR_DCPF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Branch* br;
  Bus* bus_k;
  Bus* bus_m;
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Mat* A;
  Vec* rhs;
  int* A_nnz;
  int* A_row;
  REAL b;

  // Constr data
  A = CONSTR_get_A(c);
  rhs = CONSTR_get_b(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Out of service
    if (!GEN_is_in_service(gen))
      continue;
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pk
      MAT_set_j(A,*A_nnz,GEN_get_index_P(gen,t)); // Pg
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,-GEN_get_P(gen,t));
    }
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pk
      MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t)); // Pl
      MAT_set_d(A,*A_nnz,-1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,LOAD_get_P(load,t));
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

    // Out of service
    if (!VARGEN_is_in_service(vargen))
      continue;
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pk
      MAT_set_j(A,*A_nnz,VARGEN_get_index_P(vargen,t)); // Pg
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,-VARGEN_get_P(vargen,t));
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

    // Out of service
    if (!BAT_is_in_service(bat))
      continue;
    
    //*****************************
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pk
      MAT_set_j(A,*A_nnz,BAT_get_index_Pc(bat,t)); // Pc
      MAT_set_d(A,*A_nnz,-1.);
      (*A_nnz)++;
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pk
      MAT_set_j(A,*A_nnz,BAT_get_index_Pd(bat,t)); // Pd
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,BAT_get_P(bat,t));
    }
  }
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;

    // Bus data
    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);
    
    // Branch data
    b = BRANCH_get_b(br);
          
    //***********
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pkm
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus_k,t)); // wk
      MAT_set_d(A,*A_nnz,b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,-b*BUS_get_v_ang(bus_k,t));
    }
    
    //***********
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pkm
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus_m,t)); // wm
      MAT_set_d(A,*A_nnz,-b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,b*BUS_get_v_ang(bus_m,t));
    }
    
    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pkm
      MAT_set_j(A,*A_nnz,BRANCH_get_index_phase(br,t)); // phi
      MAT_set_d(A,*A_nnz,-b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,b*BRANCH_get_phase(br,t));
    }
  }
  for (br = BUS_get_branch_m(bus); br != NULL; br = BRANCH_get_next_m(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;

    // Bus data
    bus_k = BRANCH_get_bus_k(br);
    bus_m = BRANCH_get_bus_m(br);
    
    // Branch data
    b = BRANCH_get_b(br);
          
    //***********
    if (BUS_has_flags(bus_m,FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pmk
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus_m,t)); // wm
      MAT_set_d(A,*A_nnz,b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,-b*BUS_get_v_ang(bus_m,t));
    }
    
    //***********
    if (BUS_has_flags(bus_k,FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pmk
      MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus_k,t)); // wk
      MAT_set_d(A,*A_nnz,-b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,b*BUS_get_v_ang(bus_k,t));
    }
    
    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // A
      MAT_set_i(A,*A_nnz,*A_row); // Pmk
      MAT_set_j(A,*A_nnz,BRANCH_get_index_phase(br,t)); // phi
      MAT_set_d(A,*A_nnz,b);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(rhs,*A_row,-b*BRANCH_get_phase(br,t));
    }
  }
  
  // Rows
  (*A_row)++;
}

void CONSTR_DCPF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_DCPF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Locals
  int* A_row;

  // Data
  A_row = CONSTR_get_A_row_ptr(c);

  // Check
  if (!A_row || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;
  
  BUS_set_sens_P_balance(bus,VEC_get(sA,*A_row),t);

  // Rows
  (*A_row)++;
}
