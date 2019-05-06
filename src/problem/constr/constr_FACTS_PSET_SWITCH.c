/** @file constr_FACTS_PSET_SWITCH.c
 *  @brief This file defines the data structure and routines associated with the constraint of type FACTS_PSET_SWITCH.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_FACTS_PSET_SWITCH.h>

Constr* CONSTR_FACTS_PSET_SWITCH_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_FACTS_PSET_SWITCH_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_FACTS_PSET_SWITCH_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_FACTS_PSET_SWITCH_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_FACTS_PSET_SWITCH_store_sens_step);
  CONSTR_set_name(c,"switching FACTS active power control");
  return c;
}

void CONSTR_FACTS_PSET_SWITCH_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Facts* facts;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts !=NULL; facts = FACTS_get_next_k(facts)) {

    if (FACTS_is_in_normal_series_mode(facts) &&
        FACTS_get_P_max_dc(facts) > 0. &&
        FACTS_is_in_service(facts) &&
        FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P)) {
      (*A_nnz) += 1;
      (*A_row) += 1;
    }
  }
}

void CONSTR_FACTS_PSET_SWITCH_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* A_nnz;
  int* A_row;
  Facts* facts;
  Vec* b;
  Mat* A;

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts !=NULL; facts = FACTS_get_next_k(facts)) {

    if (FACTS_is_in_normal_series_mode(facts) &&
        FACTS_get_P_max_dc(facts) > 0. &&
        FACTS_is_in_service(facts) &&
        FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P)) {
            
      VEC_set(b,*A_row,FACTS_get_P_set(facts,t));
          
      // P_m
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,FACTS_get_index_P_m(facts,t));
      MAT_set_d(A,*A_nnz,1.);	    
      (*A_nnz)++;
          
      (*A_row)++;
    }
  }
}

void CONSTR_FACTS_PSET_SWITCH_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_FACTS_PSET_SWITCH_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
