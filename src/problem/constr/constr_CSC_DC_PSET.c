/** @file constr_CSC_DC_PSET.c
 *  @brief This file defines the data structure and routines associated with the constraint of type dc_pset.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_CSC_DC_PSET.h>

Constr* CONSTR_CSC_DC_PSET_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_CSC_DC_PSET_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_CSC_DC_PSET_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_CSC_DC_PSET_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_CSC_DC_PSET_store_sens_step);
  CONSTR_set_name(c,"CSC DC power control");
  return c;
}

void CONSTR_CSC_DC_PSET_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvCSC* conv;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;

  // CSC
  for (conv = BUS_get_csc_conv(bus); conv !=NULL; conv = CONVCSC_get_next_ac(conv)) {
    
    // CSC converter with DC Power mode
    if (CONVCSC_is_in_P_dc_mode(conv) && CONVCSC_has_flags(conv,FLAG_VARS,CONVCSC_VAR_PDC)) {

      (*A_nnz) += 1;
      
      // Update rows
      (*A_row) += 1;
    }
  }
}

void CONSTR_CSC_DC_PSET_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* A_nnz;
  int* A_row;
  ConvCSC* conv;
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

  // CSC
  for (conv = BUS_get_csc_conv(bus); conv !=NULL; conv = CONVCSC_get_next_ac(conv)) {

    // CSC converter with DC Power mode
    if (CONVCSC_is_in_P_dc_mode(conv) && CONVCSC_has_flags(conv,FLAG_VARS,CONVCSC_VAR_PDC)) {
      
      VEC_set(b,*A_row,CONVCSC_get_P_dc_set(conv,t));
      
      // v
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,CONVCSC_get_index_P_dc(conv,t));
      MAT_set_d(A,*A_nnz,1.);	    
      (*A_nnz)++;

      // Update rows
      (*A_row)++;
    }
  }
}

void CONSTR_CSC_DC_PSET_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_CSC_DC_PSET_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
