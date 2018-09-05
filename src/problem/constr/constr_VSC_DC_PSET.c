/** @file constr_VSC_DC_PSET.c
 *  @brief This file defines the data structure and routines associated with the constraint of type dc_vset.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_VSC_DC_PSET.h>

Constr* CONSTR_VSC_DC_PSET_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_VSC_DC_PSET_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_VSC_DC_PSET_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_VSC_DC_PSET_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_VSC_DC_PSET_store_sens_step);
  CONSTR_set_name(c,"VSC DC power control");
  return c;
}

void CONSTR_VSC_DC_PSET_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* conv;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;

  // VSC
  for (conv = BUS_get_vsc_conv(bus); conv !=NULL; conv = CONVVSC_get_next_ac(conv)) {
    
    // VSC converter with DC Power mode
    if (CONVVSC_is_in_P_dc_mode(conv)) {
      
      if (CONVVSC_has_flags(conv,FLAG_VARS,CONVVSC_VAR_P))
        (*A_nnz) += 1;
      
      // update rows
      (*A_row) += 1;
    }
  }
}

void CONSTR_VSC_DC_PSET_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* A_nnz;
  int* A_row;
  ConvVSC* conv;
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

  // VSC
  for (conv = BUS_get_vsc_conv(bus); conv !=NULL; conv = CONVVSC_get_next_ac(conv)) {

    // VSC converter with DC Power mode
    if (CONVVSC_is_in_P_dc_mode(conv)) {
      
      if (CONVVSC_has_flags(conv,FLAG_VARS,CONVVSC_VAR_P)) {
        
        VEC_set(b,*A_row,CONVVSC_get_P_dc_set(conv,t));
        
        // v
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,CONVVSC_get_index_P(conv,t));
        MAT_set_d(A,*A_nnz,-1.);	    
        (*A_nnz)++;
      }
      else
        VEC_set(b,*A_row,CONVVSC_get_P_dc_set(conv,t)-CONVVSC_get_P(conv,t));
      
      // update rows
      (*A_row)++;
    }
  }
}

void CONSTR_VSC_DC_PSET_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_VSC_DC_PSET_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
