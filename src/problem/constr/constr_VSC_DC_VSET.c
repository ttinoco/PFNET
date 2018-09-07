/** @file constr_VSC_DC_VSET.c
 *  @brief This file defines the data structure and routines associated with the constraint of type dc_vset.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_VSC_DC_VSET.h>

Constr* CONSTR_VSC_DC_VSET_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_VSC_DC_VSET_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_VSC_DC_VSET_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_VSC_DC_VSET_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_VSC_DC_VSET_store_sens_step);
  CONSTR_set_name(c,"VSC DC voltage control");
  return c;
}

void CONSTR_VSC_DC_VSET_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* conv;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !busdc)
    return;

  // VSC
  for (conv = BUSDC_get_vsc_conv(busdc); conv !=NULL; conv = CONVVSC_get_next_dc(conv)) {
    
    // VSC converter with DC voltage mode
    if (CONVVSC_is_in_v_dc_mode(conv) && BUSDC_has_flags(busdc,FLAG_VARS,BUSDC_VAR_V)) {
      
      (*A_nnz) += 1;
      
      // Update rows
      (*A_row) += 1;
    }
  }
}

void CONSTR_VSC_DC_VSET_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

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
  if (!A_nnz || !A_row || !busdc)
    return;

  // VSC
  for (conv = BUSDC_get_vsc_conv(busdc); conv !=NULL; conv = CONVVSC_get_next_dc(conv)) {

    // VSC converter with DC voltage mode
    if (CONVVSC_is_in_v_dc_mode(conv) && BUSDC_has_flags(busdc,FLAG_VARS,BUSDC_VAR_V)) {
        
      VEC_set(b,*A_row,CONVVSC_get_v_dc_set(conv,t));
      
      // v
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,BUSDC_get_index_v(busdc,t));
      MAT_set_d(A,*A_nnz,1.);	    
      (*A_nnz)++;
      
      // Update rows
      (*A_row)++;
    }
  }
}

void CONSTR_VSC_DC_VSET_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_VSC_DC_VSET_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
