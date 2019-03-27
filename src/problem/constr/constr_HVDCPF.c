/** @file constr_HVDCPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type HVDCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2019, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_HVDCPF.h>

Constr* CONSTR_HVDCPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_HVDCPF_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_HVDCPF_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_HVDCPF_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_HVDCPF_store_sens_step);
  CONSTR_set_name(c,"HVDC power balance");
  return c;
}

void CONSTR_HVDCPF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {
  
  // Local variables
  ConvVSC* vsc_conv;
  ConvCSC* csc_conv;
  BranchDC* br;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointers
  if (!A_nnz || !A_row || !busdc)
    return;

  // VSC converters
  for (vsc_conv = BUSDC_get_vsc_conv(busdc); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_dc(vsc_conv)) {
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_PDC)) { // P_dc var
      
      // A
      (*A_nnz)++;
    }
  }

  // CSC converters
  for (csc_conv = BUSDC_get_csc_conv(busdc); csc_conv != NULL; csc_conv = CONVCSC_get_next_dc(csc_conv)) {
    
    //*****************************
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_PDC)) { // P_dc var
      
      // A
      (*A_nnz)++;
    }
  }

  // Branches
  for(br = BUSDC_get_branch_k(busdc); br != NULL; br = BRANCHDC_get_next_k(br)) {
    
    //***********
    if (BUSDC_has_flags(BRANCHDC_get_bus_k(br),FLAG_VARS,BUSDC_VAR_V)) { // vk var

      // A
      (*A_nnz)++;
      (*A_nnz)++;
    }

    //***********
    if (BUSDC_has_flags(BRANCHDC_get_bus_m(br),FLAG_VARS,BUSDC_VAR_V)) { // vm var

      // A
      (*A_nnz)++;
      (*A_nnz)++;
    }
  }

  // Rows
  (*A_row)++;
}

void CONSTR_HVDCPF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  BranchDC* br;
  ConvVSC* vsc_conv;
  ConvCSC* csc_conv;
  Mat* A;
  Vec* b;
  int* A_nnz;
  int* A_row;
  REAL r;

  // Constr data
  A = CONSTR_get_A(c);
  b = CONSTR_get_b(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointers
  if (!A_nnz || !A_row || !busdc)
    return;

  // VSC
  for (vsc_conv = BUSDC_get_vsc_conv(busdc); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_dc(vsc_conv)) {
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_PDC)) { // P_dc var
      
      // A
      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(busdc,t));
      MAT_set_j(A,*A_nnz,CONVVSC_get_index_i_dc(vsc_conv,t)); // i_dc
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(b,BUSDC_get_index_t(busdc,t),-CONVVSC_get_i_dc(vsc_conv,t));
    }
  }

  // CSC
  for (csc_conv = BUSDC_get_csc_conv(busdc); csc_conv != NULL; csc_conv = CONVCSC_get_next_dc(csc_conv)) {
    
    //*****************************
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_PDC)) { // P_dc var
      
      // A
      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(busdc,t));
      MAT_set_j(A,*A_nnz,CONVCSC_get_index_i_dc(csc_conv,t)); // i_dc
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(b,BUSDC_get_index_t(busdc,t),-CONVCSC_get_i_dc(csc_conv,t));
    }
  }

  // Branches
  for(br = BUSDC_get_branch_k(busdc); br != NULL; br = BRANCHDC_get_next_k(br)) {
    
    // Branch data
    r = BRANCHDC_get_r(br); // p.u.
    if (r < CONSTR_HVDCPF_MINR)
      r = CONSTR_HVDCPF_MINR;
    
    //***********
    if (BUSDC_has_flags(BRANCHDC_get_bus_k(br),FLAG_VARS,BUSDC_VAR_V)) { // vk var
      
      // A
      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(BRANCHDC_get_bus_k(br),t)); // dik
      MAT_set_j(A,*A_nnz,BUSDC_get_index_v(BRANCHDC_get_bus_k(br),t)); // vk
      MAT_set_d(A,*A_nnz,-1./r);
      (*A_nnz)++;

      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(BRANCHDC_get_bus_m(br),t)); // dim
      MAT_set_j(A,*A_nnz,BUSDC_get_index_v(BRANCHDC_get_bus_k(br),t)); // vk
      MAT_set_d(A,*A_nnz,1./r);
      (*A_nnz)++;

      
    }
    else {

      // b
      VEC_add_to_entry(b,BUSDC_get_index_t(BRANCHDC_get_bus_k(br),t),BUSDC_get_v(BRANCHDC_get_bus_k(br),t)/r);
      VEC_add_to_entry(b,BUSDC_get_index_t(BRANCHDC_get_bus_m(br),t),-BUSDC_get_v(BRANCHDC_get_bus_k(br),t)/r);
    }

    //***********
    if (BUSDC_has_flags(BRANCHDC_get_bus_m(br),FLAG_VARS,BUSDC_VAR_V)) { // vm var

      // A
      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(BRANCHDC_get_bus_k(br),t)); // dik
      MAT_set_j(A,*A_nnz,BUSDC_get_index_v(BRANCHDC_get_bus_m(br),t)); // vm
      MAT_set_d(A,*A_nnz,1./r);
      (*A_nnz)++;

      MAT_set_i(A,*A_nnz,BUSDC_get_index_t(BRANCHDC_get_bus_m(br),t)); // dim
      MAT_set_j(A,*A_nnz,BUSDC_get_index_v(BRANCHDC_get_bus_m(br),t)); // vm
      MAT_set_d(A,*A_nnz,-1./r);
      (*A_nnz)++;
    }
    else {

      // b
      VEC_add_to_entry(b,BUSDC_get_index_t(BRANCHDC_get_bus_k(br),t),-BUSDC_get_v(BRANCHDC_get_bus_m(br),t)/r);
      VEC_add_to_entry(b,BUSDC_get_index_t(BRANCHDC_get_bus_m(br),t),BUSDC_get_v(BRANCHDC_get_bus_m(br),t)/r);
    }
  }

  // Rows
  (*A_row)++;
}

void CONSTR_HVDCPF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_HVDCPF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
