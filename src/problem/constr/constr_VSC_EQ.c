/** @file constr_VSC_EQ.c
 *  @brief This file defines the data structure and routines associated with the constraint of type VSC_EQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_VSC_EQ.h>

Constr* CONSTR_VSC_EQ_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_VSC_EQ_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_VSC_EQ_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_VSC_EQ_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_VSC_EQ_store_sens_step);
  CONSTR_set_name(c, "VSC converter equations");
  return c;
}

void CONSTR_VSC_EQ_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* conv;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !J_row || !H_nnz || !busdc)
    return;

  // VSC
  for (conv = BUSDC_get_vsc_conv(busdc); conv != NULL; conv = CONVVSC_get_next_dc(conv)) {

    // Linear
    //*******
    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_P)) {
      
      // A
      (*A_nnz)++; // Pac
    }
    
    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {
      
      // A
      (*A_nnz)++; // Pdc
      (*A_nnz)++; // idc
    }
    
    // Count row
    (*A_row)++;
    
    // Nonlinear
    //**********
    
    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {
      
      // J
      (*J_nnz)++; // dP/dPdc
      (*J_nnz)++; // dP/didc
    }
    
    if (BUSDC_has_flags(busdc, FLAG_VARS, BUSDC_VAR_V)) {
      
      // J
      (*J_nnz)++; //v
      
      // H
      if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC))
        H_nnz[*J_row]++;     // v and idc (dP)
    }
    
    // Count
    (*J_row)++; // dP
  }  
}

void CONSTR_VSC_EQ_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* conv;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hp;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;

  // Constr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !H_array || !J_row || !H_nnz || !busdc)
    return;

  // VSC
  for (conv = BUSDC_get_vsc_conv(busdc); conv != NULL; conv = CONVVSC_get_next_dc(conv)) {

    // Linear
    //*******

    // b
    VEC_set(b, *A_row, -CONVVSC_get_loss_coeff_A(conv));

    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_P)) {

      // A
      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, CONVVSC_get_index_P(conv, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // Pac
    }

    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {

      // A
      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, CONVVSC_get_index_P_dc(conv, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // Pdc

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, CONVVSC_get_index_i_dc(conv, t));
      if (CONVVSC_get_P_dc_set(conv, t)<=0)
        MAT_set_d(A, *A_nnz, -CONVVSC_get_loss_coeff_B(conv));
      else
        MAT_set_d(A, *A_nnz, CONVVSC_get_loss_coeff_B(conv));
      (*A_nnz)++; // idc
    }

    // Count row
    (*A_row)++;

    // Nonlinear
    //**********

    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {

      // J
      MAT_set_i(J, *J_nnz, *J_row);
      MAT_set_j(J, *J_nnz, CONVVSC_get_index_P_dc(conv, t));
      (*J_nnz)++; // dP/dPdc

      MAT_set_i(J, *J_nnz, *J_row);
      MAT_set_j(J, *J_nnz, CONVVSC_get_index_i_dc(conv, t));
      (*J_nnz)++; // dP/didc
    }

    if (BUSDC_has_flags(busdc, FLAG_VARS, BUSDC_VAR_V)) {

      // J
      MAT_set_i(J, *J_nnz, *J_row);
      MAT_set_j(J, *J_nnz, BUSDC_get_index_v(busdc, t));
      (*J_nnz)++; // v

      // H
      if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {
        Hp = MAT_array_get(H_array, *J_row);
        MAT_set_i(Hp, H_nnz[*J_row], BUSDC_get_index_v(busdc, t));
        MAT_set_j(Hp, H_nnz[*J_row], CONVVSC_get_index_i_dc(conv, t));
        H_nnz[*J_row]++;     // v and idc (dP)
      }
    }

    // Count
    (*J_row)++; // dP
  }
}

void CONSTR_VSC_EQ_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {

  // Local variables
  ConvVSC* conv;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* Hp;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL Pdc;
  REAL idc;
  REAL v;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!f || !J || !J_nnz || !J_row || !H_nnz || !busdc)
    return;
  
  // VSC
  for (conv = BUSDC_get_vsc_conv(busdc); conv != NULL; conv = CONVVSC_get_next_dc(conv)) {

    // Nonlinear
    //**********
    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {
      Pdc = VEC_get(values, CONVVSC_get_index_P_dc(conv, t));
      idc = VEC_get(values, CONVVSC_get_index_i_dc(conv, t));
    }
    else {
      Pdc = CONVVSC_get_P_dc(conv, t);
      idc = CONVVSC_get_i_dc(conv, t);
    }

    if (BUSDC_has_flags(busdc, FLAG_VARS, BUSDC_VAR_V))
      v = VEC_get(values, BUSDC_get_index_v(busdc, t));
    else
      v = BUSDC_get_v(busdc, t);

    f[*J_row]  = Pdc - idc * v; // P

    if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {

      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dP/dPdc

      J[*J_nnz] = -v;
      (*J_nnz)++; // dP/dIdc
    }
	
    if (BUSDC_has_flags(busdc, FLAG_VARS, BUSDC_VAR_V)) {
	  
      // J
      J[*J_nnz] = -idc;
      (*J_nnz)++; //v

      // H
      if (CONVVSC_has_flags(conv, FLAG_VARS, CONVVSC_VAR_PDC)) {
        Hp = MAT_get_data_array(MAT_array_get(H_array, *J_row));
        Hp[H_nnz[*J_row]] = -1.;
        H_nnz[*J_row]++;     // v and Idc (dP)
      }
    }

    // Count
    (*J_row)++; // dP
  }
}

void CONSTR_VSC_EQ_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
