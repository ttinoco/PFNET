/** @file constr_LOAD_VDEP.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LOAD_VDEP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LOAD_VDEP.h>

Constr* CONSTR_LOAD_VDEP_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_LOAD_VDEP_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_LOAD_VDEP_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_LOAD_VDEP_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_LOAD_VDEP_store_sens_step);
  CONSTR_set_name(c,"load voltage dependence");
  return c;
}

void CONSTR_LOAD_VDEP_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Load* load;
  int* J_nnz;
  int* J_row;
  int* H_nnz;

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!J_nnz || !J_row || !H_nnz || !bus)
    return;

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;
    
    // P and Q var
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P|LOAD_VAR_Q)) {
      
      // J
      (*J_nnz)++; // dSp/dp
      (*J_nnz)++; // dSq/dq
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
        
        // J
        (*J_nnz)++; // dSp/dv
        (*J_nnz)++; // dSq/dv
        
        // H
        (*(H_nnz+(*J_row)))++;   // dv and dv (Sp)
        (*(H_nnz+(*J_row+1)))++; // dv and dv (Sq) 
      }
      
      // Update rows
      (*J_row)++;
      (*J_row)++;
    }	
  }
}

void CONSTR_LOAD_VDEP_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Load* load;
  Mat* HSP;
  Mat* HSQ;
  Mat* J;
  Mat* H_array;
  int* J_nnz;
  int* J_row;
  int* H_nnz;

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!J || !J_nnz || !H_array || !J_row || !H_nnz || !bus)
    return;

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;
    
    // P and Q var
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P|LOAD_VAR_Q)) {
      
      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,LOAD_get_index_P(load,t));
      (*J_nnz)++; // dSp/dp
      
      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,LOAD_get_index_Q(load,t));
      (*J_nnz)++; // dSq/dq
      
      HSP = MAT_array_get(H_array,*J_row);
      HSQ = MAT_array_get(H_array,*J_row+1);
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
        
        // J
        MAT_set_i(J,*J_nnz,*J_row);
        MAT_set_j(J,*J_nnz,BUS_get_index_v_mag(bus,t));
        (*J_nnz)++; // dSp/dv
        
        MAT_set_i(J,*J_nnz,*J_row+1);
        MAT_set_j(J,*J_nnz,BUS_get_index_v_mag(bus,t));
        (*J_nnz)++; // dSq/dv
        
        // H
        MAT_set_i(HSP,*(H_nnz+(*J_row)),BUS_get_index_v_mag(bus,t));
        MAT_set_j(HSP,*(H_nnz+(*J_row)),BUS_get_index_v_mag(bus,t));
        (*(H_nnz+(*J_row)))++;   // dv and dv (Sp)
        
        MAT_set_i(HSQ,*(H_nnz+(*J_row+1)),BUS_get_index_v_mag(bus,t));
        MAT_set_j(HSQ,*(H_nnz+(*J_row+1)),BUS_get_index_v_mag(bus,t));
        (*(H_nnz+(*J_row+1)))++; // dv and dv (Sq) 
      }
      
      // Update rows
      (*J_row)++;
      (*J_row)++;
    }	
  }
}

void CONSTR_LOAD_VDEP_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Load* load;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* HSP;
  REAL* HSQ;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL P;
  REAL Q;
  REAL v;
  REAL cp;
  REAL cq;
  REAL ci;
  REAL cj;
  REAL cg;
  REAL cb;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!H_array || !f || !J || !J_nnz || !J_row || !H_nnz || !bus)
    return;

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;
	
    // P and Q var
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P|LOAD_VAR_Q)) {

      P = VEC_get(values,LOAD_get_index_P(load,t));
      Q = VEC_get(values,LOAD_get_index_Q(load,t));
      cp = LOAD_get_comp_cp(load,t);
      cq = LOAD_get_comp_cq(load,t);
      ci = LOAD_get_comp_ci(load,t);
      cj = LOAD_get_comp_cj(load,t);
      cg = LOAD_get_comp_cg(load);
      cb = LOAD_get_comp_cb(load);

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
        v = VEC_get(values,BUS_get_index_v_mag(bus,t));
      else
        v = BUS_get_v_mag(bus,t);

      HSP = MAT_get_data_array(MAT_array_get(H_array,*J_row));
      HSQ = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));

      // f
      f[*J_row] = P - cp - ci*v - cg*v*v;   // Sp
      f[*J_row+1] = Q - cq - cj*v + cb*v*v; // Sq

      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dSp/dp

      J[*J_nnz] = 1.;
      (*J_nnz)++; // dSq/dq

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

        // J
        J[*J_nnz] = -ci - 2.*cg*v;
        (*J_nnz)++; // dSp/dv

        J[*J_nnz] = -cj + 2.*cb*v;
        (*J_nnz)++; // dSq/dv

        // H
        HSP[*(H_nnz+(*J_row))] = -2.*cg;
        (*(H_nnz+(*J_row)))++;   // dv and dv (Sp)

        HSQ[*(H_nnz+(*J_row+1))] = 2.*cb;
        (*(H_nnz+(*J_row+1)))++; // dv and dv (Sq) 
      }

      // Update rows
      (*J_row)++;
      (*J_row)++;
    }	
  }      
}

void CONSTR_LOAD_VDEP_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
