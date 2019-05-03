/** @file constr_REG_TRAN.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_TRAN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_TRAN.h>

Constr* CONSTR_REG_TRAN_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_REG_TRAN_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_REG_TRAN_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_REG_TRAN_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_REG_TRAN_store_sens_step);
  CONSTR_set_name(c,"voltage regulation by transformers");
  return c;
}

void CONSTR_REG_TRAN_count_step(Constr* c, Bus* bus, BusDC* busdc, int tau) {
  
  // Local variables
  Branch* br;
  Bus* reg_bus;
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
  if (!A_nnz || !J_nnz || !A_row || !J_row || !H_nnz || !bus)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    if (BRANCH_is_tap_changer_v(br)) {

      reg_bus = BRANCH_get_reg_bus(br);
      
      // Linear
      //*******
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
        
        // A
        (*A_nnz)++; // t
      }
      
      (*A_nnz)++; // y
      (*A_nnz)++; // z
      
      (*A_row)++;    
      
      // Nonlinear constraints 1 (vmax,vmin)
      //************************************
      
      // J
      (*J_nnz)++; // dCompVmin/dy
      (*J_nnz)++; // dCompVmax/dz
      
      // H
      H_nnz[*J_row]++;   // y and y (vmin)
      H_nnz[*J_row+1]++; // z and z (vmax)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
        H_nnz[*J_row]++;   // y and v (vmin)
        H_nnz[*J_row+1]++; // z and v (vmax)
      }
      H_nnz[*J_row]++;   // y and vl (vmin)
      H_nnz[*J_row+1]++; // z and vh (vmax)
      
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	
        // J
        (*J_nnz)++; // dCompVmin/dv
        (*J_nnz)++; // dCompVmax/dv
	
        // H
        H_nnz[*J_row]++;   // v and v (vmin)
        H_nnz[*J_row+1]++; // v and v (vmax)
        H_nnz[*J_row]++;   // v and vl (vmin)
        H_nnz[*J_row+1]++; // v and vh (vmax)
      }
      
      // J
      (*J_nnz)++; // dCompVmin/dvl
      (*J_nnz)++; // dCompVmax/dvh
      
      // H 
      H_nnz[*J_row]++;   // vl and vl (vmin)
      H_nnz[*J_row+1]++; // vh and vh (vmax)
      
      // Nonlinear constraints 2 (tmax,tmin)
      //************************************
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
	
        // J
        (*J_nnz)++; // dCompTmax/dt
        (*J_nnz)++; // dCompTmin/dt
	
        // H
        H_nnz[*J_row+2]++; // t and t (tmax)
        H_nnz[*J_row+3]++; // t and t (tmin)
        H_nnz[*J_row+2]++; // t and vl (tmax)
        H_nnz[*J_row+3]++; // t and vh (tmin)
      }
      
      // J
      (*J_nnz)++; // dCompTmax/dvl
      (*J_nnz)++; // dCompTmin/dvh
      
      // H 
      H_nnz[*J_row+2]++; // vl and vl (tmax)
      H_nnz[*J_row+3]++; // vh and vh (tmin)
      
      // Count
      (*J_row)++; // CompVmin
      (*J_row)++; // CompVmax
      (*J_row)++; // CompTmax
      (*J_row)++; // CompTmin
      
      // Num extra vars
      CONSTR_set_num_extra_vars(c,*J_row);
    }
  }
}

void CONSTR_REG_TRAN_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int tau) {
  
  // Local variables
  Branch* br;
  Bus* reg_bus;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hvmin;
  Mat* Hvmax;
  Mat* Htmin;
  Mat* Htmax;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;
  int index_yz_vmin;
  int index_yz_vmax;
  int index_vvio_tmax;
  int index_vvio_tmin;
  int index_v;
  int index_vl;
  int index_vh;
  int index_t;
  int index_y;
  int index_z;
  int num_vars;

  // Number of vars
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

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
  if (!A_nnz || !J_nnz || !A_row || !H_array || !J_row || !H_nnz || !bus)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    if (BRANCH_is_tap_changer_v(br)) {
      
      reg_bus = BRANCH_get_reg_bus(br);
      
      // Hessians (NOTE ORDER!!!)
      Hvmin = MAT_array_get(H_array,*J_row);
      Hvmax = MAT_array_get(H_array,*J_row+1);
      Htmax = MAT_array_get(H_array,*J_row+2);
      Htmin = MAT_array_get(H_array,*J_row+3);
      
      // Indices
      if (BRANCH_has_pos_ratio_v_sens(br)) {
        index_yz_vmin = num_vars+(*J_row);
        index_yz_vmax = num_vars+(*J_row+1);
        index_vvio_tmax = num_vars+(*J_row+2);
        index_vvio_tmin = num_vars+(*J_row+3);
      }
      else {
        index_yz_vmin = num_vars+(*J_row+1);
        index_yz_vmax = num_vars+(*J_row);
        index_vvio_tmax = num_vars+(*J_row+3);
        index_vvio_tmin = num_vars+(*J_row+2);
      }
      index_v = BUS_get_index_v_mag(reg_bus,tau);
      index_y = num_vars+(*J_row);
      index_z = num_vars+(*J_row+1);
      index_vl = num_vars+(*J_row+2);
      index_vh = num_vars+(*J_row+3);
      index_t = BRANCH_get_index_ratio(br,tau);
      
      // Linear (t = t_0 + y - z)
      //*************************
      
      // b
      VEC_set(b,*A_row,BRANCH_get_ratio(br,tau)); // current ratio value
      
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
        
        // A
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,index_t);
        MAT_set_d(A,*A_nnz,1.);
        (*A_nnz)++; // t
      }
      else
        VEC_add_to_entry(b,*A_row,-BRANCH_get_ratio(br,tau));
      
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_y);
      MAT_set_d(A,*A_nnz,-1.);
      (*A_nnz)++; // y
      
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_z);
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++; // z
      
      (*A_row)++;
      
      // Nonlinear constraints 1 (vmin,vmax)
      //************************************
      
      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,index_yz_vmin);
      (*J_nnz)++; // dCompVmin/dy
      
      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,index_yz_vmax);
      (*J_nnz)++; // dCompVmax/dz
      
      // H	
      MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
      MAT_set_j(Hvmin,H_nnz[*J_row],index_yz_vmin);
      H_nnz[*J_row]++;   // y and y (vmin)
      
      MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
      MAT_set_j(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
      H_nnz[*J_row+1]++; // z and z (vmax)
      
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
        MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
        MAT_set_j(Hvmin,H_nnz[*J_row],index_v);
        H_nnz[*J_row]++;   // y and v (vmin)
	
        MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
        MAT_set_j(Hvmax,H_nnz[*J_row+1],index_v);
        H_nnz[*J_row+1]++; // z and v (vmax)
      }
      
      MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
      MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
      H_nnz[*J_row]++;   // y and vl (vmin)
      
      MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
      MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
      H_nnz[*J_row+1]++; // z and vh (vmax)
      
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	
        // J
        MAT_set_i(J,*J_nnz,*J_row);
        MAT_set_j(J,*J_nnz,index_v);
        (*J_nnz)++; // dCompVmin/dv
	
        MAT_set_i(J,*J_nnz,*J_row+1);
        MAT_set_j(J,*J_nnz,index_v);
        (*J_nnz)++; // dCompVmax/dv
	
        // H
        MAT_set_i(Hvmin,H_nnz[*J_row],index_v);
        MAT_set_j(Hvmin,H_nnz[*J_row],index_v);
        H_nnz[*J_row]++;   // v and v (vmin)
	
        MAT_set_i(Hvmax,H_nnz[*J_row+1],index_v);
        MAT_set_j(Hvmax,H_nnz[*J_row+1],index_v);
        H_nnz[*J_row+1]++; // v and v (vmax)
	
        MAT_set_i(Hvmin,H_nnz[*J_row],index_v);
        MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
        H_nnz[*J_row]++;   // v and vl (vmin)
	
        MAT_set_i(Hvmax,H_nnz[*J_row+1],index_v);
        MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
        H_nnz[*J_row+1]++; // v and vh (vmax)
      }
      
      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,index_vl);
      (*J_nnz)++; // dCompVmin/dvl
      
      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,index_vh);
      (*J_nnz)++; // dCompVmax/dvh
      
      // H 
      MAT_set_i(Hvmin,H_nnz[*J_row],index_vl);
      MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
      H_nnz[*J_row]++;   // vl and vl (vmin)
      
      MAT_set_i(Hvmax,H_nnz[*J_row+1],index_vh);
      MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
      H_nnz[*J_row+1]++; // vh and vh (vmax)
      
      // Nonlinear constraints 2 (tmax,tmin)
      //************************************
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
	
        // J
        MAT_set_i(J,*J_nnz,*J_row+2);
        MAT_set_j(J,*J_nnz,index_t);
        (*J_nnz)++; // dCompTmax/dt
	
        MAT_set_i(J,*J_nnz,*J_row+3);
        MAT_set_j(J,*J_nnz,index_t);
        (*J_nnz)++; // dCompTmin/dt
	
        // H
        MAT_set_i(Htmax,H_nnz[*J_row+2],index_t);
        MAT_set_j(Htmax,H_nnz[*J_row+2],index_t);
        H_nnz[*J_row+2]++; // t and t (tmax)
	
        MAT_set_i(Htmin,H_nnz[*J_row+3],index_t);
        MAT_set_j(Htmin,H_nnz[*J_row+3],index_t);
        H_nnz[*J_row+3]++; // t and t (tmin)
	
        MAT_set_i(Htmax,H_nnz[*J_row+2],index_t);
        MAT_set_j(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
        H_nnz[*J_row+2]++; // t and vl (tmax)
	
        MAT_set_i(Htmin,H_nnz[*J_row+3],index_t);
        MAT_set_j(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
        H_nnz[*J_row+3]++; // t and vh (tmin)
      }
      
      // J
      MAT_set_i(J,*J_nnz,*J_row+2);
      MAT_set_j(J,*J_nnz,index_vvio_tmax);
      (*J_nnz)++; // dCompTmax/dvl
      
      MAT_set_i(J,*J_nnz,*J_row+3);
      MAT_set_j(J,*J_nnz,index_vvio_tmin);
      (*J_nnz)++; // dCompTmin/dvh
      
      // H 
      MAT_set_i(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
      MAT_set_j(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
      H_nnz[*J_row+2]++; // vl and vl (tmax)
      
      MAT_set_i(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
      MAT_set_j(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
      H_nnz[*J_row+3]++; // vh and vh (tmin)
      
      // Extra var limits
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row,-CONSTR_REG_TRAN_MAX_YZ);      // y
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row+1,-CONSTR_REG_TRAN_MAX_YZ);    // z
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row+2,0.-CONSTR_REG_TRAN_MAX_VLH); // vl
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row+3,0.-CONSTR_REG_TRAN_MAX_VLH); // vh
      
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row,CONSTR_REG_TRAN_MAX_YZ);    // y
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row+1,CONSTR_REG_TRAN_MAX_YZ);  // z
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row+2,CONSTR_REG_TRAN_MAX_VLH); // vl
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row+3,CONSTR_REG_TRAN_MAX_VLH); // vh
      
      // Count
      (*J_row)++; // CompVmin
      (*J_row)++; // CompVmax
      (*J_row)++; // CompTmax
      (*J_row)++; // CompTmin
    }
  }
}

void CONSTR_REG_TRAN_eval_step(Constr* c, Bus* bus, BusDC* busdc, int tau, Vec* values, Vec* values_extra) {
  
  // Local variables
  Branch* br;
  Bus* reg_bus;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* Hvmin;
  REAL* Hvmax;
  REAL* Htmin;
  REAL* Htmax;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL v;
  REAL vl;
  REAL vh;
  REAL vmin;
  REAL vmax;
  REAL t;
  REAL tmax;
  REAL tmin;
  REAL y;
  REAL z;
  REAL yz_vmin;
  REAL yz_vmax;
  REAL vvio_tmax;
  REAL vvio_tmin;  
  REAL sqrtermVmin;
  REAL sqrtermVmax;
  REAL sqrtermTmax;
  REAL sqrtermTmin;
  REAL norm = CONSTR_REG_TRAN_NORM;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!f || !J || !J_nnz || !J_row || !H_nnz || !bus)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;  
    
    if (BRANCH_is_tap_changer_v(br)) {
      
      reg_bus = BRANCH_get_reg_bus(br);
      
      // v values
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG))
        v = VEC_get(values,BUS_get_index_v_mag(reg_bus,tau));
      else
        v = BUS_get_v_mag(reg_bus,tau);
      vmax = BUS_get_v_max_reg(reg_bus);
      vmin = BUS_get_v_min_reg(reg_bus);
      if (VEC_get_size(values_extra) > 0) {
        vl = VEC_get(values_extra,*J_row+2);
        vh = VEC_get(values_extra,*J_row+3);
      }
      else {
        vl = 0;
        vh = 0;
      }
      
      // t values
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO))
        t = VEC_get(values,BRANCH_get_index_ratio(br,tau));
      else
        t = BRANCH_get_ratio(br,tau);
      tmax = BRANCH_get_ratio_max(br);
      tmin = BRANCH_get_ratio_min(br);
      if (VEC_get_size(values_extra) > 0) {
        y = VEC_get(values_extra,*J_row);
        z = VEC_get(values_extra,*J_row+1);
      }
      else {
        y = 0;
        z = 0;
      }
      
      // Values that depend on sensitivity
      if (BRANCH_has_pos_ratio_v_sens(br)) {
        yz_vmin = y;
        yz_vmax = z;
        vvio_tmax = vl;
        vvio_tmin = vh;
      }
      else {
        yz_vmin = z;
        yz_vmax = y;
        vvio_tmax = vh;
        vvio_tmin = vl;
      }
      
      // Terms
      sqrtermVmin = sqrt( (v+vl-vmin)*(v+vl-vmin) + yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM );
      sqrtermVmax = sqrt( (vmax-v+vh)*(vmax-v+vh) + yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM );
      sqrtermTmax = sqrt( (tmax-t)*(tmax-t) + vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM );
      sqrtermTmin = sqrt( (t-tmin)*(t-tmin) + vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM );
      
      // Hessians (NOTE ORDER!!!)
      Hvmin = MAT_get_data_array(MAT_array_get(H_array,*J_row));
      Hvmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));
      Htmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+2));
      Htmin = MAT_get_data_array(MAT_array_get(H_array,*J_row+3));
      
      // f
      f[*J_row] = ((v+vl-vmin) + yz_vmin - sqrtermVmin)*norm;   // vmin
      f[*J_row+1] = ((vmax-v+vh) + yz_vmax - sqrtermVmax)*norm; // vmax
      f[*J_row+2] = ((tmax-t) + vvio_tmax - sqrtermTmax)*norm;  // tmax
      f[*J_row+3] = ((t-tmin) + vvio_tmin - sqrtermTmin)*norm;  // tmin
      
      // Nonlinear constraints 1 (vmin,vmax)
      //************************************
      
      // J
      J[*J_nnz] = (1.-yz_vmin/sqrtermVmin)*norm;
      (*J_nnz)++; // dcompVmin/dy
      
      J[*J_nnz] = (1.-yz_vmax/sqrtermVmax)*norm;
      (*J_nnz)++; // dcompVmax/dz
      
      // H	
      Hvmin[H_nnz[*J_row]] = -(((v+vl-vmin)*(v+vl-vmin)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // y and y (vmin)
      
      Hvmax[H_nnz[*J_row+1]] = -(((vmax-v+vh)*(vmax-v+vh)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // z and z (vmax)
      
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
        Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
        H_nnz[*J_row]++;   // y and v (vmin)
	
        Hvmax[H_nnz[*J_row+1]] = -((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
        H_nnz[*J_row+1]++; // z and v (vmax)
      }
      
      Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // y and vl (vmin)
      
      Hvmax[H_nnz[*J_row+1]] = ((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // z and vh (vmax)
      
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	
        // J
        J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
        (*J_nnz)++; // dcompVmin/dv
	
        J[*J_nnz] = -((1.-(vmax-v+vh)/sqrtermVmax))*norm;
        (*J_nnz)++; // dcompVmax/dv
	
        // H
        Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
        H_nnz[*J_row]++;   // v and v (vmin)
	
        Hvmax[H_nnz[*J_row+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
        H_nnz[*J_row+1]++; // v and v (vmax)
	
        Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
        H_nnz[*J_row]++;   // v and vl (vmin)
	
        Hvmax[H_nnz[*J_row+1]] = ((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
        H_nnz[*J_row+1]++; // v and vh (vmax)
      }
      
      // J
      J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
      (*J_nnz)++; // dcompVmin/dvl
      
      J[*J_nnz] = (1.-(vmax-v+vh)/sqrtermVmax)*norm;
      (*J_nnz)++; // dcompVmax/dvh
      
      // H 
      Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // vl and vl (vmin)
      
      Hvmax[H_nnz[*J_row+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // vh and vh (vmax)
      
      // Nonlinear constraints 2 (tmax,tmin)
      //************************************
      if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
	
        // J
        J[*J_nnz] = -(1.-(tmax-t)/sqrtermTmax)*norm;
        (*J_nnz)++; // dcompTmax/dt
	
        J[*J_nnz] = (1.-(t-tmin)/sqrtermTmin)*norm;
        (*J_nnz)++; // dcompTmin/dt
	
        // H
        Htmax[H_nnz[*J_row+2]] = -((vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
        H_nnz[*J_row+2]++; // t and t (tmax)
	
        Htmin[H_nnz[*J_row+3]] = -((vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
        H_nnz[*J_row+3]++; // t and t (tmin)
	
        Htmax[H_nnz[*J_row+2]] = -(vvio_tmax*(tmax-t)/pow(sqrtermTmax,3.))*norm;
        H_nnz[*J_row+2]++; // t and vl (tmax)
	
        Htmin[H_nnz[*J_row+3]] = (vvio_tmin*(t-tmin)/pow(sqrtermTmin,3.))*norm;
        H_nnz[*J_row+3]++; // t and vh (tmin)
      }
      
      // J
      J[*J_nnz] = (1.-vvio_tmax/sqrtermTmax)*norm;
      (*J_nnz)++; // dcompTmax/dvl
      
      J[*J_nnz] = (1.-vvio_tmin/sqrtermTmin)*norm;
      (*J_nnz)++; // dcompTmin/dvh
      
      // H 
      Htmax[H_nnz[*J_row+2]] = -(((tmax-t)*(tmax-t) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
      H_nnz[*J_row+2]++; // vl and vl (tmax)
      
      Htmin[H_nnz[*J_row+3]] = -(((t-tmin)*(t-tmin) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
      H_nnz[*J_row+3]++; // vh and vh (tmin)
      
      // Count
      (*J_row)++; // compVmin
      (*J_row)++; // compVmax
      (*J_row)++; // compTmax
      (*J_row)++; // compTmin
    }
  }
}

void CONSTR_REG_TRAN_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int tau, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Branch* br;
  Bus* reg_bus;
  int* J_row;
  REAL lamCompVmin;
  REAL lamCompVmax;
  REAL lamCompTmax;
  REAL lamCompTmin;
  
  // Constr data
  J_row = CONSTR_get_J_row_ptr(c);

  // Check pointer
  if (!J_row || !bus)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;
    
    if (BRANCH_is_tap_changer_v(br)) {
      
      reg_bus = BRANCH_get_reg_bus(br);
      
      lamCompVmin = VEC_get(sf,*J_row);
      (*J_row)++; // compVmin
      lamCompVmax = VEC_get(sf,*J_row);
      (*J_row)++; // compVmax
      lamCompTmax = VEC_get(sf,*J_row);
      (*J_row)++; // compTmax
      lamCompTmin = VEC_get(sf,*J_row);
      (*J_row)++; // compTmin
      
      if (fabs(lamCompVmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
        BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmin,tau);
      if (fabs(lamCompVmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
        BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmax,tau);
      if (fabs(lamCompTmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
        BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmax,tau);
      if (fabs(lamCompTmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
        BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmin,tau);
    }
  }
}
