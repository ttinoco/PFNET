/** @file constr_AC_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_AC_FLOW_LIM.h>

void CONSTR_AC_FLOW_LIM_init(Constr* c) {
  
  // Local variables
  Net* net;
  int num_branches;
  int num_periods;
  int max_num_constr;

  // Init
  net = CONSTR_get_network(c);
  num_branches = NET_get_num_branches(net);
  num_periods = NET_get_num_periods(net);
  max_num_constr = 2*num_branches*num_periods;
  CONSTR_set_H_nnz(c,(int*)calloc(max_num_constr,sizeof(int)),max_num_constr);
  CONSTR_set_data(c,NULL);
}

void CONSTR_AC_FLOW_LIM_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_J_nnz(c,0);
  CONSTR_set_J_row(c,0);
  CONSTR_clear_H_nnz(c);
}

void CONSTR_AC_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  int* J_row;
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int k;
  int m;

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  J_row = CONSTR_get_J_row_ptr(c);
  
  // Check pointers
  if (!J_nnz || !H_nnz || !J_row)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
  }
  
  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);

  // Branch
  //*******
  
  for (k = 0; k < 2; k++) {
    
    if (k == 0)
      m = 1;
    else
      m = 0;

    //***********
    if (var_w[k]) { // wk var
      
      // J 
      (*J_nnz)++; // d|ikm|/dwk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // wk and wk
      if (var_v[k]) 
	H_nnz_val++; // wk and vk
      if (var_w[m]) 
	H_nnz_val++; // wk and wm
      if (var_v[m]) 
	H_nnz_val++; // wk and vm
      if (var_a)    
	H_nnz_val++; // wk and a
      if (var_phi)  
	H_nnz_val++; // wk and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[k]) { // vk var

      // J 
      (*J_nnz)++; // d|ikm|/dvk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // vk and vk
      if (var_w[m]) 
	H_nnz_val++; // vk and wm
      if (var_v[m]) 
	H_nnz_val++; // vk and vm
      if (var_a)    
	H_nnz_val++; // vk and a
      if (var_phi)  
	H_nnz_val++; // vk and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J 
      (*J_nnz)++; // d|ikm|/dwm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // wm and wm
      if (var_v[m]) 
	H_nnz_val++; // wm and vm
      if (var_a)    
	H_nnz_val++; // wm and a
      if (var_phi)  
	H_nnz_val++; // wm and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var
      
      // J 
      (*J_nnz)++; // d|ikm|/dvm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // vm and vm
      if (var_a)    
	H_nnz_val++; // vm and a
      if (var_phi)  
	H_nnz_val++; // vm and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //********
    if (var_a) { // a var
      
      // J 
      (*J_nnz)++; // d|ikm|/da
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // a and a
      if (var_phi)  
	H_nnz_val++; // a and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }
    
    //**********
    if (var_phi) { // phi var
      
      // J 
      (*J_nnz)++; // d|ikm|/dphi
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      H_nnz_val++;   // phi and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }
    
    // Constraint counter
    (*J_row)++;
    
    // Num local extra vars
    CONSTR_set_num_local_extra_vars(c,*J_row);
  }
}

void CONSTR_AC_FLOW_LIM_allocate(Constr* c) {
  
  // Local variables
  Net* net;
  int num_vars;
  int num_extra_vars;
  int J_nnz;
  int* H_nnz;
  int J_row;
  Mat* H_array;
  Mat* Hi;
  int H_comb_nnz;
  int* row;
  int* col;
  int i;

  // Data
  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);
  num_extra_vars = CONSTR_get_num_extra_vars(c);
  J_nnz = CONSTR_get_J_nnz(c);
  H_nnz = CONSTR_get_H_nnz(c);
  J_row = CONSTR_get_J_row(c);

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // G Gbar l u
  CONSTR_set_G(c,MAT_new(J_row,    // rows
			 num_vars, // columns
			 0));      // nnz
  CONSTR_set_Gbar(c,MAT_new(J_row,          // rows
			    num_extra_vars, // columns (total, not just local)
			    J_row));        // nnz
  CONSTR_set_l(c,VEC_new(J_row));
  CONSTR_set_u(c,VEC_new(J_row));

  // f J Jbar
  CONSTR_set_f(c,VEC_new(J_row));
  CONSTR_set_J(c,MAT_new(J_row,    // size1 (rows)
			 num_vars, // size2 (cols)
			 J_nnz));  // nnz
  CONSTR_set_Jbar(c,MAT_new(J_row,          // size1 (rows)
			    num_extra_vars, // size2 (cols) (total, not just local)
			    J_row));        // nnz

  // H
  H_comb_nnz = 0;
  H_array = MAT_array_new(J_row);
  CONSTR_set_H_array(c,H_array,J_row);
  for (i = 0; i < J_row; i++) {
    Hi = MAT_array_get(H_array,i);
    MAT_set_nnz(Hi,H_nnz[i]);
    MAT_set_size1(Hi,num_vars);
    MAT_set_size2(Hi,num_vars);
    MAT_set_owns_rowcol(Hi,TRUE);
    ARRAY_zalloc(row,int,H_nnz[i]);
    ARRAY_zalloc(col,int,H_nnz[i]);
    MAT_set_row_array(Hi,row);
    MAT_set_col_array(Hi,col);
    MAT_set_data_array(Hi,(REAL*)malloc(H_nnz[i]*sizeof(REAL)));
    H_comb_nnz += H_nnz[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
}

void CONSTR_AC_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  int* J_row;
  Mat* J;
  Mat* Jbar;
  Mat* Gbar;
  Mat* H_array;
  Mat* H; 
  Vec* l;
  Vec* u;
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int H_nnz_comb;
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int v_index[2];
  int w_index[2];
  int a_index;
  int phi_index;
  int offset;
  int temp;
  int k;
  int m;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J = CONSTR_get_J(c);
  Jbar = CONSTR_get_Jbar(c);
  Gbar = CONSTR_get_Gbar(c);
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  J_row = CONSTR_get_J_row_ptr(c);
  offset = CONSTR_get_local_extra_vars_offset(c);
 
  // Check pointers
  if (!J_nnz || !H_nnz || !J_row || !J || !Jbar || 
      !H_array || !Gbar || !l || !u)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    w_index[k] = BUS_get_index_v_ang(bus[k],t);
    v_index[k] = BUS_get_index_v_mag(bus[k],t);
  }
  
  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);
  a_index = BRANCH_get_index_ratio(br,t);
  phi_index = BRANCH_get_index_phase(br,t);

  // Branch
  //*******
  
  for (k = 0; k < 2; k++) {
    
    if (k == 0)
      m = 1;
    else
      m = 0;

    H = MAT_array_get(H_array,*J_row);

    //***********
    if (var_w[k]) { // wk var
      
      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,w_index[k]);
      (*J_nnz)++; // d|ikm|/dwk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,w_index[k]);
      MAT_set_j(H,H_nnz_val,w_index[k]);
      H_nnz_val++;   // wk and wk
      if (var_v[k]) {
	MAT_set_i(H,H_nnz_val,w_index[k]);
	MAT_set_j(H,H_nnz_val,v_index[k]);
	H_nnz_val++; // wk and vk
      }
      if (var_w[m]) {
	MAT_set_i(H,H_nnz_val,w_index[k]);
	MAT_set_j(H,H_nnz_val,w_index[m]);
	H_nnz_val++; // wk and wm
      }
      if (var_v[m]) {
	MAT_set_i(H,H_nnz_val,w_index[k]);
	MAT_set_j(H,H_nnz_val,v_index[m]);
	H_nnz_val++; // wk and vm
      }
      if (var_a) {
	MAT_set_i(H,H_nnz_val,w_index[k]);
	MAT_set_j(H,H_nnz_val,a_index);
	H_nnz_val++; // wk and a
      }
      if (var_phi) {
	MAT_set_i(H,H_nnz_val,w_index[k]);
	MAT_set_j(H,H_nnz_val,phi_index);
	H_nnz_val++; // wk and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[k]) { // vk var

      // J 
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,v_index[k]);
      (*J_nnz)++; // d|ikm|/dvk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,v_index[k]);
      MAT_set_j(H,H_nnz_val,v_index[k]);
      H_nnz_val++;   // vk and vk
      if (var_w[m]) {
	MAT_set_i(H,H_nnz_val,v_index[k]);
	MAT_set_j(H,H_nnz_val,w_index[m]);
	H_nnz_val++; // vk and wm
      }
      if (var_v[m]) { 
	MAT_set_i(H,H_nnz_val,v_index[k]);
	MAT_set_j(H,H_nnz_val,v_index[m]);
	H_nnz_val++; // vk and vm
      }
      if (var_a) {
	MAT_set_i(H,H_nnz_val,v_index[k]);
	MAT_set_j(H,H_nnz_val,a_index);
	H_nnz_val++; // vk and a
      }
      if (var_phi) {  
	MAT_set_i(H,H_nnz_val,v_index[k]);
	MAT_set_j(H,H_nnz_val,phi_index);
	H_nnz_val++; // vk and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J 
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,w_index[m]);
      (*J_nnz)++; // d|ikm|/dwm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,w_index[m]);
      MAT_set_j(H,H_nnz_val,w_index[m]);
      H_nnz_val++;   // wm and wm
      if (var_v[m]) {
	MAT_set_i(H,H_nnz_val,w_index[m]);
	MAT_set_j(H,H_nnz_val,v_index[m]);
	H_nnz_val++; // wm and vm
      }
      if (var_a) {
	MAT_set_i(H,H_nnz_val,w_index[m]);
	MAT_set_j(H,H_nnz_val,a_index);
	H_nnz_val++; // wm and a
      }
      if (var_phi) {
	MAT_set_i(H,H_nnz_val,w_index[m]);
	MAT_set_j(H,H_nnz_val,phi_index);
	H_nnz_val++; // wm and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var
      
      // J 
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,v_index[m]);
      (*J_nnz)++; // d|ikm|/dvm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,v_index[m]);
      MAT_set_j(H,H_nnz_val,v_index[m]);
      H_nnz_val++;   // vm and vm
      if (var_a) {
	MAT_set_i(H,H_nnz_val,v_index[m]);
	MAT_set_j(H,H_nnz_val,a_index);
	H_nnz_val++; // vm and a
      }
      if (var_phi) {
	MAT_set_i(H,H_nnz_val,v_index[m]);
	MAT_set_j(H,H_nnz_val,phi_index);
	H_nnz_val++; // vm and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //********
    if (var_a) { // a var
      
      // J 
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,a_index);
      (*J_nnz)++; // d|ikm|/da
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,a_index);
      MAT_set_j(H,H_nnz_val,a_index);
      H_nnz_val++;   // a and a
      if (var_phi) {
	MAT_set_i(H,H_nnz_val,a_index);
	MAT_set_j(H,H_nnz_val,phi_index);
	H_nnz_val++; // a and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }
    
    //**********
    if (var_phi) { // phi var
      
      // J 
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,phi_index);
      (*J_nnz)++; // d|ikm|/dphi
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      MAT_set_i(H,H_nnz_val,phi_index);
      MAT_set_j(H,H_nnz_val,phi_index);
      H_nnz_val++;   // phi and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    // Jbar
    MAT_set_i(Jbar,*J_row,*J_row);
    MAT_set_j(Jbar,*J_row,offset+*J_row);
    MAT_set_d(Jbar,*J_row,-1.);

    // Gbar, l, u
    MAT_set_i(Gbar,*J_row,*J_row);
    MAT_set_j(Gbar,*J_row,offset+*J_row);
    MAT_set_d(Gbar,*J_row,1.);
    VEC_set(l,*J_row,0.);
    VEC_set(u,*J_row,BRANCH_get_ratingA(br));
    
    // Constraint counter
    (*J_row)++;
  }

  // Done
  if ((t == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {

    // Ensure lower triangular and save struct of H comb
    H_nnz_comb = 0;
    Hi_comb = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hj_comb = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (k = 0; k < CONSTR_get_H_array_size(c); k++) {
      Hi = MAT_get_row_array(MAT_array_get(H_array,k));
      Hj = MAT_get_col_array(MAT_array_get(H_array,k));
      for (m = 0; m < MAT_get_nnz(MAT_array_get(H_array,k)); m++) {
	if (Hi[m] < Hj[m]) {
	  temp = Hi[m];
	  Hi[m] = Hj[m];
	  Hj[m] = temp;
	}
	Hi_comb[H_nnz_comb] = Hi[m];
	Hj_comb[H_nnz_comb] = Hj[m];
	H_nnz_comb++;
      }
    }
  }
}

void CONSTR_AC_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* extra_values) {
  // Nothing
}

void CONSTR_AC_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_AC_FLOW_LIM_free(Constr* c) {
  // Nothing
}
