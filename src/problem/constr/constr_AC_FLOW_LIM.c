/** @file constr_AC_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_AC_FLOW_LIM.h>

#define HESSIAN_VAL() -(R*dRdx + I*dIdx)*(R*dRdy + I*dIdy)/sqrterm3+(dRdy*dRdx+dIdy*dIdx+R*d2Rdydx+I*d2Idydx)/sqrterm

Constr* CONSTR_AC_FLOW_LIM_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_AC_FLOW_LIM_init);
  CONSTR_set_func_count_step(c, &CONSTR_AC_FLOW_LIM_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_AC_FLOW_LIM_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_AC_FLOW_LIM_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_AC_FLOW_LIM_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_AC_FLOW_LIM_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_AC_FLOW_LIM_init(Constr* c) {

  // Init
  CONSTR_set_name(c,"AC branch flow limits");
}

void CONSTR_AC_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  int* G_nnz;
  int* G_row;
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
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!G_nnz || !G_row || !J_nnz || !J_row || !H_nnz )
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Check zero rating
  if (BRANCH_get_ratingA(br) == 0.)
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

    //**********
    (*J_nnz)++;  // extra var
    
    // Nonlinear constraint counter
    (*J_row)++;

    // G constraint
    (*G_row)++;
    (*G_nnz)++;
    
    // Num extra vars
    CONSTR_set_num_extra_vars(c,*J_row);
  }
}

void CONSTR_AC_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Mat* J;
  int* G_row;
  int* G_nnz;
  int* J_row;
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  Mat* H_array;
  Mat* H;    
  Mat* G; 
  Vec* l;
  Vec* u;
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int v_index[2];
  int w_index[2];
  int a_index;
  int phi_index;
  int k;
  int m;
  int num_vars;
  
  // Num net vars
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // Constr data
  J = CONSTR_get_J(c);
  G = CONSTR_get_G(c);
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
 
  // Check pointers
  if (!G_nnz || !G_row || !J_nnz || !J_row || !H_nnz || !J || !H_array || !l || !u)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Check zero rating
  if (BRANCH_get_ratingA(br) == 0.)
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

    //**********

    // Extra vars
    VEC_set(CONSTR_get_l_extra_vars(c),*J_row,-BRANCH_get_ratingA(br));
    VEC_set(CONSTR_get_u_extra_vars(c),*J_row,BRANCH_get_ratingA(br));
    
    // J
    MAT_set_i(J,*J_nnz,*J_row);
    MAT_set_j(J,*J_nnz,num_vars+(*J_row));    
    (*J_nnz)++;  // extra var

    // G, l, u
    MAT_set_i(G,*G_nnz,*G_row);
    MAT_set_j(G,*G_nnz,num_vars+(*G_row));
    MAT_set_d(G,*G_nnz,1.);
    VEC_set(l,*G_row,-BRANCH_get_ratingA(br));
    VEC_set(u,*G_row,BRANCH_get_ratingA(br));

    // Row info
    CONSTR_set_J_row_info_string(c,
				 *J_row,
				 "branch",               // object
				 BRANCH_get_index(br),   // object id
				 (k == 0) ? "km" : "mk", // constraint info
				 t);                     // time
    CONSTR_set_G_row_info_string(c,
				 *G_row,
				 "branch",               // object
				 BRANCH_get_index(br),   // object id
				 (k == 0) ? "km" : "mk", // constraint info
				 t);                     // time
    
    // Nonlinear onstraint counter
    (*J_row)++;

    // G constraint counters
    (*G_row)++;
    (*G_nnz)++;
  }
}

void CONSTR_AC_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  
  // Local variables
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  int* J_row;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* H; 
  Bus* bus[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int k;
  int m;

  REAL extra_var;

  REAL w[2];
  REAL v[2];

  REAL a;
  REAL a_temp;
  REAL phi;
  REAL phi_temp;

  REAL b;
  REAL b_sh[2];
  
  REAL g;
  REAL g_sh[2];

  REAL R;
  REAL I;
  REAL dRdx;
  REAL dIdx;
  REAL dRdy;
  REAL dIdy;
  REAL d2Rdydx;
  REAL d2Idydx;
  REAL sqrterm;
  REAL sqrterm3;

  REAL costheta;
  REAL sintheta;

  REAL indicator_a;
  REAL indicator_phi;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  J_row = CONSTR_get_J_row_ptr(c);
 
  // Check pointers
  if (!J_nnz || !H_nnz || !J_row || !f || !J || !H_array)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Check zero rating
  if (BRANCH_get_ratingA(br) == 0.)
    return;
  
  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    if (var_w[k])
      w[k] = VEC_get(values,BUS_get_index_v_ang(bus[k],t));
    else
      w[k] = BUS_get_v_ang(bus[k],t);
    if (var_v[k])
      v[k] = VEC_get(values,BUS_get_index_v_mag(bus[k],t));
    else
      v[k] = BUS_get_v_mag(bus[k],t);
  }
  
  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);
  if (var_a)
    a = VEC_get(values,BRANCH_get_index_ratio(br,t));
  else
    a = BRANCH_get_ratio(br,t);
  if (var_phi)
    phi = VEC_get(values,BRANCH_get_index_phase(br,t));
  else
    phi = BRANCH_get_phase(br,t);
  b = BRANCH_get_b(br);
  b_sh[0] = BRANCH_get_b_k(br);
  b_sh[1] = BRANCH_get_b_m(br);
  g = BRANCH_get_g(br);
  g_sh[0] = BRANCH_get_g_k(br);
  g_sh[1] = BRANCH_get_g_m(br);

  // Branch
  //*******
  
  for (k = 0; k < 2; k++) {
    
    if (k == 0) {
      m = 1;
      a_temp = a;
      phi_temp = phi;
      indicator_a = 1.;
      indicator_phi = 1.;
    }
    else {
      m = 0;
      a_temp = 1;
      phi_temp = -phi;
      indicator_a = 0.;
      indicator_phi = -1.;
    }

    // Trigs
    costheta = cos(-w[k]+w[m]+phi_temp);
    sintheta = sin(-w[k]+w[m]+phi_temp);

    // |ikm| = |R + j I|
    R = a_temp*a_temp*(g_sh[k]+g)*v[k]-a*v[m]*(g*costheta-b*sintheta);
    I = a_temp*a_temp*(b_sh[k]+b)*v[k]-a*v[m]*(g*sintheta+b*costheta);
    sqrterm = sqrt(R*R+I*I+CONSTR_AC_FLOW_LIM_PARAM);
    sqrterm3 = sqrterm*sqrterm*sqrterm;
    
    H = MAT_get_data_array(MAT_array_get(H_array,*J_row));

    // f
    f[*J_row] = sqrterm;
    
    //***********
    if (var_w[k]) { // wk var
      
      dRdx = -a*v[m]*(g*sintheta+b*costheta);  // dRdwk
      dIdx = -a*v[m]*(-g*costheta+b*sintheta); // dIdwk 
	
      // J
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/dwk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];

      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = -a*v[m]*(-g*costheta+b*sintheta);
      d2Idydx = -a*v[m]*(-g*sintheta-b*costheta);
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // wk and wk

      if (var_v[k]) {
	dRdy = a_temp*a_temp*(g_sh[k]+g);
	dIdy = a_temp*a_temp*(b_sh[k]+b);
	d2Rdydx = 0;
	d2Idydx = 0;
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wk and vk
      }
      if (var_w[m]) {
	dRdy = -a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = -a*v[m]*(g*costheta-b*sintheta);
	d2Idydx = -a*v[m]*(g*sintheta+b*costheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wk and wm
      }
      if (var_v[m]) {
	dRdy = -a*(g*costheta-b*sintheta);
	dIdy = -a*(g*sintheta+b*costheta);
	d2Rdydx = -a*(g*sintheta+b*costheta);
	d2Idydx = -a*(-g*costheta+b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wk and vm
      }
      if (var_a) {
	dRdy = indicator_a*2.*a_temp*(g_sh[k]+g)*v[k]-v[m]*(g*costheta-b*sintheta);
	dIdy = indicator_a*2.*a_temp*(b_sh[k]+b)*v[k]-v[m]*(g*sintheta+b*costheta);
	d2Rdydx = -v[m]*(g*sintheta+b*costheta);
	d2Idydx = -v[m]*(-g*costheta+b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wk and a
      }
      if (var_phi) {
	dRdy = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Idydx = -indicator_phi*a*v[m]*(g*sintheta+b*costheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wk and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[k]) { // vk var

      dRdx = a_temp*a_temp*(g_sh[k]+g);
      dIdx = a_temp*a_temp*(b_sh[k]+b);

      // J 
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/dvk
      
      // H
      H_nnz_val = H_nnz[(*J_row)];
      
      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = 0;
      d2Idydx = 0;
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // vk and vk

      if (var_w[m]) {
	dRdy = -a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = 0;
	d2Idydx = 0;
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vk and wm
      }

      if (var_v[m]) { 
	dRdy = -a*(g*costheta-b*sintheta);
	dIdy = -a*(g*sintheta+b*costheta);
	d2Rdydx = 0;
	d2Idydx = 0;
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vk and vm
      }

      if (var_a) {
	dRdy = indicator_a*2.*a_temp*(g_sh[k]+g)*v[k]-v[m]*(g*costheta-b*sintheta);
	dIdy = indicator_a*2.*a_temp*(b_sh[k]+b)*v[k]-v[m]*(g*sintheta+b*costheta);
	d2Rdydx = indicator_a*2.*a_temp*(g_sh[k]+g);
	d2Idydx = indicator_a*2.*a_temp*(b_sh[k]+b);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vk and a
      }

      if (var_phi) {  
	dRdy = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = 0;
	d2Idydx = 0;
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vk and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      dRdx = -a*v[m]*(-g*sintheta-b*costheta);
      dIdx = -a*v[m]*(g*costheta-b*sintheta);
      
      // J 
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/dwm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];

      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = -a*v[m]*(-g*costheta+b*sintheta);
      d2Idydx = -a*v[m]*(-g*sintheta-b*costheta);
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // wm and wm

      if (var_v[m]) {
	dRdy = -a*(g*costheta-b*sintheta);
	dIdy = -a*(g*sintheta+b*costheta);
	d2Rdydx = -a*(-g*sintheta-b*costheta);
	d2Idydx = -a*(g*costheta-b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wm and vm
      }

      if (var_a) {
	dRdy = indicator_a*2.*a_temp*(g_sh[k]+g)*v[k]-v[m]*(g*costheta-b*sintheta);
	dIdy = indicator_a*2.*a_temp*(b_sh[k]+b)*v[k]-v[m]*(g*sintheta+b*costheta);
	d2Rdydx = -v[m]*(-g*sintheta-b*costheta);
	d2Idydx = -v[m]*(g*costheta-b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wm and a
      }

      if (var_phi) {
	dRdy = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = -indicator_phi*a*v[m]*(-g*costheta+b*sintheta);
	d2Idydx = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // wm and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var
      
      dRdx = -a*(g*costheta-b*sintheta);
      dIdx = -a*(g*sintheta+b*costheta);

      // J 
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/dvm
      
      // H
      H_nnz_val = H_nnz[(*J_row)];

      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = 0;
      d2Idydx = 0;
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // vm and vm

      if (var_a) {
	dRdy = indicator_a*2.*a_temp*(g_sh[k]+g)*v[k]-v[m]*(g*costheta-b*sintheta);
	dIdy = indicator_a*2.*a_temp*(b_sh[k]+b)*v[k]-v[m]*(g*sintheta+b*costheta);
	d2Rdydx = -(g*costheta-b*sintheta);
	d2Idydx = -(g*sintheta+b*costheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vm and a
      }

      if (var_phi) {
	dRdy = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = -indicator_phi*a*(-g*sintheta-b*costheta);
	d2Idydx = -indicator_phi*a*(g*costheta-b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // vm and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //********
    if (var_a) { // a var
      
      dRdx = indicator_a*2.*a_temp*(g_sh[k]+g)*v[k]-v[m]*(g*costheta-b*sintheta);
      dIdx = indicator_a*2.*a_temp*(b_sh[k]+b)*v[k]-v[m]*(g*sintheta+b*costheta);
      
      // J 
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/da
      
      // H
      H_nnz_val = H_nnz[(*J_row)];

      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = indicator_a*2.*(g_sh[k]+g)*v[k];
      d2Idydx = indicator_a*2.*(b_sh[k]+b)*v[k];
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // a and a

      if (var_phi) {
	dRdy = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
	dIdy = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
	d2Rdydx = -indicator_phi*v[m]*(-g*sintheta-b*costheta);
	d2Idydx = -indicator_phi*v[m]*(g*costheta-b*sintheta);
	H[H_nnz_val] = HESSIAN_VAL();
	H_nnz_val++; // a and phi
      }
      H_nnz[(*J_row)] = H_nnz_val;
    }
    
    //**********
    if (var_phi) { // phi var
     
      dRdx = -indicator_phi*a*v[m]*(-g*sintheta-b*costheta);
      dIdx = -indicator_phi*a*v[m]*(g*costheta-b*sintheta);
      
      // J 
      J[*J_nnz] = (R*dRdx + I*dIdx)/sqrterm;
      (*J_nnz)++; // d|ikm|/dphi
      
      // H
      H_nnz_val = H_nnz[(*J_row)];

      dRdy = dRdx;
      dIdy = dIdx;
      d2Rdydx = -a*v[m]*(-g*costheta+b*sintheta);
      d2Idydx = -a*v[m]*(-g*sintheta-b*costheta);
      H[H_nnz_val] = HESSIAN_VAL();
      H_nnz_val++;   // phi and phi
      H_nnz[(*J_row)] = H_nnz_val;
    }

    //**********

    if (VEC_get_size(values_extra) > 0)
      extra_var = VEC_get(values_extra,*J_row);
    else
      extra_var = 0;
   
    // f
    f[*J_row] -= extra_var;
 
    // J 
    J[*J_nnz] = -1.;
    (*J_nnz)++;      // extra var
    
    // Constraint counter
    (*J_row)++;
  }  
}

void CONSTR_AC_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  int* J_row;
  int k;
  REAL mu_km;

  // Constr data
  J_row = CONSTR_get_J_row_ptr(c);
  
  // Check pointers
  if (!J_row)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Check zero rating
  if (BRANCH_get_ratingA(br) == 0.)
    return;
    
  // Branch
  //*******
  
  for (k = 0; k < 2; k++) {

    mu_km = VEC_get(sGu,*J_row);
    
    // Sensitivity
    if (fabs(mu_km) > fabs(BRANCH_get_sens_i_mag_u_bound(br, t)))
      BRANCH_set_sens_i_mag_u_bound(br, mu_km, t);
    
    // Constraint counter
    (*J_row)++;
  }
}

void CONSTR_AC_FLOW_LIM_free(Constr* c) {
  // Nothing
}
