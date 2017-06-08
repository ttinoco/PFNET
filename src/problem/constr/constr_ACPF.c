/** @file constr_ACPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type ACPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_ACPF.h>

struct Constr_ACPF_Data {

  int size;

  // Key indices for Jacobian
  int* dPdw_indices;
  int* dQdw_indices;
  int* dPdv_indices;
  int* dQdv_indices;

  // Key indices for Hessian
  int* dwdw_indices;
  int* dwdv_indices;
  int* dvdv_indices;
};

Constr* CONSTR_ACPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_ACPF_init);
  CONSTR_set_func_count_step(c, &CONSTR_ACPF_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_ACPF_allocate);
  CONSTR_set_func_clear(c, &CONSTR_ACPF_clear);
  CONSTR_set_func_analyze_step(c, &CONSTR_ACPF_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_ACPF_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_ACPF_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_ACPF_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_ACPF_init(Constr* c) {

  // Local variables
  Net* net;
  int num_buses;
  int num_periods;
  Constr_ACPF_Data* data;

  // Init
  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_buses*num_periods,sizeof(int)),num_buses*num_periods);
  data = (Constr_ACPF_Data*)malloc(sizeof(Constr_ACPF_Data));
  data->size = num_buses*num_periods;
  ARRAY_zalloc(data->dPdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dQdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dPdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dQdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dwdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dwdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dvdv_indices,int,num_buses*num_periods);
  CONSTR_set_name(c,"AC power balance");
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_ACPF_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_J_nnz(c,0);
  CONSTR_clear_H_nnz(c);
  
  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_ACPF_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  char* bus_counted;
  int* dPdw_indices;
  int* dQdw_indices;
  int* dPdv_indices;
  int* dQdv_indices;
  int* dwdw_indices;
  int* dwdv_indices;
  int* dvdv_indices;
  int bus_index_t[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  Constr_ACPF_Data* data;
  int k;
  int m;
  int num_buses;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_ACPF_Data*)CONSTR_get_data(c);

  // Check pointers
  if (!J_nnz || !H_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  // Key J indices
  dPdw_indices = data->dPdw_indices;
  dQdw_indices = data->dQdw_indices;
  dPdv_indices = data->dPdv_indices;
  dQdv_indices = data->dQdv_indices;

  // Key H indices
  dwdw_indices = data->dwdw_indices;
  dwdv_indices = data->dwdv_indices;
  dvdv_indices = data->dvdv_indices;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index(bus[k])+t*num_buses;
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
      (*J_nnz)++; // dPm/dwk
      (*J_nnz)++; // dQm/dwk

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_w[m]) // wk and wm
	H_nnz_val++;
      if (var_v[m]) // wk and vm
	H_nnz_val++;
      if (var_a)    // wk and a
	H_nnz_val++;
      if (var_phi)  // wk and phi
	H_nnz_val++;
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //**********
    if (var_v[k]) { // vk var

      // J
      (*J_nnz)++; // dPm/dvk
      (*J_nnz)++; // dQm/dvk

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_w[m]) // vk and wm
	H_nnz_val++;
      if (var_v[m]) // vk and vm
	H_nnz_val++;
      if (var_a)    // vk and a
	H_nnz_val++;
      if (var_phi)  // vk and phi
	H_nnz_val++;
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      H_nnz_val++; // wm and wm
      if (var_v[m])   // wm and vm
	H_nnz_val++;
      if (var_a)      // wm and a
	H_nnz_val++;
      if (var_phi)    // wm and phi
	H_nnz_val++;
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      if (var_a)   // vm and a
	H_nnz[bus_index_t[k]]++;
      if (var_phi) // vm and phi
	H_nnz[bus_index_t[k]]++;
    }

    //********
    if (var_a) { // a var

      // J
      (*J_nnz)++; // dPk/da
      (*J_nnz)++; // dQk/da

      // H
      if (k == 0)  // a and a (important check k==0)
	H_nnz[bus_index_t[k]]++;
      if (var_phi) // a and phi
	H_nnz[bus_index_t[k]]++;
    }

    //**********
    if (var_phi) { // phi var

      // J
      (*J_nnz)++; // dPk/dphi
      (*J_nnz)++; // dQk/dphi

      // H
      H_nnz[bus_index_t[k]]++; // phi and phi
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      //***********
      if (var_w[k]) { // wk var

	// J
	dPdw_indices[bus_index_t[k]] = *J_nnz; // dPk/dwk
	(*J_nnz)++;
	dQdw_indices[bus_index_t[k]] = *J_nnz; // dQk/dwk
	(*J_nnz)++;

	// H
	H_nnz_val = H_nnz[bus_index_t[k]];
	dwdw_indices[bus_index_t[k]] = H_nnz_val; // wk and wk
	H_nnz_val++;
	if (var_v[k]) { // wk and vk
	  dwdv_indices[bus_index_t[k]] = H_nnz_val;
	  H_nnz_val++;
	}
	H_nnz[bus_index_t[k]] = H_nnz_val;
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	dPdv_indices[bus_index_t[k]] = *J_nnz; // dPk/dvk
	(*J_nnz)++;
	dQdv_indices[bus_index_t[k]] = *J_nnz; // dQk/dvk
	(*J_nnz)++;

	// H
	dvdv_indices[bus_index_t[k]] = H_nnz[bus_index_t[k]]; // vk and vk
	H_nnz[bus_index_t[k]]++;
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  (*J_nnz)++; // dPk/dPg
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  (*J_nnz)++; // dQk/dQg
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  (*J_nnz)++; // dPk/dPg
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  (*J_nnz)++; // dQk/dQg
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	//*****************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  (*J_nnz)++; // dQk/db

	  // H
	  if (var_v[k])
	    H_nnz[bus_index_t[k]]++; // b an vk
	}
      }
      
      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var

	  // J
	  (*J_nnz)++; // dPk/dPl
	}

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var

	  // J
	  (*J_nnz)++; // dQk/dQl
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // Pc and Pd var

	  // J
	  (*J_nnz)++; // Pc
	  (*J_nnz)++; // Pd
	}
      }  
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_ACPF_allocate(Constr* c) {

  // Local variables
  Net* net;
  int num_buses;
  int num_periods;
  int num_constr;
  int num_vars;
  int J_nnz;
  int* H_nnz;
  int P_index;
  int Q_index;
  Mat* HP;
  Mat* HQ;
  int* row;
  int* col;
  int i;
  int t;
  int bus_index_t;

  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  num_vars = NET_get_num_vars(net);
  num_constr = 2*num_buses*num_periods;
  J_nnz = CONSTR_get_J_nnz(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // G l u
  CONSTR_set_G(c,MAT_new(0,num_vars,0));
  CONSTR_set_l(c,VEC_new(0));
  CONSTR_set_u(c,VEC_new(0));

  // f
  CONSTR_set_f(c,VEC_new(num_constr));

  // J
  CONSTR_set_J(c,MAT_new(num_constr,  // size1 (rows)
			 num_vars,    // size2 (cols)
			 J_nnz));  // nnz

  // H array
  CONSTR_allocate_H_array(c,num_constr);
  for (t = 0; t < num_periods; t++) {
    for (i = 0; i < num_buses; i++) {
      bus_index_t = i+t*num_buses;
      P_index = BUS_get_index_P(NET_get_bus(net,i))+t*2*num_buses;
      Q_index = BUS_get_index_Q(NET_get_bus(net,i))+t*2*num_buses;
      HP = CONSTR_get_H_single(c,P_index);
      HQ = CONSTR_get_H_single(c,Q_index);
      MAT_set_nnz(HP,H_nnz[bus_index_t]);
      MAT_set_nnz(HQ,H_nnz[bus_index_t]);
      MAT_set_size1(HP,num_vars);
      MAT_set_size1(HQ,num_vars);
      MAT_set_size2(HP,num_vars);
      MAT_set_size2(HQ,num_vars);

      MAT_set_owns_rowcol(HP,TRUE);
      MAT_set_owns_rowcol(HQ,FALSE);

      ARRAY_zalloc(row,int,H_nnz[bus_index_t]);
      ARRAY_zalloc(col,int,H_nnz[bus_index_t]);

      MAT_set_row_array(HP,row); // same row array
      MAT_set_row_array(HQ,row);

      MAT_set_col_array(HP,col); // same col array
      MAT_set_col_array(HQ,col);

      MAT_set_data_array(HP,(REAL*)malloc(H_nnz[bus_index_t]*sizeof(REAL))); // different data array
      MAT_set_data_array(HQ,(REAL*)malloc(H_nnz[bus_index_t]*sizeof(REAL)));
    }
  }
}

void CONSTR_ACPF_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  Mat* J;
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  char* bus_counted;
  Mat* H_array;
  Mat* H[2];
  int bus_index_t[2];
  int w_index[2];
  int v_index[2];
  int P_index[2];
  int Q_index[2];
  BOOL var_v[2];
  BOOL var_w[2];
  BOOL var_a;
  BOOL var_phi;
  int a_index;
  int phi_index;
  int k;
  int m;
  int num_buses;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !H_nnz || !H_array || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index(bus[k])+t*num_buses;
    P_index[k] = BUS_get_index_P(bus[k])+t*2*num_buses;
    Q_index[k] = BUS_get_index_Q(bus[k])+t*2*num_buses;
    w_index[k] = BUS_get_index_v_ang(bus[k],t);
    v_index[k] = BUS_get_index_v_mag(bus[k],t);
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    H[k] = MAT_array_get(H_array,P_index[k]);
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

    //***********
    if (var_w[k]) { // wk var

      // J
      MAT_set_i(J,*J_nnz,P_index[m]); // dPm/dwk
      MAT_set_j(J,*J_nnz,w_index[k]);
      (*J_nnz)++;

      MAT_set_i(J,*J_nnz,Q_index[m]); // dQm/dwk
      MAT_set_j(J,*J_nnz,w_index[k]);
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_w[m]) { // wk and wm
	MAT_set_i(H[k],H_nnz_val,w_index[k]);
	MAT_set_j(H[k],H_nnz_val,w_index[m]);
	H_nnz_val++;
      }
      if (var_v[m]) { // wk and vm
	MAT_set_i(H[k],H_nnz_val,w_index[k]);
	MAT_set_j(H[k],H_nnz_val,v_index[m]);
	H_nnz_val++;
      }
      if (var_a) {  // wk and a
	MAT_set_i(H[k],H_nnz_val,w_index[k]);
	MAT_set_j(H[k],H_nnz_val,a_index);
	H_nnz_val++;
      }
      if (var_phi) { // wk and phi
	MAT_set_i(H[k],H_nnz_val,w_index[k]);
	MAT_set_j(H[k],H_nnz_val,phi_index);
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_v[k]) { // vk var

      // J
      MAT_set_i(J,*J_nnz,P_index[m]); // dPm/dvk
      MAT_set_j(J,*J_nnz,v_index[k]);
      (*J_nnz)++;

      MAT_set_i(J,*J_nnz,Q_index[m]); // dQm/dvk
      MAT_set_j(J,*J_nnz,v_index[k]);
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_w[m]) { // vk and wm
	MAT_set_i(H[k],H_nnz_val,v_index[k]);
	MAT_set_j(H[k],H_nnz_val,w_index[m]);
	H_nnz_val++;
      }
      if (var_v[m]) { // vk and vm
	MAT_set_i(H[k],H_nnz_val,v_index[k]);
	MAT_set_j(H[k],H_nnz_val,v_index[m]);
	H_nnz_val++;
      }
      if (var_a) {   // vk and a
	MAT_set_i(H[k],H_nnz_val,v_index[k]);
	MAT_set_j(H[k],H_nnz_val,a_index);
	H_nnz_val++;
      }
      if (var_phi) { // vk and phi
	MAT_set_i(H[k],H_nnz_val,v_index[k]);
	MAT_set_j(H[k],H_nnz_val,phi_index);
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      MAT_set_i(H[k],H_nnz_val,w_index[m]); // wm and wm
      MAT_set_j(H[k],H_nnz_val,w_index[m]);
      H_nnz_val++;
      if (var_v[m]) {   // wm and vm
	MAT_set_i(H[k],H_nnz_val,w_index[m]);
	MAT_set_j(H[k],H_nnz_val,v_index[m]);
	H_nnz_val++;
      }
      if (var_a) {      // wm and a
	MAT_set_i(H[k],H_nnz_val,w_index[m]);
	MAT_set_j(H[k],H_nnz_val,a_index);
	H_nnz_val++;
      }
      if (var_phi) {    // wm and phi
	MAT_set_i(H[k],H_nnz_val,w_index[m]);
	MAT_set_j(H[k],H_nnz_val,phi_index);
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_a) {   // vm and a
	MAT_set_i(H[k],H_nnz_val,v_index[m]);
	MAT_set_j(H[k],H_nnz_val,a_index);
	H_nnz_val++;
      }
      if (var_phi) { // vm and phi
	MAT_set_i(H[k],H_nnz_val,v_index[m]);
	MAT_set_j(H[k],H_nnz_val,phi_index);
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //********
    if (var_a) { // a var

      // J
      MAT_set_i(J,*J_nnz,P_index[k]); // dPk/da
      MAT_set_j(J,*J_nnz,a_index);
      (*J_nnz)++;

      MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/da
      MAT_set_j(J,*J_nnz,a_index);
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (k == 0) { // a and a (important check k==0)
	MAT_set_i(H[k],H_nnz_val,a_index);
	MAT_set_j(H[k],H_nnz_val,a_index);
	H_nnz_val++;
      }
      if (var_phi) { // a and phi
	MAT_set_i(H[k],H_nnz_val,a_index);
	MAT_set_j(H[k],H_nnz_val,phi_index);
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //**********
    if (var_phi) { // phi var

      // J
      MAT_set_i(J,*J_nnz,P_index[k]); // dPk/dphi
      MAT_set_j(J,*J_nnz,phi_index);
      (*J_nnz)++;

      MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/dphi
      MAT_set_j(J,*J_nnz,phi_index);
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      MAT_set_i(H[k],H_nnz_val,phi_index);
      MAT_set_j(H[k],H_nnz_val,phi_index);
      H_nnz_val++; // phi and phi
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      //***********
      if (var_w[k]) { // wk var

	// J
	MAT_set_i(J,*J_nnz,P_index[k]); // dPk/dwk
	MAT_set_j(J,*J_nnz,w_index[k]);
	(*J_nnz)++;

	MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/dwk
	MAT_set_j(J,*J_nnz,w_index[k]);
	(*J_nnz)++;

	// H
	H_nnz_val = H_nnz[bus_index_t[k]];
	MAT_set_i(H[k],H_nnz_val,w_index[k]); // wk and wk
	MAT_set_j(H[k],H_nnz_val,w_index[k]);
	H_nnz_val++;
	if (var_v[k]) { // wk and vk
	  MAT_set_i(H[k],H_nnz_val,w_index[k]);
	  MAT_set_j(H[k],H_nnz_val,v_index[k]);
	  H_nnz_val++;
	}
	H_nnz[bus_index_t[k]] = H_nnz_val;
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	MAT_set_i(J,*J_nnz,P_index[k]); // dPk/dvk
	MAT_set_j(J,*J_nnz,v_index[k]);
	(*J_nnz)++;

	MAT_set_i(J,*J_nnz,Q_index[k]); // dQk/dvk
	MAT_set_j(J,*J_nnz,v_index[k]);
	(*J_nnz)++;

	// H
	H_nnz_val = H_nnz[bus_index_t[k]];
	MAT_set_i(H[k],H_nnz_val,v_index[k]); // vk and vk
	MAT_set_j(H[k],H_nnz_val,v_index[k]);
	H_nnz_val++;
	H_nnz[bus_index_t[k]] = H_nnz_val;
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  MAT_set_i(J,*J_nnz,P_index[k]);             // dPk/dPg
	  MAT_set_j(J,*J_nnz,GEN_get_index_P(gen,t));
	  (*J_nnz)++;
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  MAT_set_i(J,*J_nnz,Q_index[k]);             // dQk/dQg
	  MAT_set_j(J,*J_nnz,GEN_get_index_Q(gen,t));
	  (*J_nnz)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  MAT_set_i(J,*J_nnz,P_index[k]);             // dPk/dPg
	  MAT_set_j(J,*J_nnz,VARGEN_get_index_P(vargen,t));
	  (*J_nnz)++;
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  MAT_set_i(J,*J_nnz,Q_index[k]);             // dQk/dQg
	  MAT_set_j(J,*J_nnz,VARGEN_get_index_Q(vargen,t));
	  (*J_nnz)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	//**************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  MAT_set_i(J,*J_nnz,Q_index[k]);         // dQk/db
	  MAT_set_j(J,*J_nnz,SHUNT_get_index_b(shunt,t));
	  (*J_nnz)++;

	  // H
	  if (var_v[k]) {
	    H_nnz_val = H_nnz[bus_index_t[k]];
	    MAT_set_i(H[k],H_nnz_val,SHUNT_get_index_b(shunt,t)); // b and vk
	    MAT_set_j(H[k],H_nnz_val,v_index[k]);
	    H_nnz_val++;
	    H_nnz[bus_index_t[k]] = H_nnz_val;
	  }
	}
      }
      
      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var

	  // J
	  MAT_set_i(J,*J_nnz,P_index[k]);                // dPk/dPl
	  MAT_set_j(J,*J_nnz,LOAD_get_index_P(load,t));
	  (*J_nnz)++;
	}

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var

	  // J
	  MAT_set_i(J,*J_nnz,Q_index[k]);                // dQk/dQl
	  MAT_set_j(J,*J_nnz,LOAD_get_index_Q(load,t));
	  (*J_nnz)++;
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {  // Pc and Pd var

	  // J
	  MAT_set_i(J,*J_nnz,P_index[k]);              // Pc
	  MAT_set_j(J,*J_nnz,BAT_get_index_Pc(bat,t));
	  (*J_nnz)++;

	  MAT_set_i(J,*J_nnz,P_index[k]);              // Pd
	  MAT_set_j(J,*J_nnz,BAT_get_index_Pd(bat,t));
	  (*J_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_ACPF_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Shunt* shunt;
  REAL* f;
  REAL* J;
  int* J_nnz;
  int* H_nnz;
  int H_nnz_val;
  char* bus_counted;
  Mat* H_array;
  REAL* HP[2];
  REAL* HQ[2];

  int bus_index_t[2];
  BOOL var_v[2];
  BOOL var_w[2];
  int P_index[2];
  int Q_index[2];

  REAL w[2];
  REAL v[2];

  REAL a;
  REAL a_temp;
  REAL phi;
  REAL phi_temp;

  BOOL var_a;
  BOOL var_phi;

  REAL b;
  REAL b_sh[2];

  REAL g;
  REAL g_sh[2];

  REAL P_km[2];
  REAL P_kk[2];
  REAL Q_km[2];
  REAL Q_kk[2];

  REAL P;
  REAL Q;

  REAL shunt_b;
  REAL shunt_g;

  Constr_ACPF_Data* data;

  int k;
  int m;

  REAL indicator_a;
  REAL indicator_phi;

  int num_buses;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_ACPF_Data*)CONSTR_get_data(c);

  // Check pointers
  if (!f || !J || !J_nnz || !H_nnz || !bus_counted || !data)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index(bus[k])+t*num_buses;
    P_index[k] = BUS_get_index_P(bus[k])+t*2*num_buses; // index in f for active power mismatch
    Q_index[k] = BUS_get_index_Q(bus[k])+t*2*num_buses; // index in f for reactive power mismatch
    var_w[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG);
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    HP[k] = MAT_get_data_array(MAT_array_get(H_array,P_index[k]));
    HQ[k] = MAT_get_data_array(MAT_array_get(H_array,Q_index[k]));
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
  b = BRANCH_get_b(br);
  b_sh[0] = BRANCH_get_b_k(br);
  b_sh[1] = BRANCH_get_b_m(br);
  g = BRANCH_get_g(br);
  g_sh[0] = BRANCH_get_g_k(br);
  g_sh[1] = BRANCH_get_g_m(br);
  if (var_a)
    a = VEC_get(values,BRANCH_get_index_ratio(br,t));
  else
    a = BRANCH_get_ratio(br,t);
  if (var_phi)
    phi = VEC_get(values,BRANCH_get_index_phase(br,t));
  else
    phi = BRANCH_get_phase(br,t);

  // Branch flows
  for (k = 0; k < 2; k++) {

    if (k == 0) {
      m = 1;
      a_temp = a;
      phi_temp = phi;
    }
    else {
      m = 0;
      a_temp = 1;
      phi_temp = -phi;
    }

    /** Branch flow equations for refernce:
     *  theta = w_k-w_m-theta_km+theta_mk
     *  P_km =  a_km^2*v_k^2*(g_km + gsh_km) - a_km*a_mk*v_k*v_m*( g_km*cos(theta) + b_km*sin(theta) )
     *  Q_km = -a_km^2*v_k^2*(b_km + bsh_km) - a_km*a_mk*v_k*v_m*( g_km*sin(theta) - b_km*cos(theta) )
     */

    // Parts of the branch flow dependent on both vk, vm and the angles
    // (note that a == a_mk*a_km regardless of the direction since the other direction will always will always be 1)
    P_km[k] = -a*v[k]*v[m]*(g*cos(w[k]-w[m]-phi_temp)+b*sin(w[k]-w[m]-phi_temp));
    Q_km[k] = -a*v[k]*v[m]*(g*sin(w[k]-w[m]-phi_temp)-b*cos(w[k]-w[m]-phi_temp));

    // Parts of the branch flow dependent on only vk^2
    P_kk[k] =  a_temp*a_temp*(g_sh[k]+g)*v[k]*v[k];
    Q_kk[k] = -a_temp*a_temp*(b_sh[k]+b)*v[k]*v[k];
  }

  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    if (k == 0) {
      m = 1;
      indicator_a = 1.;
      indicator_phi = 1.;
    }
    else {
      m = 0;
      indicator_a = 0.;
      indicator_phi = -1.;
    }

    // f
    f[P_index[k]] -= P_kk[k] + P_km[k]; // Pk
    f[Q_index[k]] -= Q_kk[k] + Q_km[k]; // Qk

    //***********
    if (var_w[k]) { // wk var

      // J
      J[*J_nnz] = -Q_km[m]; // dPm/dwk
      (*J_nnz)++;

      J[*J_nnz] = P_km[m];  // dQm/dwk
      (*J_nnz)++;

      J[data->dPdw_indices[bus_index_t[k]]] += Q_km[k];  // dPk/dwk
      J[data->dQdw_indices[bus_index_t[k]]] -= P_km[k]; // dQk/dwk

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      HP[k][data->dwdw_indices[bus_index_t[k]]] += P_km[k]; // wk and wk
      HQ[k][data->dwdw_indices[bus_index_t[k]]] += Q_km[k];
      if (var_v[k]) { // wk and vk
	HP[k][data->dwdv_indices[bus_index_t[k]]] += Q_km[k]/v[k]; // wk and wk
	HQ[k][data->dwdv_indices[bus_index_t[k]]] -= P_km[k]/v[k];
      }
      if (var_w[m]) { // wk and wm
	HP[k][H_nnz_val] = -P_km[k];
	HQ[k][H_nnz_val] = -Q_km[k];
	H_nnz_val++;
      }
      if (var_v[m]) { // wk and vm
	HP[k][H_nnz_val] = Q_km[k]/v[m];
	HQ[k][H_nnz_val] = -P_km[k]/v[m];
	H_nnz_val++;
      }
      if (var_a) {  // wk and a
	HP[k][H_nnz_val] = Q_km[k]/a;
	HQ[k][H_nnz_val] = -P_km[k]/a;
	H_nnz_val++;
      }
      if (var_phi) { // wk and phi
	HP[k][H_nnz_val] = -P_km[k]*indicator_phi;
	HQ[k][H_nnz_val] = -Q_km[k]*indicator_phi;
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //************
    if (var_v[k]) { // vk var

      // J
      J[*J_nnz] = -P_km[m]/v[k]; // dPm/dvk
      (*J_nnz)++;

      J[*J_nnz] = -Q_km[m]/v[k]; // dQm/dvk
      (*J_nnz)++;

      J[data->dPdv_indices[bus_index_t[k]]] -= 2*P_kk[k]/v[k] + P_km[k]/v[k]; // dPk/dvk
      J[data->dQdv_indices[bus_index_t[k]]] -= 2*Q_kk[k]/v[k] + Q_km[k]/v[k]; // dQk/dvk

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      HP[k][data->dvdv_indices[bus_index_t[k]]] -= 2.*P_kk[k]/(v[k]*v[k]); // vk and vk
      HQ[k][data->dvdv_indices[bus_index_t[k]]] -= 2.*Q_kk[k]/(v[k]*v[k]);
      if (var_w[m]) { // vk and wm
	HP[k][H_nnz_val] = -Q_km[k]/v[k];
	HQ[k][H_nnz_val] = P_km[k]/v[k];
	H_nnz_val++;
      }
      if (var_v[m]) { // vk and vm
	HP[k][H_nnz_val] = -P_km[k]/(v[k]*v[m]);
	HQ[k][H_nnz_val] = -Q_km[k]/(v[k]*v[m]);
	H_nnz_val++;
      }
      if (var_a) {   // vk and a
	HP[k][H_nnz_val] = -indicator_a*P_kk[k]*4/(a*v[k]) - P_km[k]/(a*v[k]);
	HQ[k][H_nnz_val] = -indicator_a*Q_kk[k]*4/(a*v[k]) - Q_km[k]/(a*v[k]);
	H_nnz_val++;
      }
      if (var_phi) { // vk and phi
	HP[k][H_nnz_val] = -indicator_phi*Q_km[k]/v[k];
	HQ[k][H_nnz_val] = indicator_phi*P_km[k]/v[k];
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      HP[k][H_nnz_val] = P_km[k]; // wm and wm
      HQ[k][H_nnz_val] = Q_km[k];
      H_nnz_val++;
      if (var_v[m]) {   // wm and vm
	HP[k][H_nnz_val] = -Q_km[k]/v[m];
	HQ[k][H_nnz_val] = P_km[k]/v[m];
	H_nnz_val++;
      }
      if (var_a) {      // wm and a
	HP[k][H_nnz_val] = -Q_km[k]/a;
	HQ[k][H_nnz_val] = P_km[k]/a;
	H_nnz_val++;
      }
      if (var_phi) {    // wm and phi
	HP[k][H_nnz_val] = P_km[k]*indicator_phi;
	HQ[k][H_nnz_val] = Q_km[k]*indicator_phi;;
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (var_a) {   // vm and a
	HP[k][H_nnz_val] = -P_km[k]/(a*v[m]);
	HQ[k][H_nnz_val] = -Q_km[k]/(a*v[m]);
	H_nnz_val++;
      }
      if (var_phi) { // vm and phi
	HP[k][H_nnz_val] = -indicator_phi*Q_km[k]/v[m];
	HQ[k][H_nnz_val] = indicator_phi*P_km[k]/v[m];
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //********
    if (var_a) { // a var

      // J
      J[*J_nnz] = indicator_a*(-2.*P_kk[k]/a) - P_km[k]/a; // dPk/da
      (*J_nnz)++;

      J[*J_nnz] = indicator_a*(-2.*Q_kk[k]/a) - Q_km[k]/a; // dQk/da
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      if (k == 0) { // a and a (important check k==0)
	HP[k][H_nnz_val] = -P_kk[k]*2./(a*a);
	HQ[k][H_nnz_val] = -Q_kk[k]*2./(a*a);
	H_nnz_val++;
      }
      if (var_phi) { // a and phi
	HP[k][H_nnz_val] = -indicator_phi*Q_km[k]/a;
	HQ[k][H_nnz_val] = indicator_phi*P_km[k]/a;
	H_nnz_val++;
      }
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }

    //**********
    if (var_phi) { // phi var

      // J
      J[*J_nnz] = -indicator_phi*Q_km[k]; // dPk/dphi
      (*J_nnz)++;

      J[*J_nnz] = indicator_phi*P_km[k]; // dQk/dphi
      (*J_nnz)++;

      // H
      H_nnz_val = H_nnz[bus_index_t[k]];
      HP[k][H_nnz_val] = P_km[k];
      HQ[k][H_nnz_val] = Q_km[k];
      H_nnz_val++; // phi and phi
      H_nnz[bus_index_t[k]] = H_nnz_val;
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      //***********
      if (var_w[k]) { // wk var

	// J
	// Nothing // dPk/dwk
	(*J_nnz)++;

	// Nothing // dQk/dwk
	(*J_nnz)++;

	// H
	H_nnz[bus_index_t[k]]++;   // wk and wk
	if (var_v[k]) {
	  H_nnz[bus_index_t[k]]++; // wk and vk
	}
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	// Nothing // dPk/dvk
	(*J_nnz)++;

	// Nothing // dQk/dvk
	(*J_nnz)++;

	// H
	H_nnz[bus_index_t[k]]++; // vk and vk
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Var values
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  P = VEC_get(values,GEN_get_index_P(gen,t)); // p.u.
	else
	  P = GEN_get_P(gen,t);                       // p.u.
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	  Q = VEC_get(values,GEN_get_index_Q(gen,t)); // p.u.
	else
	  Q = GEN_get_Q(gen,t);                       // p.u.

	// f
	f[P_index[k]] += P;
	f[Q_index[k]] += Q;

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  J[*J_nnz] = 1.; // dPk/dPg
	  (*J_nnz)++;
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  J[*J_nnz] = 1.; // dQk/dQg
	  (*J_nnz)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Var values
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
	  P = VEC_get(values,VARGEN_get_index_P(vargen,t)); // p.u.
	else
	  P = VARGEN_get_P(vargen,t);                       // p.u.
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q))
	  Q = VEC_get(values,VARGEN_get_index_Q(vargen,t)); // p.u.
	else
	  Q = VARGEN_get_Q(vargen,t);                       // p.u.

	// f
	f[P_index[k]] += P;
	f[Q_index[k]] += Q;

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  J[*J_nnz] = 1.; // dPk/dPg
	  (*J_nnz)++;
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  J[*J_nnz] = 1.; // dQk/dQg
	  (*J_nnz)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Values
	shunt_g = SHUNT_get_g(shunt);
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC))
	  shunt_b = VEC_get(values,SHUNT_get_index_b(shunt,t)); // p.u.
	else
	  shunt_b = SHUNT_get_b(shunt,t);

	// f
	f[P_index[k]] -= shunt_g*v[k]*v[k]; // p.u.
	f[Q_index[k]] += shunt_b*v[k]*v[k];  // p.u.

	//***********
	if (var_v[k]) { // var v

	  // J
	  J[data->dPdv_indices[bus_index_t[k]]] -= 2*shunt_g*v[k]; // dPk/dvk
	  J[data->dQdv_indices[bus_index_t[k]]] += 2*shunt_b*v[k]; // dQk/dvk

	  // H
	  HP[k][data->dvdv_indices[bus_index_t[k]]] -= 2*shunt_g; // vk and vk
	  HQ[k][data->dvdv_indices[bus_index_t[k]]] += 2*shunt_b;
	}

	//**************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  J[*J_nnz] = v[k]*v[k]; // dQk/db
	  (*J_nnz)++;

	  // H
	  if (var_v[k]) {
	    H_nnz_val = H_nnz[bus_index_t[k]];
	    HP[k][H_nnz_val] = 0;
	    HQ[k][H_nnz_val] = 2*v[k];
	    H_nnz_val++; // b and vk
	    H_nnz[bus_index_t[k]] = H_nnz_val;
	  }
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	// Var values
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  P = VEC_get(values,LOAD_get_index_P(load,t)); // p.u.
	else
	  P = LOAD_get_P(load,t);                       // p.u.
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q))
	  Q = VEC_get(values,LOAD_get_index_Q(load,t)); // p.u.
	else
	  Q = LOAD_get_Q(load,t);                       // p.u.

	// f
	f[P_index[k]] -= P; // p.u.
	f[Q_index[k]] -= Q; // p.u.

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var

	  // J
	  J[*J_nnz] = -1.; // dPk/dPl
	  (*J_nnz)++;
	}

	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var

	  // J
	  J[*J_nnz] = -1.; // dQk/dQl
	  (*J_nnz)++;
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	// var values
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P))
	  P = VEC_get(values,BAT_get_index_Pc(bat,t))-VEC_get(values,BAT_get_index_Pd(bat,t)); // p.u.
	else
	  P = BAT_get_P(bat,t);                                                                // p.u.
	
	// f
	f[P_index[k]] -= P; // p.u.

	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {  // Pc and Pd var

	  // J
	  J[*J_nnz] = -1.; // Pc
	  (*J_nnz)++;

	  J[*J_nnz] = 1.; // Pd
	  (*J_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_ACPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* bus[2];
  char* bus_counted;
  int k;
  int num_buses;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Buses
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    if (!bus_counted[BUS_get_index(bus[k])+t*num_buses]) {
      BUS_set_sens_P_balance(bus[k],VEC_get(sf,BUS_get_index_P(bus[k])+t*2*num_buses),t); // sens of P balance
      BUS_set_sens_Q_balance(bus[k],VEC_get(sf,BUS_get_index_Q(bus[k])+t*2*num_buses),t); // sens of Q balance
    }
    bus_counted[BUS_get_index(bus[k])+t*num_buses] = TRUE;
  }
}

void CONSTR_ACPF_free(Constr* c) {

  // Local variables
  Constr_ACPF_Data* data;

  // Get data
  data = (Constr_ACPF_Data*)CONSTR_get_data(c);

  // Free
  if (data) {
    free(data->dPdw_indices);
    free(data->dQdw_indices);
    free(data->dPdv_indices);
    free(data->dQdv_indices);
    free(data->dwdw_indices);
    free(data->dwdv_indices);
    free(data->dvdv_indices);
    free(data);
  }

  // Set data
  CONSTR_set_data(c,NULL);
}
