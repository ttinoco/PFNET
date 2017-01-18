/** @file constr_PF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_PF.h>

struct Constr_PF_Data {

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

void CONSTR_PF_init(Constr* c) {

  // Local variables
  Net* net;
  int num_buses;
  int num_periods;
  Constr_PF_Data* data;

  // Init
  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  CONSTR_set_Hcounter(c,(int*)calloc(num_buses*num_periods,sizeof(int)),num_buses*num_periods);
  data = (Constr_PF_Data*)malloc(sizeof(Constr_PF_Data));
  data->size = num_buses*num_periods;
  ARRAY_zalloc(data->dPdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dQdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dPdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dQdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dwdw_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dwdv_indices,int,num_buses*num_periods);
  ARRAY_zalloc(data->dvdv_indices,int,num_buses*num_periods);
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_PF_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_Jcounter(c,0);
  CONSTR_clear_Hcounter(c);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_PF_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  int* Jcounter;
  int* Hcounter;
  int Hcounter_val;
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
  Constr_PF_Data* data;
  int k;
  int m;
  int num_buses;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_PF_Data*)CONSTR_get_data(c);

  // Check pointers
  if (!Jcounter || !Hcounter || !bus_counted)
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
      (*Jcounter)++; // dPm/dwk
      (*Jcounter)++; // dQm/dwk

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_w[m]) // wk and wm
	Hcounter_val++;
      if (var_v[m]) // wk and vm
	Hcounter_val++;
      if (var_a)    // wk and a
	Hcounter_val++;
      if (var_phi)  // wk and phi
	Hcounter_val++;
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //**********
    if (var_v[k]) { // vk var

      // J
      (*Jcounter)++; // dPm/dvk
      (*Jcounter)++; // dQm/dvk

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_w[m]) // vk and wm
	Hcounter_val++;
      if (var_v[m]) // vk and vm
	Hcounter_val++;
      if (var_a)    // vk and a
	Hcounter_val++;
      if (var_phi)  // vk and phi
	Hcounter_val++;
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      Hcounter_val++; // wm and wm
      if (var_v[m])   // wm and vm
	Hcounter_val++;
      if (var_a)      // wm and a
	Hcounter_val++;
      if (var_phi)    // wm and phi
	Hcounter_val++;
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      if (var_a)   // vm and a
	Hcounter[bus_index_t[k]]++;
      if (var_phi) // vm and phi
	Hcounter[bus_index_t[k]]++;
    }

    //********
    if (var_a) { // a var

      // J
      (*Jcounter)++; // dPk/da
      (*Jcounter)++; // dQk/da

      // H
      if (k == 0)  // a and a (important check k==0)
	Hcounter[bus_index_t[k]]++;
      if (var_phi) // a and phi
	Hcounter[bus_index_t[k]]++;
    }

    //**********
    if (var_phi) { // phi var

      // J
      (*Jcounter)++; // dPk/dphi
      (*Jcounter)++; // dQk/dphi

      // H
      Hcounter[bus_index_t[k]]++; // phi and phi
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      //***********
      if (var_w[k]) { // wk var

	// J
	dPdw_indices[bus_index_t[k]] = *Jcounter; // dPk/dwk
	(*Jcounter)++;
	dQdw_indices[bus_index_t[k]] = *Jcounter; // dQk/dwk
	(*Jcounter)++;

	// H
	Hcounter_val = Hcounter[bus_index_t[k]];
	dwdw_indices[bus_index_t[k]] = Hcounter_val; // wk and wk
	Hcounter_val++;
	if (var_v[k]) { // wk and vk
	  dwdv_indices[bus_index_t[k]] = Hcounter_val;
	  Hcounter_val++;
	}
	Hcounter[bus_index_t[k]] = Hcounter_val;
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	dPdv_indices[bus_index_t[k]] = *Jcounter; // dPk/dvk
	(*Jcounter)++;
	dQdv_indices[bus_index_t[k]] = *Jcounter; // dQk/dvk
	(*Jcounter)++;

	// H
	dvdv_indices[bus_index_t[k]] = Hcounter[bus_index_t[k]]; // vk and vk
	Hcounter[bus_index_t[k]]++;
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  (*Jcounter)++; // dPk/dPg
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  (*Jcounter)++; // dQk/dQg
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  (*Jcounter)++; // dPk/dPg
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  (*Jcounter)++; // dQk/dQg
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	//*****************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  (*Jcounter)++; // dQk/db

	  // H
	  if (var_v[k])
	    Hcounter[bus_index_t[k]]++; // b an vk
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_PF_allocate(Constr* c) {

  // Local variables
  Net* net;
  int num_buses;
  int num_periods;
  int num_constr;
  int num_vars;
  int Jcounter;
  int* Hcounter;
  int P_index;
  int Q_index;
  Mat* H_array;
  Mat* HP;
  Mat* HQ;
  int* row;
  int* col;
  int i;
  int t;
  int H_comb_nnz;
  int bus_index_t;

  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  num_vars = NET_get_num_vars(net);
  num_constr = 2*num_buses*num_periods;
  Jcounter = CONSTR_get_Jcounter(c);
  Hcounter = CONSTR_get_Hcounter(c);

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
			 Jcounter));  // nnz

  // H array
  H_comb_nnz = 0;
  H_array = MAT_array_new(num_constr);
  CONSTR_set_H_array(c,H_array,num_constr);
  for (t = 0; t < num_periods; t++) {
    for (i = 0; i < num_buses; i++) {
      bus_index_t = i+t*num_buses;
      P_index = BUS_get_index_P(NET_get_bus(net,i))+t*2*num_buses;
      Q_index = BUS_get_index_Q(NET_get_bus(net,i))+t*2*num_buses;
      HP = MAT_array_get(H_array,P_index);
      HQ = MAT_array_get(H_array,Q_index);
      MAT_set_nnz(HP,Hcounter[bus_index_t]);
      MAT_set_nnz(HQ,Hcounter[bus_index_t]);
      MAT_set_size1(HP,num_vars);
      MAT_set_size1(HQ,num_vars);
      MAT_set_size2(HP,num_vars);
      MAT_set_size2(HQ,num_vars);

      MAT_set_owns_rowcol(HP,TRUE);
      MAT_set_owns_rowcol(HQ,FALSE);

      ARRAY_zalloc(row,int,Hcounter[bus_index_t]);
      ARRAY_zalloc(col,int,Hcounter[bus_index_t]);

      MAT_set_row_array(HP,row); // same row array
      MAT_set_row_array(HQ,row);

      MAT_set_col_array(HP,col); // same col array
      MAT_set_col_array(HQ,col);

      MAT_set_data_array(HP,(REAL*)malloc(Hcounter[bus_index_t]*sizeof(REAL))); // different data array
      MAT_set_data_array(HQ,(REAL*)malloc(Hcounter[bus_index_t]*sizeof(REAL)));

      H_comb_nnz += 2*Hcounter[bus_index_t];
    }
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
}

void CONSTR_PF_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Mat* J;
  int* Jcounter;
  int* Hcounter;
  int Hcounter_val;
  int Hcounter_comb;
  char* bus_counted;
  Mat* H_array;
  Mat* H[2];
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
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
  int temp;
  int num_buses;
  int T;

  // Num buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Jcounter || !Hcounter || !bus_counted)
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
      MAT_set_i(J,*Jcounter,P_index[m]); // dPm/dwk
      MAT_set_j(J,*Jcounter,w_index[k]);
      (*Jcounter)++;

      MAT_set_i(J,*Jcounter,Q_index[m]); // dQm/dwk
      MAT_set_j(J,*Jcounter,w_index[k]);
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_w[m]) { // wk and wm
	MAT_set_i(H[k],Hcounter_val,w_index[k]);
	MAT_set_j(H[k],Hcounter_val,w_index[m]);
	Hcounter_val++;
      }
      if (var_v[m]) { // wk and vm
	MAT_set_i(H[k],Hcounter_val,w_index[k]);
	MAT_set_j(H[k],Hcounter_val,v_index[m]);
	Hcounter_val++;
      }
      if (var_a) {  // wk and a
	MAT_set_i(H[k],Hcounter_val,w_index[k]);
	MAT_set_j(H[k],Hcounter_val,a_index);
	Hcounter_val++;
      }
      if (var_phi) { // wk and phi
	MAT_set_i(H[k],Hcounter_val,w_index[k]);
	MAT_set_j(H[k],Hcounter_val,phi_index);
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_v[k]) { // vk var

      // J
      MAT_set_i(J,*Jcounter,P_index[m]); // dPm/dvk
      MAT_set_j(J,*Jcounter,v_index[k]);
      (*Jcounter)++;

      MAT_set_i(J,*Jcounter,Q_index[m]); // dQm/dvk
      MAT_set_j(J,*Jcounter,v_index[k]);
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_w[m]) { // vk and wm
	MAT_set_i(H[k],Hcounter_val,v_index[k]);
	MAT_set_j(H[k],Hcounter_val,w_index[m]);
	Hcounter_val++;
      }
      if (var_v[m]) { // vk and vm
	MAT_set_i(H[k],Hcounter_val,v_index[k]);
	MAT_set_j(H[k],Hcounter_val,v_index[m]);
	Hcounter_val++;
      }
      if (var_a) {   // vk and a
	MAT_set_i(H[k],Hcounter_val,v_index[k]);
	MAT_set_j(H[k],Hcounter_val,a_index);
	Hcounter_val++;
      }
      if (var_phi) { // vk and phi
	MAT_set_i(H[k],Hcounter_val,v_index[k]);
	MAT_set_j(H[k],Hcounter_val,phi_index);
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      MAT_set_i(H[k],Hcounter_val,w_index[m]); // wm and wm
      MAT_set_j(H[k],Hcounter_val,w_index[m]);
      Hcounter_val++;
      if (var_v[m]) {   // wm and vm
	MAT_set_i(H[k],Hcounter_val,w_index[m]);
	MAT_set_j(H[k],Hcounter_val,v_index[m]);
	Hcounter_val++;
      }
      if (var_a) {      // wm and a
	MAT_set_i(H[k],Hcounter_val,w_index[m]);
	MAT_set_j(H[k],Hcounter_val,a_index);
	Hcounter_val++;
      }
      if (var_phi) {    // wm and phi
	MAT_set_i(H[k],Hcounter_val,w_index[m]);
	MAT_set_j(H[k],Hcounter_val,phi_index);
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_a) {   // vm and a
	MAT_set_i(H[k],Hcounter_val,v_index[m]);
	MAT_set_j(H[k],Hcounter_val,a_index);
	Hcounter_val++;
      }
      if (var_phi) { // vm and phi
	MAT_set_i(H[k],Hcounter_val,v_index[m]);
	MAT_set_j(H[k],Hcounter_val,phi_index);
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //********
    if (var_a) { // a var

      // J
      MAT_set_i(J,*Jcounter,P_index[k]); // dPk/da
      MAT_set_j(J,*Jcounter,a_index);
      (*Jcounter)++;

      MAT_set_i(J,*Jcounter,Q_index[k]); // dQk/da
      MAT_set_j(J,*Jcounter,a_index);
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (k == 0) { // a and a (important check k==0)
	MAT_set_i(H[k],Hcounter_val,a_index);
	MAT_set_j(H[k],Hcounter_val,a_index);
	Hcounter_val++;
      }
      if (var_phi) { // a and phi
	MAT_set_i(H[k],Hcounter_val,a_index);
	MAT_set_j(H[k],Hcounter_val,phi_index);
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //**********
    if (var_phi) { // phi var

      // J
      MAT_set_i(J,*Jcounter,P_index[k]); // dPk/dphi
      MAT_set_j(J,*Jcounter,phi_index);
      (*Jcounter)++;

      MAT_set_i(J,*Jcounter,Q_index[k]); // dQk/dphi
      MAT_set_j(J,*Jcounter,phi_index);
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      MAT_set_i(H[k],Hcounter_val,phi_index);
      MAT_set_j(H[k],Hcounter_val,phi_index);
      Hcounter_val++; // phi and phi
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      //***********
      if (var_w[k]) { // wk var

	// J
	MAT_set_i(J,*Jcounter,P_index[k]); // dPk/dwk
	MAT_set_j(J,*Jcounter,w_index[k]);
	(*Jcounter)++;

	MAT_set_i(J,*Jcounter,Q_index[k]); // dQk/dwk
	MAT_set_j(J,*Jcounter,w_index[k]);
	(*Jcounter)++;

	// H
	Hcounter_val = Hcounter[bus_index_t[k]];
	MAT_set_i(H[k],Hcounter_val,w_index[k]); // wk and wk
	MAT_set_j(H[k],Hcounter_val,w_index[k]);
	Hcounter_val++;
	if (var_v[k]) { // wk and vk
	  MAT_set_i(H[k],Hcounter_val,w_index[k]);
	  MAT_set_j(H[k],Hcounter_val,v_index[k]);
	  Hcounter_val++;
	}
	Hcounter[bus_index_t[k]] = Hcounter_val;
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	MAT_set_i(J,*Jcounter,P_index[k]); // dPk/dvk
	MAT_set_j(J,*Jcounter,v_index[k]);
	(*Jcounter)++;

	MAT_set_i(J,*Jcounter,Q_index[k]); // dQk/dvk
	MAT_set_j(J,*Jcounter,v_index[k]);
	(*Jcounter)++;

	// H
	Hcounter_val = Hcounter[bus_index_t[k]];
	MAT_set_i(H[k],Hcounter_val,v_index[k]); // vk and vk
	MAT_set_j(H[k],Hcounter_val,v_index[k]);
	Hcounter_val++;
	Hcounter[bus_index_t[k]] = Hcounter_val;
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  MAT_set_i(J,*Jcounter,P_index[k]);             // dPk/dPg
	  MAT_set_j(J,*Jcounter,GEN_get_index_P(gen,t));
	  (*Jcounter)++;
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  MAT_set_i(J,*Jcounter,Q_index[k]);             // dQk/dQg
	  MAT_set_j(J,*Jcounter,GEN_get_index_Q(gen,t));
	  (*Jcounter)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  MAT_set_i(J,*Jcounter,P_index[k]);             // dPk/dPg
	  MAT_set_j(J,*Jcounter,VARGEN_get_index_P(vargen,t));
	  (*Jcounter)++;
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  MAT_set_i(J,*Jcounter,Q_index[k]);             // dQk/dQg
	  MAT_set_j(J,*Jcounter,VARGEN_get_index_Q(vargen,t));
	  (*Jcounter)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	//**************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  MAT_set_i(J,*Jcounter,Q_index[k]);         // dQk/db
	  MAT_set_j(J,*Jcounter,SHUNT_get_index_b(shunt,t));
	  (*Jcounter)++; // dQk/db

	  // H
	  if (var_v[k]) {
	    Hcounter_val = Hcounter[bus_index_t[k]];
	    MAT_set_i(H[k],Hcounter_val,SHUNT_get_index_b(shunt,t)); // b and vk
	    MAT_set_j(H[k],Hcounter_val,v_index[k]);
	    Hcounter_val++;
	    Hcounter[bus_index_t[k]] = Hcounter_val;
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }

  // Done
  if ((t == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {

    // Ensure lower triangular and save struct of H comb
    Hcounter_comb = 0;
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
	Hi_comb[Hcounter_comb] = Hi[m];
	Hj_comb[Hcounter_comb] = Hj[m];
	Hcounter_comb++;
      }
    }
  }
}

void CONSTR_PF_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Shunt* shunt;
  REAL* f;
  REAL* J;
  int* Jcounter;
  int* Hcounter;
  int Hcounter_val;
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

  REAL Pg;
  REAL Qg;

  REAL shunt_b;
  REAL shunt_g;

  Constr_PF_Data* data;

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
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_PF_Data*)CONSTR_get_data(c);

  // Check pointers
  if (!f || !J || !Jcounter || !Hcounter || !bus_counted || !data)
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
      w[k] = VEC_get(var_values,BUS_get_index_v_ang(bus[k],t));
    else
      w[k] = BUS_get_v_ang(bus[k],t);
    if (var_v[k])
      v[k] = VEC_get(var_values,BUS_get_index_v_mag(bus[k],t));
    else
      v[k] = BUS_get_v_mag(bus[k],t);
  }

  // Branch data
  var_a = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO);
  var_phi = BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE);
  b = BRANCH_get_b(br);           // series susceptance
  b_sh[0] = BRANCH_get_b_k(br);   //  total shunt susceptance on bus from (i)
  b_sh[1] = BRANCH_get_b_m(br);   //  total shunt susceptance on bus to (j)
  g = BRANCH_get_g(br);           // series conductance
  g_sh[0] = BRANCH_get_g_k(br);   //  total shunt conductance on bus from (i)
  g_sh[1] = BRANCH_get_g_m(br);   //  total shunt conductance on bus to (j)
  if (var_a)
    a = VEC_get(var_values,BRANCH_get_index_ratio(br,t));
  else
    a = BRANCH_get_ratio(br,t);
  if (var_phi)
    phi = VEC_get(var_values,BRANCH_get_index_phase(br,t));
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

    /* Branch flow equations for refernce:
    *  theta = w_k-w_m-theta_km+theta_mk
    *  P_km =  a_km^2*v_k^2*(g_km + gsh_km) - a_km*a_mk*v_k*v_m*( g_km*cos(theta) + b_km*sin(theta) )
    *  Q_km = -a_km^2*v_k^2*(b_km + bsh_km) - a_km*a_mk*v_k*v_m*( g_km*sin(theta) - b_km*cos(theta) )
    */

    // parts of the branch flow dependent on both vk, vm and the angles
    // (note that a == a_mk*a_km regardless of the direction since the other direction will always will always be 1)
    P_km[k] = -a*v[k]*v[m]*(g*cos(w[k]-w[m]-phi_temp)+b*sin(w[k]-w[m]-phi_temp));
    Q_km[k] = -a*v[k]*v[m]*(g*sin(w[k]-w[m]-phi_temp)-b*cos(w[k]-w[m]-phi_temp));

    // parts of the branch flow dependent on only vk^2
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
      J[*Jcounter] = -Q_km[m]; // dPm/dwk
      (*Jcounter)++;

      J[*Jcounter] = P_km[m];  // dQm/dwk
      (*Jcounter)++;

      J[data->dPdw_indices[bus_index_t[k]]] += Q_km[k];  // dPk/dwk
      J[data->dQdw_indices[bus_index_t[k]]] -= P_km[k]; // dQk/dwk

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      HP[k][data->dwdw_indices[bus_index_t[k]]] += P_km[k]; // wk and wk
      HQ[k][data->dwdw_indices[bus_index_t[k]]] += Q_km[k];
      if (var_v[k]) { // wk and vk
	HP[k][data->dwdv_indices[bus_index_t[k]]] += Q_km[k]/v[k]; // wk and wk
	HQ[k][data->dwdv_indices[bus_index_t[k]]] -= P_km[k]/v[k];
      }
      if (var_w[m]) { // wk and wm
	HP[k][Hcounter_val] = -P_km[k];
	HQ[k][Hcounter_val] = -Q_km[k];
	Hcounter_val++;
      }
      if (var_v[m]) { // wk and vm
	HP[k][Hcounter_val] = Q_km[k]/v[m];
	HQ[k][Hcounter_val] = -P_km[k]/v[m];
	Hcounter_val++;
      }
      if (var_a) {  // wk and a
	HP[k][Hcounter_val] = Q_km[k]/a;
	HQ[k][Hcounter_val] = -P_km[k]/a;
	Hcounter_val++;
      }
      if (var_phi) { // wk and phi
	HP[k][Hcounter_val] = -P_km[k]*indicator_phi;
	HQ[k][Hcounter_val] = -Q_km[k]*indicator_phi;
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //************
    if (var_v[k]) { // vk var

      // J
      J[*Jcounter] = -P_km[m]/v[k]; // dPm/dvk
      (*Jcounter)++;

      J[*Jcounter] = -Q_km[m]/v[k]; // dQm/dvk
      (*Jcounter)++;

      J[data->dPdv_indices[bus_index_t[k]]] -= 2*P_kk[k]/v[k] + P_km[k]/v[k]; // dPk/dvk
      J[data->dQdv_indices[bus_index_t[k]]] -= 2*Q_kk[k]/v[k] + Q_km[k]/v[k]; // dQk/dvk

      // H
      Hcounter_val = Hcounter[bus_index_T[k]];
      HP[k][data->dvdv_indices[bus_index_t[k]]] -= 2.*P_kk[k]/(v[k]*v[k]); // vk and vk
      HQ[k][data->dvdv_indices[bus_index_t[k]]] -= 2.*Q_kk[k]/(v[k]*v[k]);
      if (var_w[m]) { // vk and wm
	HP[k][Hcounter_val] = -Q_km[k]/v[k];
	HQ[k][Hcounter_val] = P_km[k]/v[k];
	Hcounter_val++;
      }
      if (var_v[m]) { // vk and vm
	HP[k][Hcounter_val] = -P_km[k]/(v[k]*v[m]);
	HQ[k][Hcounter_val] = -Q_km[k]/(v[k]*v[m]);
	Hcounter_val++;
      }
      if (var_a) {   // vk and a
	HP[k][Hcounter_val] = -indicator_a*P_kk[k]*4/(a*v[k]) - P_km[k]/(a*v[k]);
	HQ[k][Hcounter_val] = -indicator_a*Q_kk[k]*4/(a*v[k]) - Q_km[k]/(a*v[k]);
	Hcounter_val++;
      }
      if (var_phi) { // vk and phi
	HP[k][Hcounter_val] = -indicator_phi*Q_km[k]/v[k];
	HQ[k][Hcounter_val] = indicator_phi*P_km[k]/v[k];
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_w[m]) { // wm var

      // J
      // Nothing

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      HP[k][Hcounter_val] = P_km[k]; // wm and wm
      HQ[k][Hcounter_val] = Q_km[k];
      Hcounter_val++;
      if (var_v[m]) {   // wm and vm
	HP[k][Hcounter_val] = -Q_km[k]/v[m];
	HQ[k][Hcounter_val] = P_km[k]/v[m];
	Hcounter_val++;
      }
      if (var_a) {      // wm and a
	HP[k][Hcounter_val] = -Q_km[k]/a;
	HQ[k][Hcounter_val] = P_km[k]/a;
	Hcounter_val++;
      }
      if (var_phi) {    // wm and phi
	HP[k][Hcounter_val] = P_km[k]*indicator_phi;
	HQ[k][Hcounter_val] = Q_km[k]*indicator_phi;;
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //***********
    if (var_v[m]) { // vm var

      // J
      // Nothing

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (var_a) {   // vm and a
	HP[k][Hcounter_val] = -P_km[k]/(a*v[m]);
	HQ[k][Hcounter_val] = -Q_km[k]/(a*v[m]);
	Hcounter_val++;
      }
      if (var_phi) { // vm and phi
	HP[k][Hcounter_val] = -indicator_phi*Q_km[k]/v[m];
	HQ[k][Hcounter_val] = indicator_phi*P_km[k]/v[m];
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //********
    if (var_a) { // a var

      // J
      J[*Jcounter] = indicator_a*(-2.*P_kk[k]/a) - P_km[k]/a; // dPk/da
      (*Jcounter)++;

      J[*Jcounter] = indicator_a*(-2.*Q_kk[k]/a) - Q_km[k]/a; // dQk/da
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      if (k == 0) { // a and a (important check k==0)
	HP[k][Hcounter_val] = -P_kk[k]*2./(a*a);
	HQ[k][Hcounter_val] = -Q_kk[k]*2./(a*a);
	Hcounter_val++;
      }
      if (var_phi) { // a and phi
	HP[k][Hcounter_val] = -indicator_phi*Q_km[k]/a;
	HQ[k][Hcounter_val] = indicator_phi*P_km[k]/a;
	Hcounter_val++;
      }
      Hcounter[bus_index_t[k]] = Hcounter_val;
    }

    //**********
    if (var_phi) { // phi var

      // J
      J[*Jcounter] = -indicator_phi*Q_km[k]; // dPk/dphi
      (*Jcounter)++;

      J[*Jcounter] = indicator_phi*P_km[k]; // dQk/dphi
      (*Jcounter)++;

      // H
      Hcounter_val = Hcounter[bus_index_t[k]];
      HP[k][Hcounter_val] = P_km[k];
      HQ[k][Hcounter_val] = Q_km[k];
      Hcounter_val++; // phi and phi
      Hcounter[bus_index_t[k]] = Hcounter_val;
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
	(*Jcounter)++;

	// Nothing // dQk/dwk
	(*Jcounter)++;

	// H
	Hcounter[bus_index_t[k]]++;   // wk and wk
	if (var_v[k]) {
	  Hcounter[bus_index_t[k]]++; // wk and vk
	}
      }

      //***********
      if (var_v[k]) { // vk var

	// J
	// Nothing // dPk/dvk
	(*Jcounter)++;

	// Nothing // dQk/dvk
	(*Jcounter)++;

	// H
	Hcounter[bus_index_t[k]]++; // vk and vk
      }

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Var values
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  Pg = VEC_get(var_values,GEN_get_index_P(gen,t)); // p.u.
	else
	  Pg = GEN_get_P(gen,t);                           // p.u.
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	  Qg = VEC_get(var_values,GEN_get_index_Q(gen,t)); // p.u.
	else
	  Qg = GEN_get_Q(gen,t);                           // p.u.

	// f
	f[P_index[k]] += Pg;
	f[Q_index[k]] += Qg;

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var

	  // J
	  J[*Jcounter] = 1.; // dPk/dPg
	  (*Jcounter)++;
	}

	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	  // J
	  J[*Jcounter] = 1.; // dQk/dQg
	  (*Jcounter)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Var values
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
	  Pg = VEC_get(var_values,VARGEN_get_index_P(vargen,t)); // p.u.
	else
	  Pg = VARGEN_get_P(vargen,t);                           // p.u.
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q))
	  Qg = VEC_get(var_values,VARGEN_get_index_Q(vargen,t)); // p.u.
	else
	  Qg = VARGEN_get_Q(vargen,t);                           // p.u.

	// f
	f[P_index[k]] += Pg;
	f[Q_index[k]] += Qg;

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var

	  // J
	  J[*Jcounter] = 1.; // dPk/dPg
	  (*Jcounter)++;
	}

	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var

	  // J
	  J[*Jcounter] = 1.; // dQk/dQg
	  (*Jcounter)++;
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {

	// f
	f[P_index[k]] -= LOAD_get_P(load,t); // p.u.
	f[Q_index[k]] -= LOAD_get_Q(load,t); // p.u.
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {

	// f
	f[P_index[k]] -= BAT_get_P(bat,t); // p.u.
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Values
	shunt_g = SHUNT_get_g(shunt);
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC))
	  shunt_b = VEC_get(var_values,SHUNT_get_index_b(shunt,t)); // p.u.
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
	  J[*Jcounter] = v[k]*v[k]; // dQk/db
	  (*Jcounter)++;

	  // H
	  if (var_v[k]) {
	    Hcounter_val = Hcounter[bus_index_t[k]];
	    HP[k][Hcounter_val] = 0;
	    HQ[k][Hcounter_val] = 2*v[k];
	    Hcounter_val++; // b and vk
	    Hcounter[bus_index_t[k]] = Hcounter_val;
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_PF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

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

void CONSTR_PF_free(Constr* c) {

  // Local variables
  Constr_PF_Data* data;

  // Get data
  data = (Constr_PF_Data*)CONSTR_get_data(c);

  // Free
  if (data) {
    free(data->dPdw_indices);
    free(data->dQdw_indices);
    free(data->dPdv_indices);
    free(data->dQdv_indices);
    free(data->dwdw_indices);
    free(data->dwdv_indices);
    free(data->dvdv_indices);
  }
  free(data);

  // Set data
  CONSTR_set_data(c,NULL);
}
