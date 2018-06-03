/** @file constr_ACPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type ACPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_ACPF.h>

Constr* CONSTR_ACPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_ACPF_init);
  CONSTR_set_func_count_step(c,&CONSTR_ACPF_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_ACPF_allocate);
  CONSTR_set_func_clear(c,&CONSTR_ACPF_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_ACPF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_ACPF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_ACPF_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_ACPF_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_ACPF_init(Constr* c) {

  // Local variables
  Net* net;
  int num_buses;
  int num_periods;

  // Init
  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_buses*num_periods,sizeof(int)),num_buses*num_periods);
  CONSTR_set_name(c,"AC power balance");
  CONSTR_set_data(c,NULL);
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
  char* bus_counted;
  int bus_index_t[2];
  BOOL var_v[2];
  int k;

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !H_nnz || !bus_counted)
    return;
  
  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index_t(bus[k],t);
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
  }

  // Branch
  //*******

  BRANCH_power_flow_count(br,
			  J_nnz,
			  H_nnz+bus_index_t[0],
			  t,
			  TRUE);  // Pkm, Qkm
  BRANCH_power_flow_count(br,
			  J_nnz,
			  H_nnz+bus_index_t[1],
			  t,
			  FALSE); // Pmk, Qmk

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

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
	if (var_v[k]) { // vk var

	  // J
	  (*J_nnz)++; // dPk/dvk
	  (*J_nnz)++; // dQk/dvk

	  // H
	  (*(H_nnz+bus_index_t[k]))++; // vk an vk
	  
	}

	//*****************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  (*J_nnz)++; // dQk/db
	  
	  // H
	  if (var_v[k])
	    (*(H_nnz+bus_index_t[k]))++; // b an vk
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
      bus_index_t = BUS_get_index_t(NET_get_bus(net,i),t);
      P_index = BUS_get_index_P(NET_get_bus(net,i),t);
      Q_index = BUS_get_index_Q(NET_get_bus(net,i),t);
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
  char* bus_counted;
  Mat* H_array;
  Mat* H[2];
  int bus_index_t[2];
  int v_index[2];
  int P_index[2];
  int Q_index[2];
  BOOL var_v[2];
  int k;
  
  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !H_nnz || !H_array || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index_t(bus[k],t);
    P_index[k] = BUS_get_index_P(bus[k],t);
    Q_index[k] = BUS_get_index_Q(bus[k],t);
    v_index[k] = BUS_get_index_v_mag(bus[k],t);
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    H[k] = MAT_array_get(H_array,P_index[k]);
  }

  // Branch
  //*******

  BRANCH_power_flow_analyze(br,
			    J_nnz,
			    J,
			    P_index[0],
			    Q_index[0],
			    H_nnz+bus_index_t[0],
			    H[0],
			    t,
			    TRUE); // Pkm, Qkm

  BRANCH_power_flow_analyze(br,
			    J_nnz,
			    J,
			    P_index[1],
			    Q_index[1],
			    H_nnz+bus_index_t[1],
			    H[1],
			    t,
			    FALSE); // Pmk, Qmk

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

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

	//*****************************
	if (var_v[k]) { // vk var
	  
	  // J
	  MAT_set_i(J,*J_nnz,P_index[k]);
	  MAT_set_j(J,*J_nnz,v_index[k]);
	  (*J_nnz)++; // dPk/dvk

	  MAT_set_i(J,*J_nnz,Q_index[k]);
	  MAT_set_j(J,*J_nnz,v_index[k]);
	  (*J_nnz)++; // dPk/dvk

	  // H
	  MAT_set_i(H[k],*(H_nnz+bus_index_t[k]),v_index[k]);
	  MAT_set_j(H[k],*(H_nnz+bus_index_t[k]),v_index[k]);
	  (*(H_nnz+bus_index_t[k]))++; // vk an vk
	  
	}

	//**************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  MAT_set_i(J,*J_nnz,Q_index[k]);
	  MAT_set_j(J,*J_nnz,SHUNT_get_index_b(shunt,t));
	  (*J_nnz)++; // dQk/db

	  // H
	  if (var_v[k]) {
	    MAT_set_i(H[k],*(H_nnz+bus_index_t[k]),SHUNT_get_index_b(shunt,t));
	    MAT_set_j(H[k],*(H_nnz+bus_index_t[k]),v_index[k]);
	    (*(H_nnz+bus_index_t[k]))++; // b and vk 
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
  char* bus_counted;
  Mat* H_array;
  REAL* HP[2];
  REAL* HQ[2];

  int bus_index_t[2];
  BOOL var_v[2];
  int P_index[2];
  int Q_index[2];

  REAL v[2];
  
  REAL P;
  REAL Q;

  REAL shunt_b;
  REAL shunt_g;

  int k;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!f || !J || !J_nnz || !H_nnz || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index_t[k] = BUS_get_index_t(bus[k],t);
    P_index[k] = BUS_get_index_P(bus[k],t); // index in f for active power mismatch
    Q_index[k] = BUS_get_index_Q(bus[k],t); // index in f for reactive power mismatch
    var_v[k] = BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG);
    HP[k] = MAT_get_data_array(MAT_array_get(H_array,P_index[k]));
    HQ[k] = MAT_get_data_array(MAT_array_get(H_array,Q_index[k]));
    if (var_v[k])
      v[k] = VEC_get(values,BUS_get_index_v_mag(bus[k],t));
    else
      v[k] = BUS_get_v_mag(bus[k],t);
  }

  // Branch
  //*******

  BRANCH_power_flow_eval(br,
			 f+P_index[0],
			 f+Q_index[0],
			 J_nnz,
			 J,
			 H_nnz+bus_index_t[0],
			 HP[0],
			 HQ[0],
			 values,
			 -1.,   // flows leaving bus are negative
			 t,
			 TRUE); // Pkm, Qkm

  BRANCH_power_flow_eval(br,
			 f+P_index[1],
			 f+Q_index[1],
			 J_nnz,
			 J,
			 H_nnz+bus_index_t[1],
			 HP[1],
			 HQ[1],
			 values,
			 -1.,    // flows leaving bus are negative
			 t,
			 FALSE); // Pmk, Qmk
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

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

	//*****************************
	if (var_v[k]) { // var v

	  // J
	  J[*J_nnz] = -2*shunt_g*v[k];
	  (*J_nnz)++; // dPk/dvk
	  J[*J_nnz] = 2*shunt_b*v[k];
	  (*J_nnz)++; // dQk/dvk
	  
	  // H
	  HP[k][*(H_nnz+bus_index_t[k])] = -2*shunt_g;
	  HQ[k][*(H_nnz+bus_index_t[k])] = 2*shunt_b;
	  (*(H_nnz+bus_index_t[k]))++; // vk and vk
	}

	//*****************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  J[*J_nnz] = v[k]*v[k]; // dQk/db
	  (*J_nnz)++;

	  // H
	  if (var_v[k]) {
	    HP[k][*(H_nnz+bus_index_t[k])] = 0;
	    HQ[k][*(H_nnz+bus_index_t[k])] = 2*v[k];
	    (*(H_nnz+bus_index_t[k]))++; // b and vk
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
  
  // Constr data
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Buses
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    if (!bus_counted[BUS_get_index_t(bus[k],t)]) {
      BUS_set_sens_P_balance(bus[k],VEC_get(sf,BUS_get_index_P(bus[k],t)),t); // sens of P balance
      BUS_set_sens_Q_balance(bus[k],VEC_get(sf,BUS_get_index_Q(bus[k],t)),t); // sens of Q balance
    }
    bus_counted[BUS_get_index_t(bus[k],t)] = TRUE;
  }
}

void CONSTR_ACPF_free(Constr* c) {
  // Nothing
}
