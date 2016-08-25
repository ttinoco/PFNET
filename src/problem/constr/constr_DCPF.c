/** @file constr_DCPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_DCPF.h>

void CONSTR_DCPF_init(Constr* c) {
    
  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_DCPF_clear(Constr* c) {
    
  // Counters
  CONSTR_set_Acounter(c,0);
   
  // Flags
  CONSTR_clear_bus_counted(c);  
}

void CONSTR_DCPF_count_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* bus[2];
  Gen* gen;
  Load* load;
  Vargen* vargen;
  Bat* bat;
  int* Acounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  int m;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Acounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
 
  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(bus[k]);
  
  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    if (k == 0)
      m = 1;
    else
      m = 0;

    //***********
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A
      (*Acounter)++; // Pk
    }

    //***********
    if (BUS_has_flags(bus[m],FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A
      (*Acounter)++; // Pk
    }
    
    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // A
      (*Acounter)++; // Pk
    }
  }
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
	  
	  // A
	  (*Acounter)++; // Pk
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {
	
	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // P var
	  
	  // A
	  (*Acounter)++; // Pk
	}
      }
      
      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // P var
	  
	  // A
	  (*Acounter)++; // Pk
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {
	
	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var
	  
	  // A
	  (*Acounter)++; // Pc
	  (*Acounter)++; // Pd
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]*T+t] = TRUE;
  }
}

void CONSTR_DCPF_allocate(Constr* c) {
  
  // Local variables
  Net* net;
  int num_buses;
  int num_vars;
  int Acounter;

  net = CONSTR_get_network(c);
  num_buses = NET_get_num_buses(net);
  num_vars = NET_get_num_vars(net);
  Acounter = CONSTR_get_Acounter(c);
  
  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));
  
  // b
  CONSTR_set_b(c,VEC_new(num_buses*NET_get_num_periods(net)));

  // A
  CONSTR_set_A(c,MAT_new(num_buses*NET_get_num_periods(net), // size1 (rows)
			 num_vars,                           // size2 (cols)
			 Acounter));                         // nnz
}

void CONSTR_DCPF_analyze_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Mat* A;
  Vec* rhs;
  int* Acounter;
  char* bus_counted;
  int bus_index[2];
  REAL b;
  REAL sign_phi;
  int k;
  int m;
  int num_buses;
  int T;
  int Pmis_index[2];

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Number of buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));
 
  // Constr data
  A = CONSTR_get_A(c);
  rhs = CONSTR_get_b(c);
  Acounter = CONSTR_get_Acounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Acounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
 
  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++) {
    bus_index[k] = BUS_get_index(bus[k]);
    Pmis_index[k] = bus_index[k]+t*num_buses;
  }
 
  // Branch data
  b = BRANCH_get_b(br);

  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    if (k == 0) {
      m = 1;
      sign_phi = 1;
    }
    else {
      m = 0;
      sign_phi = -1;
    }

    //***********
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A 
      MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
      MAT_set_j(A,*Acounter,BUS_get_index_v_ang(bus[k],t)); // wk
      MAT_set_d(A,*Acounter,b); 
      (*Acounter)++;
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,Pmis_index[k],-b*BUS_get_v_ang(bus[k],t));
    }

    //***********
    if (BUS_has_flags(bus[m],FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A 
      MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
      MAT_set_j(A,*Acounter,BUS_get_index_v_ang(bus[m],t)); // wk
      MAT_set_d(A,*Acounter,-b);
      (*Acounter)++;
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,Pmis_index[k],b*BUS_get_v_ang(bus[m],t));
    }

    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var

      // A
      MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
      MAT_set_j(A,*Acounter,BRANCH_get_index_phase(br,t)); // phi
      MAT_set_d(A,*Acounter,-b*sign_phi);
      (*Acounter)++; 
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,Pmis_index[k],b*BRANCH_get_phase(br,t)*sign_phi);
    }
  }
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
	  
	  // A
	  MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
	  MAT_set_j(A,*Acounter,GEN_get_index_P(gen,t)); // Pg
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,Pmis_index[k],-GEN_get_P(gen,t));
	}
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {
	
	//*****************************
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var
	  
	  // A
	  MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
	  MAT_set_j(A,*Acounter,LOAD_get_index_P(load,t)); // Pl
	  MAT_set_d(A,*Acounter,-1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,Pmis_index[k],LOAD_get_P(load,t));
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
	  
	  // A
	  MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
	  MAT_set_j(A,*Acounter,VARGEN_get_index_P(vargen,t)); // Pg
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,Pmis_index[k],-VARGEN_get_P(vargen,t));
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus[k]); bat != NULL; bat = BAT_get_next(bat)) {
	
	//*****************************
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // P var
	  
	  // A
	  MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
	  MAT_set_j(A,*Acounter,BAT_get_index_Pc(bat,t)); // Pc
	  MAT_set_d(A,*Acounter,-1.);
	  (*Acounter)++; 

	  // A
	  MAT_set_i(A,*Acounter,Pmis_index[k]); // Pk
	  MAT_set_j(A,*Acounter,BAT_get_index_Pd(bat,t)); // Pd
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,Pmis_index[k],BAT_get_P(bat,t));
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]*T+t] = TRUE;
  }
}

void CONSTR_DCPF_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {
  // Nothing
}

void CONSTR_DCPF_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
 
  // Local variables
  Bus* bus[2];
  int bus_index[2];
  char* bus_counted;
  int k;
  int num_buses;
  int T;
  int Pmis_index[2];

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Number of buses
  num_buses = NET_get_num_buses(CONSTR_get_network(c));

  // Constr data
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++) {
    bus_index[k] = BUS_get_index(bus[k]);
    Pmis_index[k] = bus_index[k]+t*num_buses;
  }

  // Buses
  for (k = 0; k < 2; k++) {
    
    // Store P balance sensitivity
    if (!bus_counted[bus_index[k]*T+t])
      BUS_set_sens_P_balance(bus[k],VEC_get(sA,Pmis_index[k]),t);
    
    // Update counted flag
    bus_counted[bus_index[k]*T+t] = TRUE;
  }
}

void CONSTR_DCPF_free(Constr* c) {
  // Nothing
}
