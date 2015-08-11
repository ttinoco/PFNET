/** @file constr_DCPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type DCPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
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

void CONSTR_DCPF_count_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  int* Acounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  int m;
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !bus_counted)
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
    
    if (!bus_counted[bus_index[k]]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var
	  
	  // A
	  (*Acounter)++; // Pk
	}
      }
      
      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
	  
	  // A
	  (*Acounter)++; // Pk
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
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
  CONSTR_set_b(c,VEC_new(num_buses));

  // A
  CONSTR_set_A(c,MAT_new(num_buses,   // size1 (rows)
			 num_vars,    // size2 (cols)
			 Acounter));  // nnz
}

void CONSTR_DCPF_analyze_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* bus[2];
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Mat* A;
  Vec* rhs;
  int* Acounter;
  char* bus_counted;
  int bus_index[2];

  REAL a;
  REAL a_temp;
  REAL b;
  REAL g;
  REAL g_sh[2];

  REAL sign_phi;

  int k;
  int m;
  
  // Constr data
  A = CONSTR_get_A(c);
  rhs = CONSTR_get_b(c);
  Acounter = CONSTR_get_Acounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !bus_counted)
    return;
 
  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(bus[k]);
 
  // Branch data
  a = BRANCH_get_ratio(br);
  b = BRANCH_get_b(br);
  g = BRANCH_get_g(br);
  g_sh[0] = BRANCH_get_g_from(br);
  g_sh[1] = BRANCH_get_g_to(br);

  // Branch
  //*******

  for (k = 0; k < 2; k++) {

    if (k == 0) {
      m = 1;
      a_temp = a;
      sign_phi = 1;
    }
    else {
      m = 0;
      a_temp = 1;
      sign_phi = -1;
    }

    //***********
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // A 
      MAT_set_i(A,*Acounter,bus_index[k]); // Pk
      MAT_set_j(A,*Acounter,BUS_get_index_v_ang(bus[k])); // wk
      MAT_set_d(A,*Acounter,a*b); 
      (*Acounter)++;
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,bus_index[k],-a*b*BUS_get_v_ang(bus[k]));
    }

    //***********
    if (BUS_has_flags(bus[m],FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // A 
      MAT_set_i(A,*Acounter,bus_index[k]); // Pk
      MAT_set_j(A,*Acounter,BUS_get_index_v_ang(bus[m])); // wk
      MAT_set_d(A,*Acounter,-a*b); 
      (*Acounter)++;
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,bus_index[k],a*b*BUS_get_v_ang(bus[m]));
    }

    //**********
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var

      // A
      MAT_set_i(A,*Acounter,bus_index[k]); // Pk
      MAT_set_j(A,*Acounter,BRANCH_get_index_phase(br)); // phi
      MAT_set_d(A,*Acounter,-a*b*sign_phi);
      (*Acounter)++; 
    }
    else {
      
      // b 
      VEC_add_to_entry(rhs,bus_index[k],a*b*BRANCH_get_phase(br)*sign_phi);
    }
    
    // b 
    VEC_add_to_entry(rhs,bus_index[k],a_temp*a_temp*(g+g_sh[k])-a*g);
  }
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	//*****************************
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var
	  
	  // A
	  MAT_set_i(A,*Acounter,bus_index[k]); // Pk
	  MAT_set_j(A,*Acounter,GEN_get_index_P(gen)); // Pg
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,bus_index[k],-GEN_get_P(gen));
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus[k]); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	//*****************************
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
	  
	  // A
	  MAT_set_i(A,*Acounter,bus_index[k]); // Pk
	  MAT_set_j(A,*Acounter,VARGEN_get_index_P(vargen)); // Pg
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; 
	}
	else {
	  
	  // b
	  VEC_add_to_entry(rhs,bus_index[k],-VARGEN_get_P(vargen));
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// b
	VEC_add_to_entry(rhs,bus_index[k],SHUNT_get_g(shunt));
      }

      // Loads
      for (load = BUS_get_load(bus[k]); load != NULL; load = LOAD_get_next(load)) {
	
	// b
	VEC_add_to_entry(rhs,bus_index[k],LOAD_get_P(load));
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_DCPF_eval_branch(Constr* c, Branch *br, Vec* var_values) {
  // Nothing
}

void CONSTR_DCPF_store_sens_branch(Constr* c, Branch* br, Vec* sens) {
  // Nothing
}

void CONSTR_DCPF_free(Constr* c) {
  // Nothing
}
