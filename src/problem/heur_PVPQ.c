/** @file heur_PVPQ.c
 *  @brief This file defines the data structure and routines associated with the heuristic of type PVPQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/heur_PVPQ.h>

struct Heur_PVPQ_Data {
  
  char* reg_flag; // flags for tracking regulation
};

void HEUR_PVPQ_init(Heur* h, Net* net) {

  // Local variables
  int num_buses;
  Bus* bus;
  Heur_PVPQ_Data* data;
  int i;

  // Init
  num_buses = NET_get_num_buses(net);
  HEUR_set_bus_counted(h,(char*)calloc(num_buses,sizeof(char)));
  data = (Heur_PVPQ_Data*)malloc(sizeof(Heur_PVPQ_Data));
  data->reg_flag = (char*)malloc(sizeof(char)*num_buses);
  for (i = 0; i < num_buses; i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_regulated_by_gen(bus))
      data->reg_flag[i] = TRUE;
    else
      data->reg_flag[i] = FALSE;    
  }
  HEUR_set_data(h,(void*)data);
}

void HEUR_PVPQ_clear(Heur* h, Net* net) {
  
  // Local variables
  int num_buses;

  // Clear bus counted flags
  num_buses = NET_get_num_buses(net);
  HEUR_clear_bus_counted(h,num_buses);
  
}

void HEUR_PVPQ_apply_to_branch(Heur* h, Constr* clist, Net* net, Branch* br, Vec* var_values) {

  // Local variables
  Vec* f;
  Mat* A;
  Vec* b;
  Bus* bus[2];
  Gen* gen;
  char* bus_counted;
  Heur_PVPQ_Data* data;
  char* reg_flag;
  int bus_index[2];
  int k;
  int i;
  Constr* pf;
  Constr* fix;
  REAL v;
  REAL v_set;
  REAL Q;
  REAL Qmax;
  REAL Qmin;
  char switch_flag;
  int j_old;
  int j_new;
  REAL b_new;
  
  // Heur data
  bus_counted = HEUR_get_bus_counted(h);
  data = (Heur_PVPQ_Data*)HEUR_get_data(h);
  reg_flag = data->reg_flag;

  // Bus from data
  bus[0] = BRANCH_get_bus_from(br);
  bus_index[0] = BUS_get_index(bus[0]);
  
  // Bus to data
  bus[1] = BRANCH_get_bus_to(br);
  bus_index[1] = BUS_get_index(bus[1]);

  // Power flow constraints
  for (pf = clist; pf != NULL; pf = CONSTR_get_next(pf)) {
    if (CONSTR_get_type(pf) == CONSTR_TYPE_PF)
      break;
  }
  if (!pf)
    return;

  // Fix constraints
  for (fix = clist; fix != NULL; fix = CONSTR_get_next(fix)) {
    if (CONSTR_get_type(fix) == CONSTR_TYPE_FIX)
      break;
  }
  if (!fix)
    return;

  // Constr data
  f = CONSTR_get_f(pf);
  A = CONSTR_get_A(fix);
  b = CONSTR_get_b(fix);

  // Buses
  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]] &&                     // not counted
	!BUS_is_slack(bus[k]) &&                          // not slack
	BUS_is_regulated_by_gen(bus[k]) &&                // regulated
	BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG) &&   // v mag is variable
	BUS_has_flags(bus[k],FLAG_FIXED,BUS_VAR_VMAG) &&  // v mag is fixed
	GEN_has_flags(BUS_get_reg_gen(bus[k]),FLAG_VARS,GEN_VAR_Q)) { // reg gen Q is variable

      // Voltage magnitude
      v = VEC_get(var_values,BUS_get_index_v_mag(bus[k]));
      v_set = BUS_get_v_set(bus[k]);

      // Regulating generator (first one in list of reg gens)
      gen = BUS_get_reg_gen(bus[k]);
      Q = VEC_get(var_values,GEN_get_index_Q(gen)); // per unit
      Qmax = GEN_get_Q_max(gen);                    // per unit
      Qmin = GEN_get_Q_min(gen);                    // per unit
      
      // Switch flag
      switch_flag = FALSE;
		 
      // Currently regulated
      if (reg_flag[bus_index[k]]) {
	
	// Violations
	if (Q > Qmax) {

	  // Set data
	  j_old = BUS_get_index_v_mag(bus[k]);
	  j_new = GEN_get_index_Q(gen);
	  b_new = Qmax;
	  switch_flag = TRUE;
	  reg_flag[bus_index[k]] = FALSE;
	  
	  // Update vector of var values
	  while (gen) {
	    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	      VEC_set(var_values,GEN_get_index_Q(gen),GEN_get_Q_max(gen));
	    gen = GEN_get_reg_next(gen);
	  }
	  
	  // Show
          #ifdef DEBUG
	    printf("PV-PQ switching - Q > Qmax (%d)\n",BUS_get_number(bus[k]));
	  #endif
	}
	else if (Q < Qmin) {

	  // Set data
	  j_old = BUS_get_index_v_mag(bus[k]);
	  j_new = GEN_get_index_Q(gen);
	  b_new = Qmin;
	  switch_flag = TRUE;
	  reg_flag[bus_index[k]] = FALSE;

	  // Update vector of var values
	  while (gen) {
	    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	      VEC_set(var_values,GEN_get_index_Q(gen),GEN_get_Q_min(gen));
	    gen = GEN_get_reg_next(gen);
	  }

	  // Show
	  #ifdef DEBUG
	    printf("PV-PQ switching - Q < Qmin (%d)\n",BUS_get_number(bus[k]));
	  #endif
	}
      }

      // Previously regulated
      else {
	
	// Q at Qmin and v < v_set
	if (fabs(Q-Qmin) < fabs(Q-Qmax) && v < v_set) {
	    
	  Q = Q - VEC_get(f,2*BUS_get_index(GEN_get_bus(gen))+1); // per unit (see constr_PF)

	  if (Q >= Qmax) {

	    // Set data
	    j_old = GEN_get_index_Q(gen);
	    j_new = GEN_get_index_Q(gen);
	    b_new = Qmax;
	    switch_flag = TRUE;

	    // Update vector of var values
	    while (gen) {
	      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
		VEC_set(var_values,GEN_get_index_Q(gen),GEN_get_Q_max(gen));
	      gen = GEN_get_reg_next(gen);
	    }

	    // Show
	    #ifdef DEBUG
	      printf("PQ-PQ switching - Qmin to Qmax (%d)\n",BUS_get_number(bus[k]));
	    #endif
	  }
	  else if (Qmin < Q && Q < Qmax) {

	    // Set data
	    j_old = GEN_get_index_Q(gen);
	    j_new = BUS_get_index_v_mag(bus[k]);
	    b_new = v_set;
	    switch_flag = TRUE;	      
	    reg_flag[bus_index[k]] = TRUE;

	    // Udpate vector of var values
	    VEC_set(var_values,j_new,b_new);

	    // Show
	    #ifdef DEBUG
	      printf("PQ-PV switching - Qmin to vset (%d)\n",BUS_get_number(bus[k]));
	    #endif
	  }
	}

	// Q at Qmax and v > v_set
	else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) {
	  
	  Q = Q - VEC_get(f,2*BUS_get_index(GEN_get_bus(gen))+1); // per unit (see constr_PF)

	  if (Q <= Qmin) {

	    // Set data
	    j_old = GEN_get_index_Q(gen);
	    j_new = GEN_get_index_Q(gen);
	    b_new = Qmin;
	    switch_flag = TRUE;

	    // Update vector of var values
	    while (gen) {
	      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
		VEC_set(var_values,GEN_get_index_Q(gen),GEN_get_Q_min(gen));
	      gen = GEN_get_reg_next(gen);
	    }

	    // Show
	    #ifdef DEBUG
	      printf("PQ-PQ switching - Qmax to Qmin (%d)\n",BUS_get_number(bus[k]));
	    #endif
	  }
	  else if (Qmin < Q && Q < Qmax) {

	    // Set data
	    j_old = GEN_get_index_Q(gen);
	    j_new = BUS_get_index_v_mag(bus[k]);
	    b_new = v_set;
	    switch_flag = TRUE;	      
	    reg_flag[bus_index[k]] = TRUE;

	    // Udpate vector of var values
	    VEC_set(var_values,j_new,b_new);

	    // Show 
	    #ifdef DEBUG
	      printf("PQ-PV switching - Qmax to vset (%d)\n",BUS_get_number(bus[k]));
	    #endif
	  }
	}
      }

      // Update fix constraints
      if (switch_flag) {
	for (i = 0; i < MAT_get_nnz(A); i++) {
	  if (MAT_get_j(A,i) == j_old)
	    MAT_set_d(A,i,0.);
	  if (MAT_get_j(A,i) == j_new) {
	    MAT_set_d(A,i,1.);
	    VEC_set(b,MAT_get_i(A,i),b_new);
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE; 
  }
}

void HEUR_PVPQ_free(Heur* h) {

  // Local variables
  Heur_PVPQ_Data* data;

  // Get data
  data = (Heur_PVPQ_Data*)HEUR_get_data(h);

  // Free
  if (data) {
    free(data->reg_flag);
  }
  free(data);

  // Set data
  HEUR_set_data(h,NULL);
}
