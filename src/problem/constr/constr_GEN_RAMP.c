/** @file constr_GEN_RAMP.c
 *  @brief This file defines the data structure and routines associated with the constraint of type GEN_RAMP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_GEN_RAMP.h>
#include <assert.h>

void CONSTR_GEN_RAMP_init(Constr* c) {
  
  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_GEN_RAMP_clear(Constr* c) {
  
  // Counters
  CONSTR_set_Gcounter(c,0);
  CONSTR_set_Gconstr_index(c,0);
  
  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_GEN_RAMP_count_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int* Gcounter;
  int* Gconstr_index;
  char* bus_counted;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);
  
  // Constr data
  Gcounter = CONSTR_get_Gcounter_ptr(c);
  Gconstr_index = CONSTR_get_Gconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!Gcounter || !Gconstr_index || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Variable
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // -dP_max <= P_t - P_{t-1} <= dP_max
	  if (t == 0)
	    (*Gcounter) += 1;
	  else
	    (*Gcounter) += 2;
	  (*Gconstr_index)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;    
  }
}

void CONSTR_GEN_RAMP_allocate(Constr* c) {
  
  // Local variables
  int num_constr;
  int num_vars;
  int Gcounter;
  
  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_constr = CONSTR_get_Gconstr_index(c);
  Gcounter = CONSTR_get_Gcounter(c);

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // G l u
  CONSTR_set_l(c,VEC_new(num_constr));
  CONSTR_set_u(c,VEC_new(num_constr));
  CONSTR_set_G(c,MAT_new(num_constr, // size1 (rows)
			 num_vars,   // size2 (rows)
			 Gcounter)); // nnz
}

void CONSTR_GEN_RAMP_analyze_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  int* Gcounter;
  int* Gconstr_index;
  char* bus_counted;
  Vec* u;
  Vec* l;
  Mat* G;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);
  
  // Cosntr data
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G = CONSTR_get_G(c);
  Gcounter = CONSTR_get_Gcounter_ptr(c);
  Gconstr_index = CONSTR_get_Gconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Gcounter || !Gconstr_index || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)*T+t]) {
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Variables
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // -dP_max <= P_t - P_{t-1} <= dP_max
	  
	  // G
	  MAT_set_i(G,*Gcounter,*Gconstr_index);
	  MAT_set_j(G,*Gcounter,GEN_get_index_P(gen,t));
	  MAT_set_d(G,*Gcounter,1.);
	  
	  if (t == 0) {

	    // l u
	    VEC_set(l,*Gconstr_index,-GEN_get_dP_max(gen)+GEN_get_P_prev(gen));
	    VEC_set(u,*Gconstr_index,GEN_get_dP_max(gen)+GEN_get_P_prev(gen));
	    
	    (*Gcounter) += 1;
	  }
	  else {
	    
	    // l u
	    VEC_set(l,*Gconstr_index,-GEN_get_dP_max(gen));
	    VEC_set(u,*Gconstr_index,GEN_get_dP_max(gen));
	    
	    // G
	    MAT_set_i(G,*Gcounter+1,*Gconstr_index);
	    MAT_set_j(G,*Gcounter+1,GEN_get_index_P(gen,t-1));
	    MAT_set_d(G,*Gcounter+1,-1.);
	    
	    (*Gcounter) += 2;
	  }
	  
	  (*Gconstr_index)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;   
  }  
}

void CONSTR_GEN_RAMP_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {
  // Nothing to do
}

void CONSTR_GEN_RAMP_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing for now
}

void CONSTR_GEN_RAMP_free(Constr* c) {
  // Nothing to do
}
