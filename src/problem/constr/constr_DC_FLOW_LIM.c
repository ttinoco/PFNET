/** @file constr_DC_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type DC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_DC_FLOW_LIM.h>

void CONSTR_DC_FLOW_LIM_init(Constr* c) {
  
  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_DC_FLOW_LIM_clear(Constr* c) {
    
  // Counters
  CONSTR_set_Gcounter(c,0);
  CONSTR_set_Gconstr_index(c,0);
}

void CONSTR_DC_FLOW_LIM_count_branch(Constr* c, Branch* br) {
  
  // Local variables
  int* Gcounter;
  int* Gconstr_index;
  Bus* bus[2];
  
  // Constr data
  Gcounter = CONSTR_get_Gcounter_ptr(c);
  Gconstr_index = CONSTR_get_Gconstr_index_ptr(c);
  
  // Check pointer
  if (!Gcounter || !Gconstr_index)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  
  if (BUS_has_flags(bus[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
    
    // G
    (*Gcounter)++;
  }
  
  if (BUS_has_flags(bus[1],FLAG_VARS,BUS_VAR_VANG)) { // wm var
    
    // G
    (*Gcounter)++;
  }
  
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
    
    // G
    (*Gcounter)++;
  } 

  // Constraint index
  (*Gconstr_index)++;
}

void CONSTR_DC_FLOW_LIM_allocate(Constr* c) {
  
  // Local variables
  Net* net;
  int num_vars;
  int Gcounter;
  int Gconstr_index;

  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);
  Gcounter = CONSTR_get_Gcounter(c);
  Gconstr_index = CONSTR_get_Gconstr_index(c);
  
  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));
  
  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));
  
  // h
  CONSTR_set_l(c,VEC_new(Gconstr_index));
  CONSTR_set_u(c,VEC_new(Gconstr_index));

  // G
  CONSTR_set_G(c,MAT_new(Gconstr_index, // size1 (rows)
			 num_vars,      // size2 (cols)
			 Gcounter));    // nnz
}

void CONSTR_DC_FLOW_LIM_analyze_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* bus[2];
  Mat* G;
  Vec* l;
  Vec* u;
  int* Gcounter;
  int* Gconstr_index;
  REAL b;
  double rating;
  
  // Constr data
  G = CONSTR_get_G(c);
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  Gcounter = CONSTR_get_Gcounter_ptr(c);
  Gconstr_index = CONSTR_get_Gconstr_index_ptr(c);

  // Check pointer
  if (!Gcounter || !Gconstr_index)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  
  b = BRANCH_get_b(br);
  
  if (BRANCH_get_ratingA(br) > 0)
    rating = BRANCH_get_ratingA(br);
  else
    rating = BRANCH_INF_FLOW;

  VEC_set(l,*Gconstr_index,-rating); // p.u.
  VEC_set(u,*Gconstr_index,rating);  // p.u.
  
  if (BUS_has_flags(bus[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
    
    // G
    MAT_set_i(G,*Gcounter,*Gconstr_index);
    MAT_set_j(G,*Gcounter,BUS_get_index_v_ang(bus[0])); // wk
    MAT_set_d(G,*Gcounter,-b);
    (*Gcounter)++;
  }
  else {
    
    // b 
    VEC_add_to_entry(l,*Gconstr_index,b*BUS_get_v_ang(bus[0]));
    VEC_add_to_entry(u,*Gconstr_index,b*BUS_get_v_ang(bus[0]));
  }

  if (BUS_has_flags(bus[1],FLAG_VARS,BUS_VAR_VANG)) { // wm var
    
    // G
    MAT_set_i(G,*Gcounter,*Gconstr_index);
    MAT_set_j(G,*Gcounter,BUS_get_index_v_ang(bus[1])); // wk
    MAT_set_d(G,*Gcounter,b);
    (*Gcounter)++;
  }
  else {
    
    // b 
    VEC_add_to_entry(l,*Gconstr_index,-b*BUS_get_v_ang(bus[1]));
    VEC_add_to_entry(u,*Gconstr_index,-b*BUS_get_v_ang(bus[1]));
  }

  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
    
    // G
    MAT_set_i(G,*Gcounter,*Gconstr_index);
    MAT_set_j(G,*Gcounter,BRANCH_get_index_phase(br)); // phi
    MAT_set_d(G,*Gcounter,b);
    (*Gcounter)++;
  }
  else {
    
    // b 
    VEC_add_to_entry(l,*Gconstr_index,-b*BRANCH_get_phase(br));
    VEC_add_to_entry(u,*Gconstr_index,-b*BRANCH_get_phase(br));
  }

  // Constraint index
  (*Gconstr_index)++;
}

void CONSTR_DC_FLOW_LIM_eval_branch(Constr* c, Branch* br, Vec* var_values) {
  // Nothing
}

void CONSTR_DC_FLOW_LIM_store_sens_branch(Constr* c, Branch* br, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  int* Gconstr_index;

  // Constr data
  Gconstr_index = CONSTR_get_Gconstr_index_ptr(c);

  // Check pointer
  if (!Gconstr_index)
    return;
  
  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  // Store sensitivies
  BRANCH_set_sens_P_u_bound(br,VEC_get(sGu,*Gconstr_index));
  BRANCH_set_sens_P_l_bound(br,VEC_get(sGl,*Gconstr_index));

  // Constraint index
  (*Gconstr_index)++;
}

void CONSTR_DC_FLOW_LIM_free(Constr* c) {
  // Nothing
}
