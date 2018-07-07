/** @file constr_DC_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type DC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_DC_FLOW_LIM.h>

Constr* CONSTR_DC_FLOW_LIM_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_DC_FLOW_LIM_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_DC_FLOW_LIM_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_DC_FLOW_LIM_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_DC_FLOW_LIM_store_sens_step);
  CONSTR_set_name(c,"DC branch flow limits");
  return c;
}

void CONSTR_DC_FLOW_LIM_count_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Branch* br;
  int* G_nnz;
  int* G_row;
  Bus* buses[2];

  // Constr data
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);

  // Check pointer
  if (!G_nnz || !G_row)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      continue;
    
    // Zero limits
    if (BRANCH_get_ratingA(br) == 0.)
      continue;
    
    buses[0] = BRANCH_get_bus_k(br);
    buses[1] = BRANCH_get_bus_m(br);
    
    if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // G
      (*G_nnz)++;
    }
    
    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // G
      (*G_nnz)++;
    }
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // G
      (*G_nnz)++;
    }
    
    // Constraint index
    (*G_row)++;
  }
}

void CONSTR_DC_FLOW_LIM_analyze_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Bus* buses[2];
  Mat* G;
  Vec* l;
  Vec* u;
  int* G_nnz;
  int* G_row;
  REAL b;
  double rating;

  // Constr data
  G = CONSTR_get_G(c);
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);

  // Check pointer
  if (!G_nnz || !G_row)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      continue;
    
    // Zero limits
    if (BRANCH_get_ratingA(br) == 0.)
      continue;
    
    buses[0] = BRANCH_get_bus_k(br);
    buses[1] = BRANCH_get_bus_m(br);
    
    b = BRANCH_get_b(br);
    
    rating = BRANCH_get_ratingA(br); // p.u.
    
    VEC_set(l,*G_row,-rating); // p.u.
    VEC_set(u,*G_row,rating);  // p.u.
    
    if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
      
      // G
      MAT_set_i(G,*G_nnz,*G_row);
      MAT_set_j(G,*G_nnz,BUS_get_index_v_ang(buses[0],t)); // wk
      MAT_set_d(G,*G_nnz,-b);
      (*G_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(l,*G_row,b*BUS_get_v_ang(buses[0],t));
      VEC_add_to_entry(u,*G_row,b*BUS_get_v_ang(buses[0],t));
    }
    
    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG)) { // wm var
      
      // G
      MAT_set_i(G,*G_nnz,*G_row);
      MAT_set_j(G,*G_nnz,BUS_get_index_v_ang(buses[1],t)); // wk
      MAT_set_d(G,*G_nnz,b);
      (*G_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(l,*G_row,-b*BUS_get_v_ang(buses[1],t));
      VEC_add_to_entry(u,*G_row,-b*BUS_get_v_ang(buses[1],t));
    }
    
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) { // phi var
      
      // G
      MAT_set_i(G,*G_nnz,*G_row);
      MAT_set_j(G,*G_nnz,BRANCH_get_index_phase(br,t)); // phi
      MAT_set_d(G,*G_nnz,b);
      (*G_nnz)++;
    }
    else {
      
      // b
      VEC_add_to_entry(l,*G_row,-b*BRANCH_get_phase(br,t));
      VEC_add_to_entry(u,*G_row,-b*BRANCH_get_phase(br,t));
    }
    
    // Constraint index
    (*G_row)++;
  }
}

void CONSTR_DC_FLOW_LIM_eval_step(Constr* c, Bus* bus, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_DC_FLOW_LIM_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Branch* br;
  int* G_row;

  // Constr data
  G_row = CONSTR_get_G_row_ptr(c);

  // Check pointer
  if (!G_row)
    return;

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Check outage
    if (BRANCH_is_on_outage(br))
      continue;
    
    // Zero limits
    if (BRANCH_get_ratingA(br) == 0.)
      continue;
    
    // Store sensitivies
    BRANCH_set_sens_P_u_bound(br,VEC_get(sGu,*G_row),t);
    BRANCH_set_sens_P_l_bound(br,VEC_get(sGl,*G_row),t);
    
    // Constraint index
    (*G_row)++;
  }
}
