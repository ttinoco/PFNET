/** @file constr_GEN_RAMP.c
 *  @brief This file defines the data structure and routines associated with the constraint of type GEN_RAMP.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_GEN_RAMP.h>

Constr* CONSTR_GEN_RAMP_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_GEN_RAMP_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_GEN_RAMP_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_GEN_RAMP_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_GEN_RAMP_store_sens_step);
  CONSTR_set_name(c,"generator ramp limits");
  return c;
}

void CONSTR_GEN_RAMP_count_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Gen* gen;
  int* G_nnz;
  int* G_row;

  // Constr data
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);

  // Check pointer
  if (!G_nnz || !G_row)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    // Variable
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // -dP_max <= P_t - P_{t-1} <= dP_max
      if (t == 0)
	(*G_nnz) += 1;
      else
	(*G_nnz) += 2;
      (*G_row)++;
    }
  }
}

void CONSTR_GEN_RAMP_analyze_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Gen* gen;
  int* G_nnz;
  int* G_row;
  Vec* u;
  Vec* l;
  Mat* G;
  
  // Cosntr data
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G = CONSTR_get_G(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);

  // Check pointers
  if (!G_nnz || !G_row)
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    // Variables
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // -dP_max <= P_t - P_{t-1} <= dP_max
      
      // G
      MAT_set_i(G,*G_nnz,*G_row);
      MAT_set_j(G,*G_nnz,GEN_get_index_P(gen,t));
      MAT_set_d(G,*G_nnz,1.);
      
      if (t == 0) {
	
	// l u
	VEC_set(l,*G_row,-GEN_get_dP_max(gen)+GEN_get_P_prev(gen));
	VEC_set(u,*G_row,GEN_get_dP_max(gen)+GEN_get_P_prev(gen));
	
	(*G_nnz) += 1;
      }
      else {
	
	// l u
	VEC_set(l,*G_row,-GEN_get_dP_max(gen));
	VEC_set(u,*G_row,GEN_get_dP_max(gen));
	
	// G
	MAT_set_i(G,*G_nnz+1,*G_row);
	MAT_set_j(G,*G_nnz+1,GEN_get_index_P(gen,t-1));
	MAT_set_d(G,*G_nnz+1,-1.);
	
	(*G_nnz) += 2;
      }
      
      (*G_row)++;
    }
  }
}

void CONSTR_GEN_RAMP_eval_step(Constr* c, Bus* bus, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_GEN_RAMP_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing for now
}
