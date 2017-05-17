/** @file constr_AC_LIN_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_LIN_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_AC_LIN_FLOW_LIM.h>

struct Constr_AC_LIN_FLOW_LIM_Data {
  #if HAVE_LINE_FLOW
  LINE_FLOW_Results** results; // array of pointers to LINE_FLOW results (2 per branch and time)
  int size;
  #endif  
};

Constr* CONSTR_AC_LIN_FLOW_LIM_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_AC_LIN_FLOW_LIM_init);
  CONSTR_set_func_count_step(c,&CONSTR_AC_LIN_FLOW_LIM_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_AC_LIN_FLOW_LIM_allocate);
  CONSTR_set_func_clear(c,&CONSTR_AC_LIN_FLOW_LIM_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_AC_LIN_FLOW_LIM_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_AC_LIN_FLOW_LIM_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_AC_LIN_FLOW_LIM_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_AC_LIN_FLOW_LIM_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_AC_LIN_FLOW_LIM_init(Constr* c) {
    
  // Local variables
  int i;
  int num;
  Net* net;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Init
  net = CONSTR_get_network(c);
  num = 2*NET_get_num_branches(net)*NET_get_num_periods(net); // max number of LINE_FLOW results
  data = (Constr_AC_LIN_FLOW_LIM_Data*)malloc(sizeof(Constr_AC_LIN_FLOW_LIM_Data));
  #if HAVE_LINE_FLOW
  data->results = (LINE_FLOW_Results**)malloc(sizeof(LINE_FLOW_Results*)*num);
  data->size = num;
  for (i = 0; i < num; i++)
    data->results[i] = NULL;
  #endif
  CONSTR_set_name(c,"linearized AC branch flow limits");
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_AC_LIN_FLOW_LIM_clear(Constr* c) {
  
  // Counters
  CONSTR_set_G_nnz(c,0); // number of nonzeros of G matrix (from l <= Gx <= u)
  CONSTR_set_G_row(c,0); // number of rows of G matrix (from l <= Gx <= u)
}

void CONSTR_AC_LIN_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {

  #if HAVE_LINE_FLOW

  // Local variables
  int k;
  int* G_nnz;
  int* G_row;
  Net* net;
  Bus* bus[2];
  REAL V_min[2];
  REAL V_max[2];
  int offset;
  int num_constr;
  LINE_FLOW_Results* results1;
  LINE_FLOW_Results* results2;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Constr data
  net = CONSTR_get_network(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  data = (Constr_AC_LIN_FLOW_LIM_Data*)CONSTR_get_data(c);

  // Offset
  offset = 2*NET_get_num_branches(net)*t;
  
  // Check pointer
  if (!G_nnz || !G_row || !data) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint has undefined pointers");  
    return;
  }

  // Outage (skip)
  if (BRANCH_is_on_outage(br))
    return;

  // Zero limits (skip)
  if (BRANCH_get_ratingA(br) == 0.)
    return;
  
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);

  // Check tap ratio 
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) ||
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint does not support tap ratios or phase shifts as variables");
    return;
  }

  // Check angles
  if (!BUS_has_flags(bus[0],FLAG_VARS,BUS_VAR_VANG) &&
      !BUS_has_flags(bus[1],FLAG_VARS,BUS_VAR_VANG)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires at least one voltage angle as variable accross branch");
    return;
  }
  
  // Volage limits
  for (k = 0; k < 2; k++) {
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) {
      if (!BUS_has_flags(bus[k],FLAG_BOUNDED,BUS_VAR_VMAG)) {
	CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires variable voltage magnitudes to be bounded");
	return;
      }
      V_min[k] = BUS_get_v_min_norm(bus[k]);
      V_max[k] = BUS_get_v_max_norm(bus[k]);
    }
    else {
      V_min[k] = BUS_get_v_mag(bus[k],t);
      V_max[k] = BUS_get_v_mag(bus[k],t);
    }
  }

  // Construct constraints Ikm <= Imax
  results1 = LINE_FLOW_construct(V_min[0], 
                                 V_max[0], 
                                 V_min[1], 
                                 V_max[1],
                                 BRANCH_get_g(br), 
                                 BRANCH_get_b(br),
                                 BRANCH_get_b_k(br),
                                 1./BRANCH_get_ratio(br,t),
                                 BRANCH_get_phase(br,t)*180./PI, // degrees
                                 BRANCH_get_ratingA(br),
                                 10, 
                                 5., 
                                 1,     // "k" side
                                 NULL);

  // Construct constraints Imk <= Imax
  results2 = LINE_FLOW_construct(V_min[0], 
                                 V_max[0], 
                                 V_min[1], 
                                 V_max[1],
                                 BRANCH_get_g(br), 
                                 BRANCH_get_b(br),
                                 BRANCH_get_b_k(br), // What about b shunt on the "m" side?
                                 1./BRANCH_get_ratio(br,t),
                                 BRANCH_get_phase(br,t)*180./PI, // degrees
                                 BRANCH_get_ratingA(br),
                                 10, 
                                 5.,    
                                 2,     // "m" side
                                 NULL);

  // Save results
  data->results[2*BRANCH_get_index(br)+0+offset] = results1; 
  data->results[2*BRANCH_get_index(br)+1+offset] = results2;
 
  // Check results
  if (LINE_FLOW_get_flag(results1) != 2) {
    CONSTR_set_error(c,LINE_FLOW_get_message(results1));    
    return;
  }
  if (LINE_FLOW_get_flag(results2) != 2) {
    CONSTR_set_error(c,LINE_FLOW_get_message(results2));    
    return;
  }
 
  // Count
  num_constr = LINE_FLOW_get_number_constraints(results1)+LINE_FLOW_get_number_constraints(results2);
  (*G_row) += num_constr;
  for (k = 0; k < 2; k++) {
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG))
      (*G_nnz) += num_constr;
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG))
      (*G_nnz) += num_constr;
  }

  #else
  CONSTR_set_error(c,"unable to construct AC_LIN_FLOW_LIM constraint - line_flow library not available");
  #endif  
}

void CONSTR_AC_LIN_FLOW_LIM_allocate(Constr* c) {
  

}

void CONSTR_AC_LIN_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t) {


}

void CONSTR_AC_LIN_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  

}

void CONSTR_AC_LIN_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing yet
}

void CONSTR_AC_LIN_FLOW_LIM_free(Constr* c) {
  // Nothing
}
