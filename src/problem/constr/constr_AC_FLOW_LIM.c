/** @file constr_AC_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_AC_FLOW_LIM.h>

void CONSTR_AC_FLOW_LIM_init(Constr* c) {
  
  // Local variables
  Net* net;
  int num_branches;
  int num_periods;

  // Init
  net = CONSTR_get_network(c);
  num_branches = NET_get_num_branches(net);
  num_periods = NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_branches*num_periods,sizeof(int)),num_branches*num_periods);
  CONSTR_set_data(c,NULL);

}

void CONSTR_AC_FLOW_LIM_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_J_nnz(c,0);
  CONSTR_set_Jbar_nnz(c,0);
  CONSTR_set_Jconstr_index(c,0);
  CONSTR_clear_H_nnz(c);

  // Flags
  CONSTR_clear_bus_counted(c);
  
}

void CONSTR_AC_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {

}

void CONSTR_AC_FLOW_LIM_allocate(Constr* c) {


}

void CONSTR_AC_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t) {


}

void CONSTR_AC_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {
  // Nothing
}

void CONSTR_AC_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_AC_FLOW_LIM_free(Constr* c) {
  // Nothing
}
