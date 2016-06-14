/** @file constr_LINPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LINPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LINPF.h>

void CONSTR_LINPF_init(Constr* c) {

  // Init
  Constr* acpf = CONSTR_new(CONSTR_TYPE_PF,CONSTR_get_network(c));
  CONSTR_set_data(c,(void*)acpf);
}

void CONSTR_LINPF_clear(Constr* c) {
  
  // ACPF
  Constr* acpf = (Constr*)CONSTR_get_data(c);
  CONSTR_clear(acpf);
}

void CONSTR_LINPF_count_branch(Constr* c, Branch* br) {

  // ACPF
  Constr* acpf = (Constr*)CONSTR_get_data(c);
  CONSTR_count_branch(acpf,br);
}

void CONSTR_LINPF_allocate(Constr* c) {
  
  // Local variables
  Net* net;
  int num_vars;
  Constr* acpf;
  Vec* f;
  Mat* J;
  Vec* x0;

  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);
  acpf = (Constr*)CONSTR_get_data(c);

  // PF
 
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

void CONSTR_LINPF_analyze_branch(Constr* c, Branch* br) {
  
  
}

void CONSTR_LINPF_eval_branch(Constr* c, Branch* br, Vec* var_values) {
  // Nothing
}

void CONSTR_LINPF_store_sens_branch(Constr* c, Branch* br, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
}

void CONSTR_LINPF_free(Constr* c) {
  // Nothing
}
