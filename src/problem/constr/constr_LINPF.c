/** @file constr_LINPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LINPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LINPF.h>
#include <pfnet/constr_ACPF.h>

Constr* CONSTR_LINPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_LINPF_init);
  CONSTR_set_func_count_step(c,&CONSTR_LINPF_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_LINPF_allocate);
  CONSTR_set_func_clear(c,&CONSTR_LINPF_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_LINPF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_LINPF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_LINPF_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_LINPF_free);
  CONSTR_set_name(c,"linearized AC power balance");
  CONSTR_init(c);
  return c;
}

void CONSTR_LINPF_init(Constr* c) {

  // Init
  Constr* acpf = CONSTR_ACPF_new(CONSTR_get_network(c));
  CONSTR_set_data(c,(void*)acpf);
}

void CONSTR_LINPF_clear(Constr* c) {
  
  // Local vars
  Constr* acpf = (Constr*)CONSTR_get_data(c);

  // ACPF clear
  CONSTR_clear(acpf);
}

void CONSTR_LINPF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local vars
  Constr* acpf = (Constr*)CONSTR_get_data(c);
  Net* net = CONSTR_get_network(c);

  // ACPF count
  CONSTR_count_step(acpf,bus,busdc,t);

  // Done 
  if ((t == BUS_get_num_periods(bus)-1) && (BUS_get_index(bus) == NET_get_num_buses(net,FALSE)-1)) {
    CONSTR_set_A_row(c,CONSTR_get_J_row(acpf));
    CONSTR_set_A_nnz(c,CONSTR_get_J_nnz(acpf));
  }
}

void CONSTR_LINPF_allocate(Constr* c) {
  
  // Local vars
  Constr* acpf = (Constr*)CONSTR_get_data(c);
  
  // ACPF
  CONSTR_allocate(acpf); 
}

void CONSTR_LINPF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local vars
  Constr* acpf;
  Net* net;
  Vec* f;
  Mat* J;
  Mat* A;
  Vec* x0;
  Vec* b;
  int i;

  // Data
  net = CONSTR_get_network(c);
  A = CONSTR_get_A(c);

  // ACPF
  acpf = (Constr*)CONSTR_get_data(c);
  CONSTR_analyze_step(acpf,bus,busdc,t);

  // Done 
  if ((t == BUS_get_num_periods(bus)-1) && (BUS_get_index(bus) == NET_get_num_buses(net,FALSE)-1)) {
    x0 = NET_get_var_values(net,CURRENT);
    CONSTR_eval(acpf,x0,NULL);
    J = CONSTR_get_J(acpf);
    f = CONSTR_get_f(acpf);
    b = MAT_rmul_by_vec(J,x0);
    VEC_sub_inplace(b,f);
    CONSTR_set_b(c,b);
    for (i = 0; i < MAT_get_nnz(J); i++) {
      MAT_set_i(A,i,MAT_get_i(J,i));
      MAT_set_j(A,i,MAT_get_j(J,i));
      MAT_set_d(A,i,MAT_get_d(J,i));
    }
  }
}

void CONSTR_LINPF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_LINPF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing for now
}
 
void CONSTR_LINPF_free(Constr* c) {
  
  // ACPF
  Constr* acpf = (Constr*)CONSTR_get_data(c);
  CONSTR_del(acpf);
  CONSTR_set_data(c,NULL);
}
