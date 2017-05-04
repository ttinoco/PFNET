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

Constr* CONSTR_AC_LIN_FLOW_LIM_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_AC_LIN_FLOW_LIM_init);
  CONSTR_set_func_count_step(c, &CONSTR_AC_LIN_FLOW_LIM_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_AC_LIN_FLOW_LIM_allocate);
  CONSTR_set_func_clear(c, &CONSTR_AC_LIN_FLOW_LIM_clear);
  CONSTR_set_func_analyze_step(c, &CONSTR_AC_LIN_FLOW_LIM_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_AC_LIN_FLOW_LIM_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_AC_LIN_FLOW_LIM_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_AC_LIN_FLOW_LIM_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_AC_LIN_FLOW_LIM_init(Constr* c) {
    
  // Init
  CONSTR_set_name(c,"linearized AC branch flow limits");
  CONSTR_set_data(c,NULL);
}

void CONSTR_AC_LIN_FLOW_LIM_clear(Constr* c) {

}

void CONSTR_AC_LIN_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {


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
