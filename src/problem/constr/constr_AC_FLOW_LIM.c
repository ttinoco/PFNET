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


}

void CONSTR_AC_FLOW_LIM_clear(Constr* c) {

  
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
