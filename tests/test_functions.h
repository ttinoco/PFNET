/** @file test_constraints.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/pfnet.h>

static char* test_func_GEN_COST() {

  // Local variables
  Net* net;
  Vec* x;
  Func* f;

  printf("test_func_GEN_COST ...");
  
  // Load
  net = NET_new(1);
  NET_load(net,test_case,0);
  
  // Set flags
  NET_set_flags(net,OBJ_GEN,
		FLAG_VARS,
		GEN_PROP_ANY,
		GEN_VAR_P);
  
  Assert("error - wrong number of variables",NET_get_num_gens(net) == NET_get_num_vars(net));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  f = FUNC_GEN_COST_new(1.2,net);
  Assert("error - unable to create new function",f != NULL);
  Assert("error - bad function name",strcmp(FUNC_get_name(f),"generation cost") == 0);

  VEC_del(x);
  FUNC_del(f);
  NET_del(net);
  printf("ok\n");
  return 0;
}
  
