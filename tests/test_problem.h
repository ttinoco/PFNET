/** @file test_problem.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/pfnet.h>

static char* test_problem_basic() {

  Parser* parser;
  Net* net;
  Prob* p;
  Vec* x;
  Mat* A;

  printf("test_problem_basic ...");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);

  Assert("error - invalid number of buses",NET_get_num_buses(net) > 0);

  p = PROB_new(net);
  
  Assert("error - bad prob net initialization",PROB_get_network(p) != NULL);

  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_NOT_REG_BY_GEN,
		BUS_VAR_VMAG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_NOT_SLACK,
		BUS_VAR_VANG);
  NET_set_flags(net,
		OBJ_GEN,
		FLAG_VARS,
		GEN_PROP_SLACK,
		GEN_VAR_P|GEN_VAR_Q);
  NET_set_flags(net,
		OBJ_GEN,
		FLAG_VARS,
		GEN_PROP_REG,
		GEN_VAR_Q);

  Assert("error - wrong number of vars",
	 NET_get_num_vars(net) == (NET_get_num_buses(net)-NET_get_num_buses_reg_by_gen(net) +
				   NET_get_num_buses(net)-NET_get_num_slack_buses(net)+
				   NET_get_num_slack_gens(net)+
				   NET_get_num_reg_gens(net)));

  PROB_add_constr(p,CONSTR_ACPF_new(net));
  PROB_add_constr(p,CONSTR_PAR_GEN_P_new(net));
  PROB_add_constr(p,CONSTR_PVPQ_SWITCHING_new(net));

  PROB_add_func(p,FUNC_REG_VMAG_new(3.4,net));

  Assert("error - cannot find constraint",PROB_find_constr(p,"AC power balance"));
  Assert("error - cannot find constraint",PROB_find_constr(p,"generator active power participation"));
  Assert("error - cannot find constraint",PROB_find_constr(p,"PVPQ switching"));
  Assert("error - finds nonexisting constraint",!PROB_find_constr(p,"voltage regulation by generators"));
  
  x = PROB_get_init_point(p);

  Assert("error - bad problem A init",PROB_get_A(p) == NULL);
  Assert("error - bad problem J init",PROB_get_J(p) == NULL);
  Assert("error - bad problem Hcomb init",PROB_get_H_combined(p) == NULL);
  Assert("error - bad problem f init",PROB_get_f(p) == NULL);

  PROB_analyze(p);

  Assert("error - bad problem A allocation",PROB_get_A(p) != NULL);
  Assert("error - bad problem J allocation",PROB_get_J(p) != NULL);
  Assert("error - bad problem f allocation",PROB_get_f(p) != NULL);
  Assert("error - bad problem b allocation",PROB_get_b(p) != NULL);
  Assert("error - bad problem Hcomb allocation",PROB_get_H_combined(p) != NULL);

  A = PROB_get_A(p);
  
  Assert("error - bad objective value init",PROB_get_phi(p) == 0.);

  Assert("error - bad A shape",MAT_get_size2(A) == NET_get_num_vars(net));
  
  PROB_eval(p,x);

  Assert("error - problem failed on eval",!PROB_has_error(p));
  Assert("error - bad objective value",PROB_get_phi(p) > 0.);

  VEC_del(x);
  PROB_del(p);
  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}
