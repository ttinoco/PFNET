/** @file test_constraints.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/pfnet.h>

static char* test_constr_BOUND() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  int num;
  Vec* f;
  Mat* J;

  printf("test_constr_BOUND ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set flags
  NET_set_flags(net,OBJ_BUS,
		FLAG_VARS|FLAG_BOUNDED,
		BUS_PROP_ANY,
		BUS_VAR_VMAG|BUS_VAR_VANG);
  NET_set_flags(net,OBJ_GEN,
		FLAG_VARS|FLAG_BOUNDED,
		GEN_PROP_ANY,
		GEN_VAR_P|GEN_VAR_Q);
  NET_set_flags(net,OBJ_BRANCH,
		FLAG_VARS|FLAG_BOUNDED,
		BRANCH_PROP_TAP_CHANGER,BRANCH_VAR_RATIO);
  NET_set_flags(net,OBJ_BRANCH,
		FLAG_VARS|FLAG_BOUNDED,
		BRANCH_PROP_PHASE_SHIFTER,
		BRANCH_VAR_PHASE);
  NET_set_flags(net,OBJ_SHUNT,
		FLAG_VARS|FLAG_BOUNDED,
		SHUNT_PROP_SWITCHED_V,
		SHUNT_VAR_SUSC);
  
  num = (2*NET_get_num_buses(net) +
	 2*NET_get_num_gens(net) +
	 NET_get_num_tap_changers(net) +
	 NET_get_num_phase_shifters(net) +
	 NET_get_num_switched_shunts(net));

  Assert("error - empty network",num > 0);
  Assert("error - wrong number of variables",num == NET_get_num_vars(net));
  Assert("error - wrong number of bounded quantities",num == NET_get_num_bounded(net));
  
  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  c = CONSTR_new(CONSTR_TYPE_BOUND,net);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 0);
  CONSTR_count(c);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 2*num);  
  CONSTR_allocate(c);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 2*num);
  CONSTR_analyze(c);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 2*num);
  CONSTR_eval(c,x);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 2*num);
  CONSTR_store_sens(c,NULL,NULL,NULL,NULL);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 2*num);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == 0);

  f = CONSTR_get_f(c);
  J = CONSTR_get_J(c);

  Assert("error - NULL f",f != NULL);
  Assert("error - NULL J",J != NULL);
  Assert("error - bad f size", VEC_get_size(f) == 2*num);
  Assert("error - bad J size", MAT_get_size1(J) == 2*num);
  Assert("error - bad J size", MAT_get_size2(J) == NET_get_num_vars(net));
  Assert("error - bad J size", MAT_get_nnz(J) == 2*num);
  
  CONSTR_clear(c);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 0);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_FIX() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  int num;
  Vec* b;
  Mat* A;

  printf("test_constr_FIX ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set flags
  NET_set_flags(net,OBJ_BUS,
		FLAG_VARS|FLAG_FIXED,
		BUS_PROP_ANY,
		BUS_VAR_VMAG|BUS_VAR_VANG);
  NET_set_flags(net,OBJ_GEN,
		FLAG_VARS|FLAG_FIXED,
		GEN_PROP_ANY,
		GEN_VAR_P|GEN_VAR_Q);
  NET_set_flags(net,OBJ_BRANCH,
		FLAG_VARS|FLAG_FIXED,
		BRANCH_PROP_TAP_CHANGER,
		BRANCH_VAR_RATIO);
  NET_set_flags(net,OBJ_BRANCH,
		FLAG_VARS|FLAG_FIXED,
		BRANCH_PROP_PHASE_SHIFTER,
		BRANCH_VAR_PHASE);
  NET_set_flags(net,OBJ_SHUNT,
		FLAG_VARS|FLAG_FIXED,
		SHUNT_PROP_SWITCHED_V,
		SHUNT_VAR_SUSC);

  num = (2*NET_get_num_buses(net) +
	 2*NET_get_num_gens(net) +
	 NET_get_num_tap_changers(net) +
	 NET_get_num_phase_shifters(net) +
	 NET_get_num_switched_shunts(net));

  Assert("error - empty network",num > 0);
  Assert("error - wrong number of variables",num == NET_get_num_vars(net));
  Assert("error - wrong number of fixed quantities",num == NET_get_num_fixed(net));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  c = CONSTR_new(CONSTR_TYPE_FIX,net);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == 0);
  CONSTR_count(c);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == num+NET_get_num_buses_reg_by_gen(net));  
  CONSTR_allocate(c);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == num+NET_get_num_buses_reg_by_gen(net));
  CONSTR_analyze(c);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == num+NET_get_num_buses_reg_by_gen(net));
  CONSTR_eval(c,x);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == 0);
  CONSTR_store_sens(c,NULL,NULL,NULL,NULL);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 0);

  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);

  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size1(A) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad A size", MAT_get_nnz(A) == num+NET_get_num_buses_reg_by_gen(net));
  
  CONSTR_clear(c);
  Assert("error - wrong Annz counter",CONSTR_get_Acounter(c) == 0);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_PAR_GEN_P() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  int num;
  int nnz;
  Vec* b;
  Mat* A;
  Bus* bus;
  int i;

  printf("test_constr_PAR_GEN_P ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set flags
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
  
  Assert("error - empty network",NET_get_num_vars(net) > 0);
  Assert("error - wrong number of variables",
	 NET_get_num_vars(net) == (NET_get_num_buses(net)-NET_get_num_buses_reg_by_gen(net)+
				   NET_get_num_buses(net)-NET_get_num_slack_buses(net)+
				   NET_get_num_slack_gens(net)+
				   NET_get_num_reg_gens(net)));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  num = 0;
  nnz = 0;
  for (i = 0; i < NET_get_num_buses(net); i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_slack(bus)) {
      num += BUS_get_num_gens(bus)-1;
      nnz += 2*(BUS_get_num_gens(bus)-1);
    }
  }

  c = CONSTR_new(CONSTR_TYPE_PAR_GEN_P,net);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);
  CONSTR_count(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_allocate(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_analyze(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_eval(c,x);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);  
  CONSTR_store_sens(c,NULL,NULL,NULL,NULL);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);  
  
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);

  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size1(A) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad A size", MAT_get_nnz(A) == nnz);
  
  CONSTR_clear(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_PAR_GEN_Q() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  int num;
  int nnz;
  Vec* b;
  Mat* A;
  Bus* bus;
  int i;

  printf("test_constr_PAR_GEN_Q ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set flags
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
  
  Assert("error - empty network",NET_get_num_vars(net) > 0);
  Assert("error - wrong number of variables",
	 NET_get_num_vars(net) == (NET_get_num_buses(net)-NET_get_num_buses_reg_by_gen(net)+
				   NET_get_num_buses(net)-NET_get_num_slack_buses(net)+
				   NET_get_num_slack_gens(net)+
				   NET_get_num_reg_gens(net)));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  num = 0;
  nnz = 0;
  for (i = 0; i < NET_get_num_buses(net); i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_regulated_by_gen(bus)) {
      num += BUS_get_num_reg_gens(bus)-1;
      nnz += 2*(BUS_get_num_reg_gens(bus)-1);
    }       
  }

  c = CONSTR_new(CONSTR_TYPE_PAR_GEN_Q,net);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  Assert("error - unable to create new constraint",c != NULL);
  Assert("error - bad constraint initialization",CONSTR_get_b(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_A(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_f(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_J(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_array(c) == NULL);
  Assert("error - bad constraint initialization",CONSTR_get_H_combined(c) == NULL);
  
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);
  CONSTR_count(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_allocate(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_analyze(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == nnz);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == num);  
  CONSTR_eval(c,x);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);  
  CONSTR_store_sens(c,NULL,NULL,NULL,NULL);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);  
  
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);

  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size1(A) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad A size", MAT_get_nnz(A) == nnz);
  
  CONSTR_clear(c);
  Assert("error - wrong A counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - wrong A counter",CONSTR_get_Aconstr_index(c) == 0);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_PF() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  Vec* f;
  Mat* J;
  Mat* H;
  int Jnnz_computed;
  int Hnnz;
  int Hnnz_computed;
  int* Hcounter;
  int size;
  int i;

  printf("test_constr_PF ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_ANY,
		BUS_VAR_VMAG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_ANY,
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
  NET_set_flags(net,
		OBJ_BRANCH,
		FLAG_VARS,
		BRANCH_PROP_TAP_CHANGER,
		BRANCH_VAR_RATIO);
  NET_set_flags(net,
		OBJ_BRANCH,
		FLAG_VARS,
		BRANCH_PROP_PHASE_SHIFTER,
		BRANCH_VAR_PHASE);
  NET_set_flags(net,
		OBJ_SHUNT,
		FLAG_VARS,
		SHUNT_PROP_SWITCHED_V,
		SHUNT_VAR_SUSC);

  Assert("error - empty network",NET_get_num_vars(net) > 0);
  Assert("error - wrong number of variables",
	    NET_get_num_vars(net) == (2*NET_get_num_buses(net)+
				      NET_get_num_slack_gens(net)+
				      NET_get_num_reg_gens(net)+
				      NET_get_num_tap_changers(net)+
				      NET_get_num_phase_shifters(net)+
				      NET_get_num_switched_shunts(net)));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  c = CONSTR_new(CONSTR_TYPE_PF,net);
  
  Jnnz_computed = (NET_get_num_buses(net)*4 +
		   NET_get_num_branches(net)*8 +
		   NET_get_num_tap_changers(net)*4 +
		   NET_get_num_phase_shifters(net)*4 +
		   NET_get_num_switched_shunts(net) +
		   NET_get_num_slack_gens(net) +
		   NET_get_num_reg_gens(net));

  Hnnz_computed = (NET_get_num_buses(net)*3 +
		   NET_get_num_branches(net)*12 +
		   NET_get_num_tap_changers(net)*9 +
                   NET_get_num_phase_shifters(net)*10 +
		   NET_get_num_switched_shunts(net));
		     
  CONSTR_count(c);
  
  Hcounter = CONSTR_get_Hcounter(c);
  size = CONSTR_get_Hcounter_size(c);
  Hnnz = 0;
  for (i = 0; i < size; i++)
    Hnnz += Hcounter[i];
  
  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",Jnnz_computed == CONSTR_get_Jcounter(c));
  Assert("error - bad Hnnz counter",Hnnz_computed == Hnnz);

  CONSTR_allocate(c);
  
  CONSTR_analyze(c);

  Hnnz = 0;
  for (i = 0; i < size; i++)
    Hnnz += Hcounter[i];

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",Jnnz_computed == CONSTR_get_Jcounter(c));
  Assert("error - bad Hnnz counter",Hnnz_computed == Hnnz);

  CONSTR_eval(c,x);
  Hnnz = 0;
  for (i = 0; i < size; i++)
    Hnnz += Hcounter[i];

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",Jnnz_computed == CONSTR_get_Jcounter(c));
  Assert("error - bad Hnnz counter",Hnnz_computed == Hnnz);

  f = CONSTR_get_f(c);
  J = CONSTR_get_J(c);
  H = CONSTR_get_H_single(c,0);

  Assert("error - NULL f",f != NULL);
  Assert("error - NULL J",J != NULL);
  Assert("error - bad f size", VEC_get_size(f) == 2*NET_get_num_buses(net));
  Assert("error - bad J size", MAT_get_size1(J) == 2*NET_get_num_buses(net));
  Assert("error - bad J size", MAT_get_size2(J) == NET_get_num_vars(net));
  Assert("error - bad J size", MAT_get_nnz(J) == Jnnz_computed);
  Assert("error - bad H size", MAT_get_size1(H) == NET_get_num_vars(net));
  Assert("error - bad H size", MAT_get_size2(H) == NET_get_num_vars(net));
  
  CONSTR_clear(c);
  Assert("error - wrong Jnnz counter",CONSTR_get_Jcounter(c) == 0);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_REG_GEN() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  Vec* b;
  Mat* A;
  Vec* f;
  Mat* J;
  int num_vars;
  int num_Annz;
  int num_Jnnz;
  int num;
  int i;
  Bus* bus;

  printf("test_constr_REG_GEN ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_NOT_SLACK,
		BUS_VAR_VMAG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_NOT_SLACK,
		BUS_VAR_VANG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_NOT_SLACK|BUS_PROP_REG_BY_GEN,
		BUS_VAR_VDEV);
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
  
  num_vars = (2*(NET_get_num_buses(net)-NET_get_num_slack_buses(net))+
	      2*(NET_get_num_buses_reg_by_gen(net)-NET_get_num_slack_buses(net))+
	      NET_get_num_slack_gens(net)+
	      NET_get_num_reg_gens(net));
  Assert("error - invalid number of varibles",num_vars == NET_get_num_vars(net));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));

  c = CONSTR_new(CONSTR_TYPE_REG_GEN,net);

  num = NET_get_num_buses_reg_by_gen(net)-NET_get_num_slack_buses(net);
  num_Annz = 3*num;
  num_Jnnz = 0;
  for (i = 0; i < NET_get_num_buses(net); i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus))
      num_Jnnz += 2 + 2*BUS_get_num_reg_gens(bus);
  }

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == 0);
  
  CONSTR_count(c);
  
  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 2*num);

  CONSTR_allocate(c);

  CONSTR_analyze(c);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 2*num);

  CONSTR_eval(c,x);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == 0);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 2*num);

  A = CONSTR_get_A(c);
  b = CONSTR_get_b(c);
  f = CONSTR_get_f(c);
  J = CONSTR_get_J(c);

  Assert("error - NULL f",f != NULL);
  Assert("error - NULL J",J != NULL);
  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);

  Assert("error - bad f size", VEC_get_size(f) == 2*num);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad J size", MAT_get_size2(J) == NET_get_num_vars(net));

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_REG_TRAN() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  Vec* b;
  Mat* A;
  Vec* f;
  Mat* J;
  int num_vars;
  int num_Annz;
  int num_Jnnz;
  int num;
  int i;
  Bus* bus;

  printf("test_constr_REG_TRAN ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_REG_BY_TRAN,
		BUS_VAR_VMAG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_REG_BY_TRAN,
		BUS_VAR_VVIO);
  NET_set_flags(net,
		OBJ_BRANCH,
		FLAG_VARS,
		BRANCH_PROP_TAP_CHANGER_V,
		BRANCH_VAR_RATIO);
  NET_set_flags(net,
		OBJ_BRANCH,
		FLAG_VARS,
		BRANCH_PROP_TAP_CHANGER_V,
		BRANCH_VAR_RATIO_DEV);
  
  num_vars = (3*NET_get_num_buses_reg_by_tran(net)+
	      3*NET_get_num_tap_changers_v(net));
  Assert("error - invalid number of varibles",num_vars == NET_get_num_vars(net));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));
  
  c = CONSTR_new(CONSTR_TYPE_REG_TRAN,net);

  num = NET_get_num_tap_changers_v(net);
  num_Annz = 3*num;
  num_Jnnz = 10*num;

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == 0);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == 0);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 0);
  
  CONSTR_count(c);
  
  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  CONSTR_allocate(c);

  CONSTR_analyze(c);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  CONSTR_eval(c,x);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == 0);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  A = CONSTR_get_A(c);
  b = CONSTR_get_b(c);
  f = CONSTR_get_f(c);
  J = CONSTR_get_J(c);

  Assert("error - NULL f",f != NULL);
  Assert("error - NULL J",J != NULL);
  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);

  Assert("error - bad f size", VEC_get_size(f) == 4*num);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad A size", MAT_get_nnz(A) == num_Annz);
  Assert("error - bad J size", MAT_get_size2(J) == NET_get_num_vars(net));
  Assert("error - bad J size", MAT_get_nnz(J) == num_Jnnz);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_constr_REG_SHUNT() {
  
  // Local variables
  Net *net;
  Vec* x;
  Constr* c;
  Vec* b;
  Mat* A;
  Vec* f;
  Mat* J;
  int num_vars;
  int num_Annz;
  int num_Jnnz;
  int num;
  int i;
  Bus* bus;

  printf("test_constr_REG_SHUNT ...");

  // Load
  net = NET_new();
  NET_load(net,test_case,0);

  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_REG_BY_SHUNT,
		BUS_VAR_VMAG);
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_REG_BY_SHUNT,
		BUS_VAR_VVIO);
  NET_set_flags(net,
		OBJ_SHUNT,
		FLAG_VARS,
		SHUNT_PROP_SWITCHED_V,
		SHUNT_VAR_SUSC);
  NET_set_flags(net,
		OBJ_SHUNT,
		FLAG_VARS,
		SHUNT_PROP_SWITCHED_V,
		SHUNT_VAR_SUSC_DEV);
  
  num_vars = (3*NET_get_num_buses_reg_by_shunt(net)+
	      3*NET_get_num_switched_shunts(net));
  Assert("error - invalid number of varibles",num_vars == NET_get_num_vars(net));

  x = NET_get_var_values(net,CURRENT);

  Assert("error - NULL vector of var values",x != NULL);
  Assert("error - vector of var values has wrong shape",VEC_get_size(x) == NET_get_num_vars(net));
  
  c = CONSTR_new(CONSTR_TYPE_REG_SHUNT,net);

  num = NET_get_num_switched_shunts(net);
  num_Annz = 3*num;
  num_Jnnz = 10*num;

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == 0);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == 0);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 0);
  
  CONSTR_count(c);
  
  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  CONSTR_allocate(c);

  CONSTR_analyze(c);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == num_Annz);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == num);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  CONSTR_eval(c,x);

  Assert("error - bad Annz counter",CONSTR_get_Acounter(c) == 0);
  Assert("error - bad Jnnz counter",CONSTR_get_Jcounter(c) == num_Jnnz);
  Assert("error - bad Aindex counter",CONSTR_get_Aconstr_index(c) == 0);
  Assert("error - bad Jindex counter",CONSTR_get_Jconstr_index(c) == 4*num);

  A = CONSTR_get_A(c);
  b = CONSTR_get_b(c);
  f = CONSTR_get_f(c);
  J = CONSTR_get_J(c);

  Assert("error - NULL f",f != NULL);
  Assert("error - NULL J",J != NULL);
  Assert("error - NULL b",b != NULL);
  Assert("error - NULL A",A != NULL);

  Assert("error - bad f size", VEC_get_size(f) == 4*num);
  Assert("error - bad b size", VEC_get_size(b) == num);
  Assert("error - bad A size", MAT_get_size2(A) == NET_get_num_vars(net));
  Assert("error - bad A size", MAT_get_nnz(A) == num_Annz);
  Assert("error - bad J size", MAT_get_size2(J) == NET_get_num_vars(net));
  Assert("error - bad J size", MAT_get_nnz(J) == num_Jnnz);

  VEC_del(x);
  CONSTR_del(c);
  NET_del(net);
  printf("ok\n");
  return 0;
}
