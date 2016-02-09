/** @file test_network.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/net.h>

static char* test_net_new() {

  Net* net;

  printf("test_net_new ... ");

  net = NET_new();

  Assert("error - failed to create net",net != NULL);

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_load() {

  Net* net;

  printf("test_net_load ... ");

  net = NET_new();
  Assert("error - invalid number of buses",NET_get_num_buses(net) == 0);
  NET_load(net,test_case);
  Assert("error - failed to parse case",!NET_has_error(net));
  Assert("error - invalid number of buses",NET_get_num_buses(net) > 0);

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_check() {
  
  Net* net;

  printf("test_net_check ... ");

  net = NET_new();				
  NET_load(net,test_case);
  Assert("error - net check failed",NET_check(net,0));
  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_variables() {
  
  int num = 0;
  Net* net;

  printf("test_net_variables ... ");
  
  net = NET_new();				

  NET_load(net,test_case);

  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);
  
  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_SLACK,BUS_VAR_VMAG);
  num += NET_get_num_slack_buses(net);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_ANY,BUS_VAR_VANG);
  num += NET_get_num_buses(net);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_ANY,BUS_VAR_VANG);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_VARS,GEN_PROP_SLACK,GEN_VAR_P|GEN_VAR_Q);
  num += 2*NET_get_num_slack_gens(net);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_VARS,GEN_PROP_SLACK,GEN_VAR_P|GEN_VAR_Q);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_VARS,BRANCH_PROP_TAP_CHANGER,BRANCH_VAR_RATIO);
  num += NET_get_num_tap_changers(net);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_VARS,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += NET_get_num_phase_shifters(net);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_fixed() {

  int num = 0;
  Net* net;

  printf("test_net_fixed ... ");

  net = NET_new();	
			
  NET_load(net,test_case);
  
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_FIXED,BUS_PROP_REG_BY_GEN,BUS_VAR_VMAG);
  num += NET_get_num_buses_reg_by_gen(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);
  
  NET_set_flags(net,OBJ_BUS,FLAG_FIXED,BUS_PROP_REG_BY_GEN,BUS_VAR_VMAG);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_FIXED,GEN_PROP_SLACK,GEN_VAR_P);
  num += NET_get_num_slack_gens(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_FIXED,GEN_PROP_SLACK,GEN_VAR_Q);
  num += NET_get_num_slack_gens(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += NET_get_num_phase_shifters(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_TAP_CHANGER,BRANCH_VAR_RATIO);
  num += NET_get_num_tap_changers(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_SHUNT,FLAG_FIXED,SHUNT_PROP_SWITCHED_V,SHUNT_VAR_SUSC);
  num += NET_get_num_switched_shunts(net);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_properties() {

  Bus* bus;
  Gen* gen;
  REAL busvmax = 0;
  REAL genQvio = 0;
  REAL dQ;
  Net* net;
  int i;

  printf("test_net_properties ... ");

  net = NET_new();

  NET_load(net,test_case);

  for (i = 0; i < NET_get_num_buses(net); i++) {
    
    bus = NET_get_bus(net,i);
    
    if (BUS_get_v_mag(bus) > busvmax)
      busvmax = BUS_get_v_mag(bus);

  }

  for (i = 0; i < NET_get_num_gens(net); i++) {

    gen = NET_get_gen(net,i);

    if (GEN_is_regulator(gen)) {
      
      dQ = 0;
      if (GEN_get_Q(gen) > GEN_get_Q_max(gen))
	dQ = GEN_get_Q(gen)-GEN_get_Q_max(gen);
      if (GEN_get_Q(gen) < GEN_get_Q_min(gen))
	dQ = GEN_get_Q_min(gen)-GEN_get_Q(gen);
      if (dQ > genQvio)
	genQvio = dQ;
    }
  }
  
  Assert("error - bad network property (bus_v_max)",busvmax == NET_get_bus_v_max(net));
  Assert("error - bad network property (gen_Q_vio)",genQvio*NET_get_base_power(net) == NET_get_gen_Q_vio(net));

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_init_point() {

  Bus* bus;
  Gen* gen;
  Net* net;
  Vec* point;
  int i;

  printf("test_net_init_point ... ");

  net = NET_new();

  NET_load(net,test_case);
  
  // Set variables
  NET_set_flags(net,
		OBJ_BUS,
		FLAG_VARS,
		BUS_PROP_ANY,
		BUS_VAR_VMAG|BUS_VAR_VANG);
  NET_set_flags(net,
		OBJ_GEN,
		FLAG_VARS,
		GEN_PROP_SLACK,
		GEN_VAR_P|GEN_VAR_Q);
  
  Assert("error - wrong number of vars",
	 NET_get_num_vars(net) == (2*NET_get_num_buses(net)+
				   2*NET_get_num_slack_gens(net)));

  point = NET_get_var_values(net,CURRENT);

  Assert("error - NULL init point",point != NULL);
  Assert("error - invalid vector size",VEC_get_size(point) == NET_get_num_vars(net));

  for (i = 0; i < NET_get_num_buses(net); i++) {
    bus = NET_get_bus(net,i);
    Assert("error - bad vars flag of bus",BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG|BUS_VAR_VANG));
    Assert("error - bad voltage init point",BUS_get_v_mag(bus) == VEC_get(point,BUS_get_index_v_mag(bus)));
    Assert("error - bad voltage init point",BUS_get_v_ang(bus) == VEC_get(point,BUS_get_index_v_ang(bus)));
  }
  for (i = 0; i < NET_get_num_gens(net); i++) {
    gen = NET_get_gen(net,i);
    if (GEN_is_slack(gen)) {
      Assert("error - bad vars flag of gen",GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P|GEN_VAR_Q));
      Assert("error - bad gen init point",GEN_get_P(gen) == VEC_get(point,GEN_get_index_P(gen)));
      Assert("error - bad gen init point",GEN_get_Q(gen) == VEC_get(point,GEN_get_index_Q(gen)));
    }
    else {
      Assert("error - bad vars flag of gen",!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P));
      Assert("error - bad vars flag of gen",!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q));
    }
  }

  VEC_del(point);
  NET_del(net);
  printf("ok\n");
  return 0;
}

