/** @file test_network.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/parser.h>
#include <pfnet/net.h>

static char* test_net_new() {

  Net* net;

  printf("test_net_new ... ");

  net = NET_new(1);

  Assert("error - failed to create net",net != NULL);

  NET_del(net);
  printf("ok\n");
  return 0;
}

static char* test_net_load() {
  
  Parser* parser;
  Net* net;

  printf("test_net_load ... ");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);
  Assert("error - unable to get parser",parser != NULL);
  Assert(PARSER_get_error_string(parser),!PARSER_has_error(parser));
  Assert("error - failed to parse case",!PARSER_has_error(parser));
  Assert("error - invalid number of buses",NET_get_num_buses(net,FALSE) > 0);

  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}

static char* test_net_check() {
  
  Parser* parser;
  Net* net;

  printf("test_net_check ... ");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);
  Assert("error - net check failed",NET_check(net,0));
  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}

static char* test_net_variables() {
  
  int num = 0;
  Parser* parser;
  Net* net;

  printf("test_net_variables ... ");
  
  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);
  
  Assert("error - bad number of buses",NET_get_num_buses(net,FALSE) > 0);

  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);
  
  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_SLACK,BUS_VAR_VMAG);
  num += NET_get_num_slack_buses(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_ANY,BUS_VAR_VANG);
  num += NET_get_num_buses(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_VARS,BUS_PROP_ANY,BUS_VAR_VANG);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_VARS,GEN_PROP_SLACK,GEN_VAR_P|GEN_VAR_Q);
  num += 2*NET_get_num_slack_gens(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_VARS,GEN_PROP_SLACK,GEN_VAR_P|GEN_VAR_Q);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_VARS,BRANCH_PROP_TAP_CHANGER,BRANCH_VAR_RATIO);
  num += NET_get_num_tap_changers(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_VARS,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += NET_get_num_phase_shifters(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_vars(net) == num);

  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}

static char* test_net_fixed() {

  int num = 0;
  Parser* parser;
  Net* net;

  printf("test_net_fixed ... ");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);

  Assert("error - bad number of buses",NET_get_num_buses(net,FALSE) > 0);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BUS,FLAG_FIXED,BUS_PROP_REG_BY_GEN,BUS_VAR_VMAG);
  num += NET_get_num_buses_reg_by_gen(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);
  
  NET_set_flags(net,OBJ_BUS,FLAG_FIXED,BUS_PROP_REG_BY_GEN,BUS_VAR_VMAG);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_FIXED,GEN_PROP_SLACK,GEN_VAR_P);
  num += NET_get_num_slack_gens(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_GEN,FLAG_FIXED,GEN_PROP_SLACK,GEN_VAR_Q);
  num += NET_get_num_slack_gens(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += NET_get_num_phase_shifters(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_PHASE_SHIFTER,BRANCH_VAR_PHASE);
  num += 0;
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_BRANCH,FLAG_FIXED,BRANCH_PROP_TAP_CHANGER,BRANCH_VAR_RATIO);
  num += NET_get_num_tap_changers(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_set_flags(net,OBJ_SHUNT,FLAG_FIXED,SHUNT_PROP_SWITCHED_V,SHUNT_VAR_SUSC);
  num += NET_get_num_switched_shunts(net,TRUE);
  Assert("error - wrong number of variables",NET_get_num_fixed(net) == num);

  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}

static char* test_net_properties() {

  Bus* bus;
  Gen* gen;
  REAL busvmax = 0;
  REAL genQvio = 0;
  REAL dQ;
  Parser* parser;
  Net* net;
  int i;

  printf("test_net_properties ... ");
  
  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);

  Assert("error - invalid number of buses",NET_get_num_buses(net,FALSE) > 0);

  for (i = 0; i < NET_get_num_buses(net,FALSE); i++) {
    
    bus = NET_get_bus(net,i);
    
    if (BUS_get_v_mag(bus,0) > busvmax && BUS_is_in_service(bus))
      busvmax = BUS_get_v_mag(bus,0);

  }

  for (i = 0; i < NET_get_num_gens(net,FALSE); i++) {

    gen = NET_get_gen(net,i);

    if (GEN_is_regulator(gen) && !GEN_is_slack(gen) && GEN_is_in_service(gen)) {
      
      dQ = 0;
      if (GEN_get_Q(gen,0) > GEN_get_Q_max(gen))
        dQ = GEN_get_Q(gen,0)-GEN_get_Q_max(gen);
      if (GEN_get_Q(gen,0) < GEN_get_Q_min(gen))
        dQ = GEN_get_Q_min(gen)-GEN_get_Q(gen,0);
      if (dQ > genQvio)
        genQvio = dQ;
    }
  }
  
  Assert("error - bad network property (bus_v_max)",busvmax == NET_get_bus_v_max(net,0));
  Assert("error - bad network property (gen_Q_vio)",genQvio*NET_get_base_power(net) == NET_get_gen_Q_vio(net,0));

  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}

static char* test_net_init_point() {

  Bus* bus;
  Gen* gen;
  Parser* parser;
  Net* net;
  Vec* point;
  int i;

  printf("test_net_init_point ... ");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);

  Assert("error - invalid number of buses",NET_get_num_buses(net,FALSE) > 0); 

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
         NET_get_num_vars(net) == (2*NET_get_num_buses(net,TRUE)+
                                   2*NET_get_num_slack_gens(net,TRUE)));

  point = NET_get_var_values(net,CURRENT);

  Assert("error - NULL init point",point != NULL);
  Assert("error - invalid vector size",VEC_get_size(point) == NET_get_num_vars(net));

  for (i = 0; i < NET_get_num_buses(net,FALSE); i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_in_service(bus)) {
      Assert("error - bad vars flag of bus",BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG|BUS_VAR_VANG));
      Assert("error - bad voltage init point",BUS_get_v_mag(bus,0) == VEC_get(point,BUS_get_index_v_mag(bus,0)));
      Assert("error - bad voltage init point",BUS_get_v_ang(bus,0) == VEC_get(point,BUS_get_index_v_ang(bus,0)));
    }
  }
  for (i = 0; i < NET_get_num_gens(net,FALSE); i++) {
    gen = NET_get_gen(net,i);
    if (GEN_is_slack(gen) && GEN_is_in_service(gen)) {
      Assert("error - bad vars flag of gen",GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P|GEN_VAR_Q));
      Assert("error - bad gen init point",GEN_get_P(gen,0) == VEC_get(point,GEN_get_index_P(gen,0)));
      Assert("error - bad gen init point",GEN_get_Q(gen,0) == VEC_get(point,GEN_get_index_Q(gen,0)));
    }
    else if (GEN_is_in_service(gen)) {
      Assert("error - bad vars flag of gen",!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P));
      Assert("error - bad vars flag of gen",!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q));
    }
  }

  VEC_del(point);
  PARSER_del(parser);
  NET_del(net);
  printf("ok\n");
  return 0;
}

