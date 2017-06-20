/** @file test_graph.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/parser.h>
#include <pfnet/net.h>
#include <pfnet/graph.h>

static char* test_graph_basic() {

  Parser* parser;
  Graph* g;
  Net* net;

  printf("test_graph_basic ...");

  parser = PARSER_new_for_file(test_case);
  net = PARSER_parse(parser,test_case,1);
  
  Assert("error - bad number of buses",NET_get_num_buses(net) > 0);

  g = GRAPH_new(net);
  Assert("error - unable to create graph",g != NULL);
  Assert("error - bad error flag initialization",!GRAPH_has_error(g));

  GRAPH_set_layout(g);
  Assert("error - graph layout failed",!GRAPH_has_error(g));

  GRAPH_del(g);
  NET_del(net);
  PARSER_del(parser);
  printf("ok\n");
  return 0;
}
