/** @file test_graph.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include <pfnet/net.h>
#include <pfnet/graph.h>

static char* test_graph_basic() {
  
  Graph* g;
  Net* net;

  printf("test_graph_basic ...");

  net = NET_new(1);

  NET_load(net,test_case,0);

  g = GRAPH_new(net);
  Assert("error - unable to create graph",g != NULL);
  Assert("error - bad error flag initialization",!GRAPH_has_error(g));

  GRAPH_set_layout(g);
  Assert("error - graph layout failed",!GRAPH_has_error(g));

  GRAPH_del(g);
  NET_del(net);
  printf("ok\n");
  return 0;
}
