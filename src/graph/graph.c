/** @file graph.c
 *  @brief This file defines the Graph data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/graph.h>

struct Graph {

  // Error
  BOOL error_flag;                      /**< @brief Error flags */
  char error_string[GRAPH_BUFFER_SIZE]; /**< @brief Error string */

  Net* net;
  BOOL layout_done;

  #ifndef NO_GRAPHVIZ
  Agraph_t* G;
  GVC_t* gvc;
  #endif
};

void GRAPH_del(Graph* g) {
  if (g) {

    #ifndef NO_GRAPHVIZ
    gvFreeLayout(g->gvc,g->G);
    agclose(g->G);
    gvFreeContext(g->gvc);
    #endif

    free(g);
  }
}

char* GRAPH_get_error_string(Graph* g) {
  if (!g)
    return NULL;
  else
    return g->error_string;
}

BOOL GRAPH_has_error(Graph* g) {
  if (!g)
    return FALSE;
  else
    return g->error_flag;
}

void GRAPH_clear_error(Graph* g) {
  if (g) {
    g->error_flag = FALSE;
    strcpy(g->error_string,"");
  }
}

Graph* GRAPH_new(Net* net) {

  Graph* g;

  #ifndef NO_GRAPHVIZ
  Branch* br;
  Bus* busk;
  Bus* busm;
  char buffer[100];
  Agnode_t* nodek;
  Agnode_t* nodem;
  int i;
  #endif

  // Init
  g = (Graph*)malloc(sizeof(Graph));
  g->net = net;
  g->layout_done = FALSE;
  g->error_flag = FALSE;
  strcpy(g->error_string,"");

  // Graph
  #ifndef NO_GRAPHVIZ
  g->gvc = gvContext();
  g->G = agopen("network",Agundirected,NULL);

  // Graph options
  agattr(g->G,AGRAPH,"size","10,7");
  agattr(g->G,AGRAPH,"ratio","true");
  agattr(g->G,AGRAPH,"overlap","false");
  agattr(g->G,AGRAPH,"splines","false");
  agattr(g->G,AGRAPH,"K","1.0");
  agattr(g->G,AGRAPH,"start","1");
  agattr(g->G,AGRAPH,"dpi","96");

  // Node options
  agattr(g->G,AGNODE,"shape","circle");
  agattr(g->G,AGNODE,"width","0.5");
  agattr(g->G,AGNODE,"height","0.5");
  agattr(g->G,AGNODE,"fixedsize","true");
  agattr(g->G,AGNODE,"label","");
  agattr(g->G,AGNODE,"style","filled");
  agattr(g->G,AGNODE,"fillcolor","white");
  agattr(g->G,AGNODE,"color","black");

  // Edge options
  agattr(g->G,AGEDGE,"style","filled");
  agattr(g->G,AGEDGE,"color","black");

  // Construct graph
  for (i = 0; i < NET_get_num_branches(net); i++) {

    br = NET_get_branch(net,i);

    busk = BRANCH_get_bus_k(br);
    sprintf(buffer,"%d",BUS_get_number(busk));
    nodek = agnode(g->G,strdup(buffer),TRUE);

    busm = BRANCH_get_bus_m(br);
    sprintf(buffer,"%d",BUS_get_number(busm));
    nodem = agnode(g->G,strdup(buffer),TRUE);

    sprintf(buffer,"%d",BRANCH_get_index(br));
    agedge(g->G,nodek,nodem,strdup(buffer),TRUE);
  }
  #endif

  return g;
}

BOOL GRAPH_can_viz(Graph* g) {
  if (g) {
    #ifndef NO_GRAPHVIZ
    return TRUE;
    #else
    return FALSE;
    #endif
  }
  return FALSE;
}

void GRAPH_set_layout(Graph* g) {
  if (g) {
    #ifndef NO_GRAPHVIZ
    gvFreeLayout(g->gvc,g->G);
    gvLayout(g->gvc,g->G,"sfdp");
    g->layout_done = TRUE;
    #endif
  }
}

void GRAPH_set_node_property(Graph* g, Bus* bus, char* prop, char* value) {

  #ifndef NO_GRAPHVIZ

  // Local variables
  Agsym_t* att;
  Agnode_t* node;
  char buffer[100];

  // Set property
  if (g) {
    sprintf(buffer,"%d",BUS_get_number(bus));
    node = agnode(g->G,buffer,FALSE);
    att = agattrsym(node,prop);
    if (att)
      agxset(node,att,value);
    else {
      sprintf(g->error_string,"invalid node property");
      g->error_flag = TRUE;
      return;
    }
  }
  #endif
}

void GRAPH_set_nodes_property(Graph* g, char* prop, char* value) {

  // Local variables
  int i;

  // Set property
  if (g) {
    for (i = 0; i < NET_get_num_buses(g->net); i++)
      GRAPH_set_node_property(g,NET_get_bus(g->net,i),prop,value);
  }
}

void GRAPH_set_edges_property(Graph* g, char* prop, char* value) {

  #ifndef NO_GRAPHVIZ

  // Local variables
  Agsym_t* att;
  Agedge_t* edge;
  Agnode_t* nodek;
  Bus* busk;
  Branch* branch;
  Agnode_t* nodem;
  Bus* busm;
  char buffer[100];
  int i;

  // Set property
  if (g) {
    for (i = 0; i < NET_get_num_branches(g->net); i++) {

      branch = NET_get_branch(g->net,i);

      busk = BRANCH_get_bus_k(branch);
      sprintf(buffer,"%d",BUS_get_number(busk));
      nodek = agnode(g->G,buffer,FALSE);

      busm = BRANCH_get_bus_m(branch);
      sprintf(buffer,"%d",BUS_get_number(busm));
      nodem = agnode(g->G,buffer,FALSE);

      sprintf(buffer,"%d",BRANCH_get_index(branch));
      edge = agedge(g->G,nodek,nodem,buffer,FALSE);
      att = agattrsym(edge,prop);
      if (att)
	agxset(edge,att,value);
      else {
	sprintf(g->error_string,"invalid edge property");
	g->error_flag = TRUE;
	return;
      }
    }
  }
  #endif
}

void GRAPH_color_nodes_by_mismatch(Graph* g, int mis_type, int t) {

  // Local variables
  #ifndef NO_GRAPHVIZ
  Agsym_t* att;
  Agnode_t* node;
  Bus* bus;
  REAL val;
  REAL mis;
  char buffer[100];
  char base[100];
  REAL eps = 1e-8;
  int i;
  #endif

  // Check graph and network
  if (!g || !g->net)
    return;

  // Check mismatch type
  if (mis_type < BUS_MIS_LARGEST || mis_type > BUS_MIS_REACTIVE) {
    sprintf(g->error_string,"invalid mismatch type");
    g->error_flag = TRUE;
    return;
  }

  // Update properties
  NET_update_properties(g->net,NULL);

  #ifndef NO_GRAPHVIZ

  // Color
  for (i = 0; i < NET_get_num_buses(g->net); i++) {

    bus = NET_get_bus(g->net,i);

    // Mismatch
    mis = BUS_get_quantity(bus,mis_type,t)*NET_get_base_power(g->net); // MW or MVAr

    // Value
    val = log10(fabs(mis) > eps ? fabs(mis) : eps);
    if (val > 3.)
      val = 3;
    if (val < -3)
      val = -3;
    val += 3;                // from 0 to 6
    if (mis < 0)
      val *= -1;             // from -6 to 6
    val = (5./6.)*(-val)+6.; // red = large pos, blue = large neg

    // Node
    sprintf(buffer,"%d",BUS_get_number(bus));
    node = agnode(g->G,buffer,FALSE);

    // Color scheme and value
    sprintf(buffer,"%d",(int)val);
    att = agattrsym(node,"fillcolor");
    sprintf(base,"/%s/",GRAPH_COLORMAP);
    agxset(node,att,strcat(base,buffer));
  }
  #endif
}

void GRAPH_color_nodes_by_sensitivity(Graph* g, int sens_type, int t) {

  // Local variables
  #ifndef NO_GRAPHVIZ
  Agsym_t* att;
  Agnode_t* node;
  Bus* bus;
  REAL sens;
  char buffer[100];
  char base[100];
  int i;
  int j;
  REAL max_sens = 1e-10;
  #endif

  // Check graph and network
  if (!g || !g->net)
    return;

  // Check sensitivity type
  if (sens_type < BUS_SENS_LARGEST || sens_type > BUS_SENS_V_REG_BY_SHUNT) {
    sprintf(g->error_string,"invalid sensitivity type");
    g->error_flag = TRUE;
    return;
  }

  // Update properties
  NET_update_properties(g->net,NULL);

  #ifndef NO_GRAPHVIZ

  // Color
  for (i = 0; i < 2; i++) {
    for (j = 0; j < NET_get_num_buses(g->net); j++) {

      bus = NET_get_bus(g->net,j);

      // Sensitivity
      sens = BUS_get_quantity(bus,sens_type,t);

      if (i == 0) {
	if (fabs(sens) > max_sens)
	  max_sens = fabs(sens);       // for normalizing
      }
      else {

	// Value
	sens = 5.*(-sens/max_sens)+6.; // red = large pos, blue = large neg

	// Node
	sprintf(buffer,"%d",BUS_get_number(bus));
	node = agnode(g->G,buffer,FALSE);

	// Color scheme and value
	sprintf(buffer,"%d",(int)sens);
	att = agattrsym(node,"fillcolor");
	sprintf(base,"/%s/",GRAPH_COLORMAP);
	agxset(node,att,strcat(base,buffer));
      }
    }
  }
  #endif
}

void GRAPH_write(Graph* g, char* format, char* filename) {
  if (!g) {
    sprintf(g->error_string,"graph missing");
    g->error_flag = TRUE;
    return;
  }

  if (!g->layout_done) {
    sprintf(g->error_string,"graph has no layout");
    g->error_flag = TRUE;
    return;
  }

  #ifndef NO_GRAPHVIZ
  if (gvRenderFilename(g->gvc,g->G,format,filename) == -1) {
    sprintf(g->error_string,"unable to render graph");
    g->error_flag = TRUE;
  }
  #endif
}
