/** @file graph.h
 *  @brief This file lists the constants and routines associated with the Graph data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __GRAPH_HEADER__
#define __GRAPH_HEADER__

#include <stdlib.h>
#include "net.h"

#ifndef NO_GRAPHVIZ
#include <graphviz/gvc.h>
#endif

#define GRAPH_COLORMAP "spectral11" /**< @brief Default graph colormap */
#define GRAPH_BUFFER_SIZE 1024      /**< @brief Default graph buffer size for strings */

typedef struct Graph Graph;

void GRAPH_color_nodes_by_mismatch(Graph* g, int mis_type);
void GRAPH_color_nodes_by_sensitivity(Graph* g, int sens_type);
void GRAPH_del(Graph* g);
void GRAPH_clear_error(Graph* g);
char* GRAPH_get_error_string(Graph* g);
BOOL GRAPH_has_error(Graph* g);
Graph* GRAPH_new(Net* net);
BOOL GRAPH_can_viz(Graph* g);
void GRAPH_set_layout(Graph* g);
void GRAPH_set_nodes_property(Graph* g, char* prop, char* value);
void GRAPH_set_edges_property(Graph* g, char* prop, char* value);
void GRAPH_write(Graph* g, char* format, char* filename);

#endif


