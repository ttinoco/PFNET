#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/graph.h":

    ctypedef struct Graph
    ctypedef struct Bus
    ctypedef struct Net

    void GRAPH_color_nodes_by_mismatch(Graph* g, int mis_type)
    void GRAPH_color_nodes_by_sensitivity(Graph* g, int sens_type)
    void GRAPH_clear_error(Graph* g)
    void GRAPH_del(Graph* g)
    char* GRAPH_get_error_string(Graph* g)
    bint GRAPH_has_error(Graph* g)
    bint GRAPH_can_viz(Graph* g)
    Graph* GRAPH_new(Net* net)
    void GRAPH_set_layout(Graph* g)
    void GRAPH_set_node_property(Graph* g, Bus* bus, char* prop, char* value)
    void GRAPH_set_edges_property(Graph* g, char* prop, char* value)
    void GRAPH_set_nodes_property(Graph* g, char* prop, char* value)
    void GRAPH_write(Graph* g, char* format, char* filename)


    
