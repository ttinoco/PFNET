#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cgraph

class GraphError(Exception):
    """
    Graph error exception.
    """

    pass

cdef class Graph:
    """
    Graph class.
    """

    cdef cgraph.Graph* _c_graph
    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self, net, alloc=True):
        """
        Graph class.

        Parameters
        ----------
        net : |Network|
        alloc : |TrueFalse|
        """

        pass

    def __cinit__(self, Network net, alloc=True):

        self._c_net = net._c_net
        if alloc:
            self._c_graph = cgraph.GRAPH_new(net._c_net)
        else:
            self._c_graph = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees graph C data structure.
        """

        if self.alloc:
            cgraph.GRAPH_del(self._c_graph)
            self._c_graph = NULL

    def has_viz(self):
        """
        Determines whether graph has visualization
        capabilities.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cgraph.GRAPH_can_viz(self._c_graph)

    def has_error(self):
        """
        Indicates whether the graph has the error flag set due to an
        invalid operation.
        """

        return cgraph.GRAPH_has_error(self._c_graph)

    def clear_error(self):
        """
        Clear error flag and message string.
        """

        cgraph.GRAPH_clear_error(self._c_graph);

    def set_layout(self):
        """
        Determines and saves a layout for the graph nodes.
        """

        cgraph.GRAPH_set_layout(self._c_graph)

    def set_node_property(self, bus, prop,value):
        """
        Sets property of node. See |Graphviz|.

        Parameters
        ----------
        bus : |Bus|
        prop : string
        value : string
        """

        cdef Bus b = bus
        cgraph.GRAPH_set_node_property(self._c_graph,b._c_ptr,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def set_nodes_property(self, prop, value):
        """
        Sets property of nodes. See |Graphviz|.

        Parameters
        ----------
        prop : string
        value : string
        """

        cgraph.GRAPH_set_nodes_property(self._c_graph,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def set_edges_property(self, prop, value):
        """
        Sets property of edges. See |Graphviz|. 

        Parameters
        ----------
        prop : string
        value : string
        """

        cgraph.GRAPH_set_edges_property(self._c_graph,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def color_nodes_by_mismatch(self, mis_type, t=0):
        """
        Colors the graphs nodes according to their power mismatch.

        Parameters
        ----------
        mis_type : string (mismatch attribute name)
        t : int
        """

        cgraph.GRAPH_color_nodes_by_mismatch(self._c_graph, str2mis_bus[mis_type], t)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def color_nodes_by_sensitivity(self, sens_type, t=0):
        """
        Colors the graphs nodes according to their sensitivity.

        Parameters
        ----------
        sens_type : string (sensitivity attribute name)
        t : int
        """

        cgraph.GRAPH_color_nodes_by_sensitivity(self._c_graph, str2sens_bus[sens_type], t)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def view(self, inline=False):
        """
        Displays the graph.

        Parameters
        ----------
        inline : |TrueFalse|
        """

        temp = tempfile.NamedTemporaryFile(delete=True, suffix='.png')
        try:
            self.write("png",temp.name)

            if inline is True:
                from IPython.display import Image

                self.write("png",temp.name)
                return Image(filename=temp.name)
            else:
                im = misc.imread(temp.name.encode('UTF-8'))
                misc.imshow(im)
        finally:
            temp.close()

    def write(self, format, filename):
        """
        Writes the graph to a file.

        Parameters
        ----------
        format : string (see |GraphvizOutputFormats|).
        filename : string
        """

        format = format.encode('UTF-8')
        filename = filename.encode('UTF-8')
        cgraph.GRAPH_write(self._c_graph,format,filename)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))
