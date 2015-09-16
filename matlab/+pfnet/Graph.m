classdef Graph < handle

  properties
    c_graph = libpointer; 
  end

  methods

    function graph = Graph(net)
      graph.c_graph = calllib('libpfnet','GRAPH_new',net.c_net);   
    end

    function delete(graph)
      calllib('libpfnet','GRAPH_del',graph.c_graph);
    end

    function set_layout(graph)
      calllib('libpfnet','GRAPH_set_layout',graph.c_graph);
    end

    function write(graph,format,filename)
      calllib('libpfnet','GRAPH_write',graph.c_graph,format,filename);
    end

  end
end
