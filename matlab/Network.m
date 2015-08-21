classdef Network < handle

  properties

    c_net = libpointer;
    num_buses = 0;

  end

  methods

    function net = Network()
      net.c_net = calllib('libpfnet','NET_new');
    end

    function delete(net)
      calllib('libpfnet','NET_del',net.c_net);
    end

    function bus = get_bus(net,i)
      bus = Bus(calllib('libpfnet','NET_get_bus',net.c_net,i));
    end

    function load(net,filename)
      calllib('libpfnet','NET_load',net.c_net,filename);
    end

    function show_components(net)
      calllib('libpfnet','NET_show_components',net.c_net);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function num = get.num_buses(net)
      num = calllib('libpfnet','NET_get_num_buses',net.c_net);
    end
    
  end

end

