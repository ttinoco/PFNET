classdef Network < handle

  properties

    c_net = libpointer;

    base_power = 0;

    num_buses = 0;
    num_gens = 0;
    num_branches = 0;
    num_loads = 0;
    num_shunts = 0;
    num_vars = 0;
   
  end

  methods

    function net = Network()
      net.c_net = calllib('libpfnet','NET_new');
    end

    function delete(net)
      calllib('libpfnet','NET_del',net.c_net);
    end

    function bus = get_bus(net,i)
      bus = pfnet.Bus(calllib('libpfnet','NET_get_bus',net.c_net,i));
    end

    function br = get_branch(net,i)
      br = pfnet.Branch(calllib('libpfnet','NET_get_branch',net.c_net,i));
    end

    function gen = get_gen(net,i)
      gen = pfnet.Generator(calllib('libpfnet','NET_get_gen',net.c_net,i));
    end

    function load = get_load(net,i)
      load = pfnet.Load(calllib('libpfnet','NET_get_load',net.c_net,i));
    end

    function shunt = get_shunt(net,i)
      shunt = pfnet.Shunt(calllib('libpfnet','NET_get_shunt',net.c_net,i));
    end

    function load(net,filename)
      calllib('libpfnet','NET_load',net.c_net,filename);
    end

    function show_components(net)
      calllib('libpfnet','NET_show_components',net.c_net);
    end

    function show_properties(net)
      calllib('libpfnet','NET_show_properties',net.c_net);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function num = get.num_buses(net)
      num = calllib('libpfnet','NET_get_num_buses',net.c_net);
    end

    function num = get.num_gens(net)
      num = calllib('libpfnet','NET_get_num_gens',net.c_net);
    end

    function num = get.num_branches(net)
      num = calllib('libpfnet','NET_get_num_branches',net.c_net);
    end

    function num = get.num_loads(net)
      num = calllib('libpfnet','NET_get_num_loads',net.c_net);
    end

    function num = get.num_shunts(net)
      num = calllib('libpfnet','NET_get_num_shunts',net.c_net);
    end
    
  end

end

