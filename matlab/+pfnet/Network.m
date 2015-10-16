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

    buses = [];
    branches = [];
    generators = [];
    loads = [];
    shunts = [];
    var_generators = [];
   
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

    function bus = get_bus_by_number(net,i)
      bus = pfnet.Bus(calllib('libpfnet','NET_bus_hash_find',net.c_net,i));
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

    function values = get_var_values(net,code)
      if (~exist('code','var'))
	code=0;
      end
      values = pfnet.Vector(calllib('libpfnet','NET_get_var_values',net.c_net,code));
    end

    function gbuses = get_gen_buses(net)
      gbuses = [];
      ptr = calllib('libpfnet','NET_get_gen_buses',net.c_net);
      while ~isNull(ptr)
        bus = pfnet.Bus(ptr);
        gbuses = [bus gbuses];
        ptr = calllib('libpfnet','BUS_get_next',ptr);
      end
    end

    function num = get_num_buses_reg_by_gen(net)
      num = calllib('libpfnet','NET_get_num_buses_reg_by_gen',net.c_net);
    end

    function num = get_num_lines(net)
      num = calllib('libpfnet','NET_get_num_lines',net.c_net);
    end

    function num = get_num_slack_gens(net)
      num = calllib('libpfnet','NET_get_num_slack_gens',net.c_net);
    end

    function load(net,filename)
      calllib('libpfnet','NET_load',net.c_net,filename);
    end

    function set_flags(net,obj_type,flags,props,vals);
      calllib('libpfnet','NET_set_flags',net.c_net,obj_type,flags,props,vals);
    end

    function show_components(net)
      calllib('libpfnet','NET_show_components',net.c_net);
    end

    function show_properties(net)
      calllib('libpfnet','NET_show_properties',net.c_net);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function bp = get.base_power(net)
      bp = calllib('libpfnet','NET_get_base_power',net.c_net);
    end

    function num = get.num_vars(net)
      num = calllib('libpfnet','NET_get_num_vars',net.c_net);
    end

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

    function buses = get.buses(net)
      buses = [];
      for i=1:net.num_buses
          bus = net.get_bus(i-1);
          buses = [bus buses];
      end
    end

    function branches = get.branches(net)
      branches = [];
      for i=1:net.num_branches
          br = net.get_branch(i-1);
          branches = [br branches];
      end
    end

    function gens = get.generators(net)
      gens = [];
      for i=1:net.num_gens
          g = net.get_gen(i-1);
          gens = [g gens];
      end
    end

    function loads = get.loads(net)
      loads = [];
      for i=1:net.num_loads
          l = net.get_load(i-1);
          loads = [l loads];
      end
    end

    function shunts = get.shunts(net)
      shunts = [];
      for i=1:net.num_shunts
          s = net.get_shunt(i-1);
          shunts = [s shunts];
      end
    end

     function vgens = get.var_generators(net)
      vgens = [];
      for i=1:net.num_vargens
          vg = net.get_vargen(i-1);
          vgens = [vg vgens];
      end
    end  
    
  end

end

