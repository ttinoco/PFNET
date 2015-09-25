classdef Bus < handle

  properties

    c_bus = libpointer;
    index = 0;
    index_v_mag = 0;
    index_v_ang = 0;
    number = 0;
    degree = 0;
    v_mag = 0;
    v_ang = 0;
    v_set = 0;
    v_max = 0;
    v_min = 0;
    gens = {};

  end

  methods

    function bus = Bus(ptr)
      bus.c_bus = ptr;
    end

    function flag = is_slack(bus)
      flag = calllib('libpfnet','BUS_is_slack',bus.c_bus);
    end

    function flag = is_regulated_by_gen(bus)
      flag = calllib('libpfnet','BUS_is_regulated_by_gen',bus.c_bus);
    end

    function flag = has_flags(bus,fmask,vmask)
      flag = calllib('libpfnet','BUS_has_flags',bus.c_bus,fmask,vmask);
    end

    function P = get_total_gen_P(bus)
      P = calllib('libpfnet','BUS_get_total_gen_P',bus.c_bus);
    end

    function Q = get_total_gen_Q(bus)
      Q = calllib('libpfnet','BUS_get_total_gen_Q',bus.c_bus);
    end

    function P = get_total_load_P(bus)
      P = calllib('libpfnet','BUS_get_total_load_P',bus.c_bus);
    end

    function Q = get_total_load_Q(bus)
      Q = calllib('libpfnet','BUS_get_total_load_Q',bus.c_bus);
    end

    function show(bus)
      calllib('libpfnet','BUS_show',bus.c_bus);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(bus)
      idx = calllib('libpfnet','BUS_get_index',bus.c_bus);
    end

    function idx = get.index_v_mag(bus)
      idx = calllib('libpfnet','BUS_get_index_v_mag',bus.c_bus);
    end

    function idx = get.index_v_ang(bus)
      idx = calllib('libpfnet','BUS_get_index_v_ang',bus.c_bus);
    end
    
    function num = get.number(bus)
      num = calllib('libpfnet','BUS_get_number',bus.c_bus);
    end

    function deg = get.degree(bus)
      deg = calllib('libpfnet','BUS_get_degree',bus.c_bus);
    end

    function v = get.v_mag(bus)
      v = calllib('libpfnet','BUS_get_v_mag',bus.c_bus);
    end

    function v = get.v_ang(bus)
      v = calllib('libpfnet','BUS_get_v_ang',bus.c_bus);
    end

    function v = get.v_set(bus)
      v = calllib('libpfnet','BUS_get_v_set',bus.c_bus);
    end

    function v = get.v_max(bus)
      v = calllib('libpfnet','BUS_get_v_max',bus.c_bus);
    end

    function v = get.v_min(bus)
      v = calllib('libpfnet','BUS_get_v_min',bus.c_bus);
    end

    function gens = get.gens(bus)
      g = calllib('libpfnet','BUS_get_gen',bus.c_bus);
      len = calllib('libpfnet','GEN_list_len',g);
      gens = cell(len,1);
      i = 1;
      while ~isNull(g)
        gens{i} = pfnet.Generator(g);
        g = calllib('libpfnet','GEN_get_next',g);
        i = i + 1;
      end
    end      
    
  end

end

