classdef Generator < handle

  properties

    c_gen = libpointer;
    index = 0;
    P = 0;
    Q = 0;
    bus = pfnet.Bus.empty;

  end

  methods

    function gen = Generator(ptr)
      gen.c_gen = ptr;
    end

    function flag = is_slack(gen)
      flag = calllib('libpfnet','GEN_is_slack',gen.c_gen);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(gen)
      idx = calllib('libpfnet','GEN_get_index',gen.c_gen);
    end

    function p = get.P(gen)
      p = calllib('libpfnet','GEN_get_P',gen.c_gen);
    end

    function q = get.Q(gen)
      q = calllib('libpfnet','GEN_get_Q',gen.c_gen);
    end

    function b = get.bus(gen)
      b = pfnet.Bus(calllib('libpfnet','GEN_get_bus',gen.c_gen));
    end
    
  end

end

