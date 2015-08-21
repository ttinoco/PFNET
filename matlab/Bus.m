classdef Bus < handle

  properties

    c_bus = libpointer;
    index = 0;

  end

  methods

    function bus = Bus(ptr)
      bus.c_bus = ptr;
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(bus)
      idx = calllib('libpfnet','BUS_get_index',bus.c_bus);
    end
    
  end

end

