classdef Load < handle

  properties

    c_load = libpointer;
    index = 0;

  end

  methods

    function load = Load(ptr)
      load.c_load = ptr;
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(load)
      idx = calllib('libpfnet','LOAD_get_index',load.c_load);
    end
    
  end

end

