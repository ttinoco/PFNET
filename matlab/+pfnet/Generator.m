classdef Generator < handle

  properties

    c_gen = libpointer;
    index = 0;

  end

  methods

    function gen = Generator(ptr)
      gen.c_gen = ptr;
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(gen)
      idx = calllib('libpfnet','GEN_get_index',gen.c_gen);
    end
    
  end

end

