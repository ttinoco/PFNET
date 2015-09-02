classdef Shunt < handle

  properties

    c_shunt = libpointer;
    index = 0;

  end

  methods

    function shunt = Shunt(ptr)
      shunt.c_shunt = ptr;
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(shunt)
      idx = calllib('libpfnet','SHUNT_get_index',shunt.c_shunt);
    end
    
  end

end

