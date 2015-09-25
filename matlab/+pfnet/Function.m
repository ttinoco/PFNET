classdef Function < handle

  properties

    c_func = libpointer;
    type = 0;
    weight = 0;
    phi = 0;
    gphi = [];
    Hphi = [];

  end

  methods

    function func = Function(type,weight,net)
      func.c_func = calllib('libpfnet','FUNC_new',type,weight,net.c_net);
    end

    function analyze(func)
      calllib('libpfnet','FUNC_count',func.c_func);
      calllib('libpfnet','FUNC_allocate',func.c_func);
      calllib('libpfnet','FUNC_analyze',func.c_func);
    end

    function eval(func,var_values)
      assert(size(var_values,1)==1);
      v = calllib('libpfnet','VEC_new_from_array',var_values,size(var_values,2));
      calllib('libpfnet','FUNC_eval',func.c_func,v);
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function t = get.type(func)
      t = calllib('libpfnet','FUNC_get_type',func.c_func);
    end

    function w = get.weight(func)
      w = calllib('libpfnet','FUNC_get_weight',func.c_func);
    end

    function phi = get.phi(func)
      phi = calllib('libpfnet','FUNC_get_phi',func.c_func);
    end

    function gphi = get.gphi(func)
      gphi = pfnet.Vector(calllib('libpfnet','FUNC_get_gphi',func.c_func));
    end
    
  end

end

