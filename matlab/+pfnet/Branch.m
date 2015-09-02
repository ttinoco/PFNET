classdef Branch < handle

  properties

    c_branch = libpointer;
    index = 0;

  end

  methods

    function br = Branch(ptr)
      br.c_branch = ptr;
    end

    % Getters and Setters
    %%%%%%%%%%%%%%%%%%%%%

    function idx = get.index(br)
      idx = calllib('libpfnet','BRANCH_get_index',br.c_branch);
    end
    
  end

end

