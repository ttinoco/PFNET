<<<<<<< HEAD
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
    
=======
function load
 
  root = getenv('PFNET');
  header = fullfile(root,'include','pfnet','pfnet.h');
  lib = fullfile(root,'lib','libpfnet.so');

  if ~libisloaded('libpfnet')
     [notfound,warnings] = loadlibrary(lib,header,...
				       'addheader','net.h',...
				       'addheader','bus.h')
  end

  headers = {fullfile(root,'include','pfnet','obj_types.h'),...
             fullfile(root,'include','pfnet','flag_types.h')};

  for header = headers
    f = fopen(header{1},'rt');
    s = char(fread(f)');
    v = regexp(s,'#define\s+(\S+)\s+(\S+)','tokens');
    counter = 0;
    for e = v
       if counter ~= 0
        evalin('caller',strcat(e{1}{1},'=',e{1}{2},';'));
       end
       counter = counter+1;
    end
    fclose(f);
>>>>>>> Experimented with routine for loading #defines
  end

end

