classdef Network < handle

  properties

    c_net = libpointer;

  end

  methods

    function net = Network()
      net.load_library();
      net.c_net = calllib('libpfnet','NET_new');
    end

    function show_components(net)
      calllib('libpfnet','NET_show_components',net.c_net)
    end
    
  end

  methods (Static)

    function load_library()
      root = getenv('PFNET');
      header = fullfile(root,'include','pfnet','net.h');
      lib = fullfile(root,'lib','libpfnet.so');
      if ~libisloaded('libpfnet')
         [notfound,warnings] = loadlibrary(lib,header);
      end
    end
        
  end

end

