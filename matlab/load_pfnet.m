function load_pfnet
  root = getenv('PFNET');
  header = fullfile(root,'include','pfnet','pfnet.h');
  lib = fullfile(root,'lib','libpfnet.so');
  if ~libisloaded('libpfnet')
     [notfound,warnings] = loadlibrary(lib,header,'addheader','net.h','addheader','bus.h')
  end
  %libfunctions libpfnet
end
