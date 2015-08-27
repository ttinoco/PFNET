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
  end

end

