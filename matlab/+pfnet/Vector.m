function vec = Vector(v)

  data = calllib('libpfnet','VEC_get_data',v);
  size = calllib('libpfnet','VEC_get_size',v);

  setdatatype(data,'doublePtr',1,size);
  vec = data.Value;
end
