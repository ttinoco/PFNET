#! /bin/sh

if [ ! -d "lib/pfnet" ]; then
  cd lib
  gunzip -c pfnet*.tar.gz > pfnet.tar
  tar -xvf pfnet.tar
  rm pfnet.tar
  mv pfnet*/ pfnet
  cd pfnet
  if [ "$(uname)" == "Darwin" ]; then
      ./configure --prefix=$PWD/build LDFLAGS="-Wl,-install_name,@loader_path/libpfnet.dylib"
  else
      ./configure --prefix=$PWD/build
  fi
  make clean
  make uninstall
  make
  make check
  make install
  cp build/lib/* ../../pfnet/
fi
