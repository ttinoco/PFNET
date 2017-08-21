#! /bin/sh

if [ ! -d "lib/pfnet" ]; then
  cd lib
  gunzip -c pfnet*.tar.gz > pfnet.tar
  tar -xvf pfnet.tar
  rm pfnet.tar
  mv pfnet*/ pfnet
  cd pfnet
  ./configure --prefix=$PWD/build
  make clean
  make uninstall
  make
  make check
  make install
  cp build/lib/* ../../pfnet/
fi
