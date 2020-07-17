#! /bin/bash

./clean.sh
./autogen.sh
./configure --prefix="$PWD/build"
make clean
make
make check
make install
make dist
