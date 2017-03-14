#! /bin/sh

cd ../../
./autogen.sh
./configure
make
make install
cd python/docs
python setup.py build_ext --inplace
