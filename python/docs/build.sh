#! /bin/sh

cd ../../
apt-get install autoconf
apt-get install automake
./autogen.sh
./configure
make
make install
cd python
python setup.py build_ext --inplace
cd docs
