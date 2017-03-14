#! /bin/sh

cd ../../
sudo apt-get install autoconf
sudo apt-get install automake
./autogen.sh
./configure
make
sudo make install
cd python
python setup.py build_ext --inplace
cd docs
