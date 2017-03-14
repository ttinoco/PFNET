#! /bin/sh

cd ../../
./configure
make
make install
cd python
python setup.py build_ext --inplace
cd docs
