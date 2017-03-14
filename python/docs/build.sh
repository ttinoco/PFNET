#! /bin/sh

cd ../../
./configure
make
cd python
python setup.py build_ext --inplace --rpath=../src
cd docs
