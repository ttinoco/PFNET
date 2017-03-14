#! /bin/sh

cd ../../
aclocal
autoconf
autoheader
automake --add-missing
./configure
make
cd python
python setup.py build_ext --inplace --rpath=../src
cd docs
