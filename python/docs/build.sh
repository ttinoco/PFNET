#! /bin/sh

cd _static
gunzip pfnet*
tar -xvf pfnet*
cd pfnet*
export PFNET=$PWD
./configure
make
cd ../../../
python setup.py build_ext --inplace --rpath=$PFNET/src
cd docs
