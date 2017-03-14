#! /bin/sh

cd _static
gunzip pfnet*
tar -xvf pfnet*
cd pfnet*
./configure
make
cd ../../../
python setup.py build_ext --inplace --rpath=./docs/_static
cd docs
