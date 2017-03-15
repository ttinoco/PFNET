#! /bin/sh

cd _static
gunzip pfnet*
tar -xvf pfnet*
cd pfnet*
export PFNET=$PWD
echo $PFNET
./configure --prefix=$PWD/build
make
make install
ls $PFNET/build/lib
cd ../../../
./clean.sh
python setup.py build_ext --inplace --libdirs=$PFNET/build/lib
ls pfnet
python -c "import pfnet; assert True"
cd docs
