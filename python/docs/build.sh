#! /bin/sh

cd _static
rm -rf pfnet*/
rm pfnet*.tar
gunzip pfnet*
tar -xvf pfnet*
cd pfnet*
export PFNET=$PWD
echo $PFNET
./configure --prefix=$PWD/build
make clean
make uninstall
make
make check
make install
ls $PFNET/build/lib
ls $PFNET/build/include/pfnet
cd ../../../
./clean.sh
python setup.py build_ext --inplace --libdirs=$PFNET/build/lib --incdirs=$PFNET/build/include
ls pfnet
python -c "import pfnet; assert True"
cd docs
