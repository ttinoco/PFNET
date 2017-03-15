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
cd ../../../
./clean.sh
ls $PFNET/build/lib
python setup.py build_ext --inplace --libdirs=$PFNET/build/lib
ls pfnet
python -c "import pfnet; assert True"
cd docs
