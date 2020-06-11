./clean.sh
./autogen.sh
./configure --prefix="$PWD/build"
make clean
make
make check
make install
make dist
yes | cp pfnet*.gz ../PFNET.py/lib
