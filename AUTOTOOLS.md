Procedure for using autoconf and automake
-----------------------------------------
* run autoscan to generate configure.scan
* renable to configure.ac
* run **autoheader** to generate config.h.in
* make source portable by looking at config.h.in
* create Makefile.am src/Makefile.am tests/Makefile.am
* run libtoolize
* run automake --add-missing
* run **automake**
* run **aclocal**
* run **autoconf**
* configure, make, make check, make install
