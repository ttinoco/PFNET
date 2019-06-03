cmake -DCMAKE_INSTALL_PREFIX=.\build -G"MinGW Makefiles" -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE .
mingw32-make -j
mingw32-make install
