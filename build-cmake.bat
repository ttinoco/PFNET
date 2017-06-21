mkdir ./cmake-build
cd ./cmake-build
  cmake -DCMAKE_INSTALL_PREFIX=..\cmake-install -DRAW_PARSER_SOURCE_DIR=C:\Users\pawi002\ROMCON\raw-parser\ -G"MinGW Makefiles" -DPFNET_GRAPHVIZ=OFF -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE ..
  mingw32-make -j
  mingw32-make install
cd ..
