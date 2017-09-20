IF NOT EXIST "lib\pfnet" (
  cd lib
  7z x pfnet*.tar.gz
  move pfnet*.tar pfnet.tar
  7z x pfnet.tar
  del pfnet.tar
  for /d %%G in ("pfnet*") do move "%%~G" pfnet
  cd pfnet
  cmake -DCMAKE_INSTALL_PREFIX=.\build -G"MinGW Makefiles" -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE .
  mingw32-make -j
  mingw32-make install
  copy build\lib\* ..\..\pfnet\
)
