del CMakeCache.txt
del CTestTestfile.cmake
del cmake_install.cmake
del install_manifest.txt
del pfnet_static_tests.exe
del pfnet_tests.exe
del Makefile
rmdir /s /q %~dp0build
rmdir /s /q %~dp0CMakeFiles
del /s /f %~dp0*.dll
del /s /f %~dp0*.a





