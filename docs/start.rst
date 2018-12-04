.. _start:

***************
Getting Started
***************

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_install:

Installation
============

.. _start_install_unix:

Linux and macOS
---------------

Installing PFNET on Linux or macOS requires typing the following commands in the root directory of the library::

  ./autogen.sh
  ./configure --prefix=$PWD/build
  make
  make check
  make install

For executing the command ``./autogen.sh``, Autotools is needed (m4, automake, autoconf, autoconf-archive, etc). 

.. _start_install_win:

Windows
-------

Installation PFNET on Windows requires `Cmake`_ and `MinGW`_, and typing the following commands in the root directory of the library::

  cmake -DCMAKE_INSTALL_PREFIX=.\build -G"MinGW Makefiles" -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE .
  mingw32-make -j
  mingw32-make install

.. _Graphviz: http://www.graphviz.org/
.. _Cmake: https://cmake.org/
.. _MinGW: http://www.mingw.org/
