.. _start:

***************
Getting Started
***************

.. _start_requirements:

Dependencies
============

The PFNET C library has the following *optional* dependencies:

* `Graphviz`_ (>= 2.38): for creating layouts and visualizations.
* Raw Parser (== 1.2.5): for parsing PSSE raw files.
* Line Flow: for constructing linear conservative AC branch flow limits.

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_build:

Build Instructions (Linux and Mac OS X)
=======================================

Building PFNET on Linux or Mac OS X requires typing the following commands in the root directory of the library::

  ./autogen.sh
  ./configure --prefix=$PWD/build
  make
  make check
  make install

For executing the command ``./autogen.sh``, Autotools is needed (m4, automake, autoconf, autoconf-archive, etc). 

If ``Raw Parser`` is available, the environment variable ``RAW_PARSER`` should point to the root of the source code. If ``Line Flow`` is available, the environment variable ``LINE_FLOW`` should be defined so that the library can be located in the directory ``$LINE_FLOW/lib``.

.. _start_build_win:

Build Instructions (Windows)
============================

Building PFNET on Windows requires `Cmake`_ and `MinGW`_, and typing the following commands in the root directory of the library::

  cmake -DCMAKE_INSTALL_PREFIX=.\build -G"MinGW Makefiles" -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE .
  mingw32-make -j
  mingw32-make install

.. _Graphviz: http://www.graphviz.org/
.. _Cmake: https://cmake.org/
.. _MinGW: http://www.mingw.org/
