# PFNET: Power Flow Network Library #

## Overview ##

PFNET is a C library for modeling and analyzing electric power networks. It provides the following:

* Parsers for power network data files.
* Network visualization routines.
* Fast and customizable constraint and objective function evaluators for network optimization problems.

## License ##

BSD 2-clause license.

## Documentation ##

The documentation for this library can be found in http://ttinoco.github.io/PFNET.

## Download ##

The latest version of the library can be downloaded from https://github.com/ttinoco/PFNET.

## Dependencies ##

* [Graphviz](http://www.graphviz.org/) (``libgvc``>= 2.38) (Optional)
* Raw parser (``libraw_parser``) (Optional)

## Build Instructions (Linux) ##

Building PFNET requires typing the command ``make`` in the root directory of the library, say ``$PFNET``. This creates the shared library ``libpfnet.so`` in the ``$PFNET/lib`` directory.

To build the library without visualization capabilities (no ``libgvc`` dependency), ``make`` should be passed the argument ``NO_GRAPHVIZ=1``. Otherwise, the compiler needs to find the header ``graphviz/gvc.h``, and the linker needs to find ``libgvc.so`` (to build the tests). If Graphviz has been installed in some non-standard directory, then ``$GRAPHVIZ`` should be set to that directory so that the compiler can find the above header in ``$GRAPHVIZ/include``. The directory ``$GRAPHVIZ/lib`` can be added to ``$LD_LIBRARY_PATH`` so that the linker can find the above shared library.

To build the library without raw parsing capabilities (no ``libraw_parser`` dependency), ``make`` should be passed the argument ``NO_RAW_PARSER=1``.

## Wrappers ##

Wrappers for PFNET are available for the following languages:

* [Python](http://ttinoco.github.io/PFNET/python)
* [Matlab](http://ttinoco.github.io/PFNET/matlab) (getting started)

## Contributors ##

* [Tomas Tinoco De Rubira](http://n.ethz.ch/~tomast/) (principal developer)