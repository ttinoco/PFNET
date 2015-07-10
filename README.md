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

* [Graphviz](http://www.graphviz.org/) (``libgvc``>= 2.36) (Optional)
* Raw parser (``libraw_parser``) (Optional)

## Build Instructions ##

To build the library, type the command ``make`` in the root directory of the package. This will create the file ``libpfnet.so`` in the ``lib`` directory.

To build the library without visualization capabilities (no ``libgvc`` dependency), pass ``NO_GRAPHVIZ=1`` to ``make``.

To build the library without raw parsing capabilities (no ``libraw_parser`` dependency), pass ``NO_RAW_PARSER=1`` to ``make``.

## Wrappers ##

Wrappers for PFNET are available for the following languages:

* [Python](http://ttinoco.github.io/PFNET/python)

## Citing PFNET ##

If you use PFNET in your work, please cite the software as follows:

\code
@misc{pfnet,
  author={Tomas Tinoco De Rubira},
  title={{PFNET}: A library for modeling and analyzing electric power networks},
  howpublished={\url{https://github.com/ttinoco/PFNET}},
  month={July},
  year={2015}
}
\endcode

## Contact ##

* Tomas Tinoco De Rubira (<ttinoco5687@gmail.com>)