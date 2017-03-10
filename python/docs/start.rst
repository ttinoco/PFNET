.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET in Python. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

Dependencies
============

PFNET for Python has the following dependencies:

* `Numpy <http://www.numpy.org>`_ (>=1.8.2): the fundamental package for scientific computing in Python.
* `Scipy <http://www.scipy.org>`_ (>=0.13.3): a collection of mathematical algorithms and functions built on top of Numpy.
* `Cython <http://cython.org>`_ (>=0.20.1): an optimising static compiler for both Python and the extended Cython programming language.
* `PFNET <https://github.com/ttinoco/PFNET>`_ (== 1.2.7): underlying C routines wrapped by this package (``libpfnet``).
* `Graphviz <http://www.graphviz.org/>`_ (>= 2.38): graph visualization library (``libgvc``) (Optional).
* `Raw parser <some_URL>`_ (>=1.2.1): library for parsing power flow files in PSSE raw format version 32 (``libraw_parser``) (Optional).

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_installation:

Installation
============

After building the C library ``libpfnet``, the PFNET Python module can be installed using::

  > sudo python setup.py install

from the ``python`` directory of the PFNET package.

If ``libpfnet`` was built without visualization capabilities, the argument ``--no_graphviz`` should be passed to ``setup.py``. Similarly, if ``libpfnet`` was build without raw parsing capabilities, the argument ``--no_raw_parser`` should be passed to ``setup.py``.

The installation can be tested using `nose <https://nose.readthedocs.org/en/latest/>`_ as follows::

  > python setup.py build_ext --inplace
  > nosetests -v --exe

.. _start_example:

Example
=======

As a quick example of how to use the PFNET Python module, consider the task of constructing a power network from a `MATPOWER <http://www.pserc.cornell.edu//matpower/>`_-converted power flow file and computing the average bus degree. This can be done as follows::

  >>> import numpy as np
  >>> from pfnet import Network

  >>> net = Network()
  >>> net.load('ieee14.mat')

  >>> print np.average([b.degree for b in net.buses])
  2.86

Documentation
=============

Requirements to build the PFNET Python documentation:

* `Sphinx <http://www.sphinx-doc.org/>`_ (>=1.4).

To build the documentation, the environment variable ``PFNET_DOCS`` must be set. The generated files will be placed in the directory ``PFNET_DOCS/python``. To generate the files, run ``make html`` from the ``python/docs`` directory of the PFNET package.

It may also be necessary to pass the environment variable with the path to the dynamic shared libraries using ``LD_LIBRARY_PATH`` on Linux or ``DYLD_FALLBACK_LIBRARY_PATH`` on Mac OSX. The command would then be::

  > make html DYLD_FALLBACK_LIBRARY_PATH=$PFNET/lib
