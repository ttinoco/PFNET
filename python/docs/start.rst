.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET in Python. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

Dependencies
============

PFNET for Python has the following dependencies:

* `Numpy <http://www.numpy.org>`_ (>=1.11.2): the fundamental package for scientific computing in Python.
* `Scipy <http://www.scipy.org>`_ (>=0.18.1): a collection of mathematical algorithms and functions built on top of Numpy.
* `Cython <http://cython.org>`_ (>=0.20.1): an optimizing static compiler for both Python and the extended Cython programming language.
* `PFNET <https://github.com/ttinoco/PFNET>`_ (== 1.2.8): underlying C routines wrapped by this package.

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_installation:

Installation
============

After building the C library, the PFNET Python module can be installed using::

  sudo pip install -r requirements.txt
  sudo python setup.py install

from the ``python`` directory of the PFNET library. The module can be tested using::

  python setup.py build_ext --inplace
  nosetests -s -v

For making PFNET parsers written in Python work seamlessly with the C library, the C library and Python module need to be built twice. The availability of optional features of PFNET can be checked as in the following example::

  >>> import pfnet
  >>> pfnet.info
  {'python parsers': True, 'raw parser': False, 'graphviz': True}

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

Requirements for building the PFNET Python documentation locally:

* `Sphinx <http://www.sphinx-doc.org/>`_ (>=1.4).

To build the documentation, the environment variable ``PFNET_DOCS`` must be set. The generated files will be placed in the directory ``PFNET_DOCS/python``. To generate the files, run ``make html`` from the ``python/docs`` directory of the PFNET library.
