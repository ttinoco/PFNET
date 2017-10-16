.. include:: defs.hrst

.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET in Python. In particular, it covers dependencies, installation, and provides a simple example showing how to use this package.

.. _start_dependencies:

Dependencies
============

The PFNET Python module has the following dependencies:

* |Numpy| (>=1.11.2): the fundamental package for scientific computing in Python.
* |Scipy| (>=0.18.1): a collection of mathematical algorithms and functions built on top of Numpy.
* |Cython| (>=0.20.1): an optimizing static compiler for both Python and the extended Cython programming language.

.. _start_installation:

Installation
============

The PFNET Python module can be easily installed using the following commands::

  pip install numpy cython
  pip install pfnet  

After installation, the availability of optional features and the version of PFNET can be checked using::

  >>> import pfnet
  >>> pfnet.info
  {'line_flow': True, 'raw_parser': True, 'graphviz': True, 'version': '1.3.2'}
  
.. _start_example:

Example
=======

As a simple example of how to use the PFNET Python module, consider the task of constructing a power network from a |MATPOWER|-converted power flow file and computing the average bus degree. This can be done as follows::

  >>> import pfnet
  >>> import numpy as np

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> print np.average([bus.degree for bus in net.buses])
  2.86

In this example, is it assumed that the Python interpreter was started in a directory where the sample case |ieee14| is located.
