.. include:: defs.hrst

.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET in Python. In particular, it covers prerequisites and installation, and provides a simple example that shows how to use this package.

.. _start_prerequisites:

Prerequisites
=============

Before installing the PFNET Python module, the following tools are needed:

* Linux and macOS: a C compiler, |Make|, |Python| and |pip|.
* Windows : |Anaconda| (for Python 2.7), |CMake|, |7-Zip|, and |MinGW|.

.. _start_installation:

Installation
============

After the prerequisites for the appropriate operating system have been obtained, the PFNET Python module can be easily installed by executing the following commands on the terminal or Anaconda prompt::

  pip install numpy cython
  pip install pfnet

After installation, the availability of optional features and the version of PFNET can be checked in Python using::

  >>> import pfnet
  >>> pfnet.info
  {'raw_parser': True, 'graphviz': True, 'version': '1.3.3'}

To install the module from source, the code can be obtained from `<https://github.com/ttinoco/PFNET>`_, and then the following commands can be executed on the terminal or Anaconda prompt from the ``python`` directory of the package::

    pip install numpy cython
    python setup.py install

Running the unit tests can be done with::

    python setup.py build_ext --inplace
    nosetests -s -v

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
