.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

Dependencies
============

PFNET has the following dependencies:

* `Numpy <http://www.numpy.org>`_ (>=1.8.2): the fundamental package for scientific computing in Python.
* `Scipy <http://www.scipy.org>`_ (>=0.13.3): a collection of mathematical algorithms and functions built on top of Numpy.
* `PFNET <http://some_URL>`_: underlying C routines wrapped by this package (``libpfnet``).
* `Graphviz <http://www.graphviz.org/>`_ (>= 2.36): graph visualization library (``libgvc``) (Optional). 
* `Raw parser <some_URL>`_ (>=1.0): library for parsing power flow files in PSSE raw format version 32 (``libraw_parser``) (Optional).

.. _start_download:

Download
========

The latest version of PFNET can be downloaded `here <some_URL>`_.

.. _start_installation:

Installation
============

After building the C library ``libpfnet``, the PFNET Python module can be installed using::

  > sudo python setup.py install

from the ``python`` directory of the PFNET package.

If ``libpfnet`` was built without visualization capabilities, the argument ``--no_graphviz`` should be passed to ``setup.py``. Similarly, if ``libpfnet`` was build without raw parsing capabilities, the argument ``--no_raw_parser`` should be passed to ``setup.py``.

The installation can be tested using `nose <https://nose.readthedocs.org/en/latest/>`_ as follows::

  > nosetests -v

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
