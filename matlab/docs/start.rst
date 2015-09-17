.. _start:

***************
Getting Started
***************

This section describes how to get started with PFNET in Matlab. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

Dependencies
============

PFNET for Matlab has the following dependencies:

* `PFNET <http://some_URL>`_: underlying C routines wrapped by this package (``libpfnet``).
* `Graphviz <http://www.graphviz.org/>`_ (>= 2.38): graph visualization library (``libgvc``) (Optional).
* `Raw parser <some_URL>`_ (>=1.0): library for parsing power flow files in PSSE raw format version 32 (``libraw_parser``) (Optional).

.. _start_download:

Download
========

The latest version of PFNET can be downloaded from `<https://github.com/ttinoco/PFNET>`_. Right now the Matlab wrapper is on the branch ``tomas-matlab``. 

.. _start_installation:

Installation (Linux)
====================

To use PFNET from Matlab, the location of the package needs to be added to Matlab's search path. This can be done within Matlab using::

  >>> addpath(strcat(getenv('PFNET'),'/matlab'));

where ``$PFNET`` is the root directory of the PFNET library. Then, the library needs to be loaded using::

  >>> pfnet.load_library

The environment variable ``$PFNET`` is needed by this routine to find the required header files and shared library. The command::

  >>> libfunctions libpfnet

can be used to list all the loaded functions of PFNET and hence to check whether the library was loaded successfully.

If PFNET was built with visualization capabilities, then Matlab needs to be loaded with the Graphviz shared libraries ``libcgraph.so`` and ``libgvc.so``. This can be done by starting Matlab using the command::

  > LD_PRELOAD=${GRAPHVIZ}/lib/libcgraph.so:${GRAPHVIZ}/lib/libgvc.so matlab

where ``$GRAPHVIZ`` is the Graphviz installation directory, or more conveniently, by defining::

  alias matlab='LD_PRELOAD=${GRAPHVIZ}/lib/libcgraph.so:${GRAPHVIZ}/lib/libgvc.so matlab'

and then starting Matlab with the redefined command ``matlab``. 

.. _start_example:

Example
=======

As a quick example of how to use the PFNET Matlab package, consider the task of constructing a power network from a `MATPOWER <http://www.pserc.cornell.edu//matpower/>`_-converted power flow file and computing the average bus degree. This can be done as follows::
  
  >>> addpath(strcat(getenv('PFNET'),'/matlab'));

  >>> pfnet.load_library

  >>> net = pfnet.Network();
  >>> net.load('ieee14.mat');

  >>> deg = 0;
  >>> for i=1:net.num_buses
  >>>     bus = net.get_bus(i-1);
  >>>     deg = deg + bus.degree/net.num_buses;
  >>> end
  >>> disp(deg);
  2.86







