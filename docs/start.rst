.. _start:

***************
Getting Started
***************

.. _start_requirements:

Dependencies
============

The PFNET C library has the following optional dependencies:

* `Graphviz`_ (libgvc >= 2.38): for creating network layouts and visualizing networks.
* Raw Parser (libraw_parser >= 1.2.4) : for parsing PSSE raw files.
* Line Flow (libline_flow): for constructing linear conservative AC branch flow limits.

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_install:

Installation (Linux and Mac OS X)
=================================

Installing PFNET requires typing the following commands in the root directory of the library::

  ./autogen.sh
  ./configure
  make
  make check
  sudo make install

For executing the command ``./autogen.sh`` you need Autotools (m4, automake, autoconf, autoconf-archive, etc). 

If ``Raw Parser`` is available, the library should be placed in the ``lib`` directory of PFNET. If ``Line Flow`` is available, the environment variable ``LINE_FLOW`` should be defined so that the library can be located in the directory ``$LINE_FLOW/lib``.

.. _start_docs:

Documentation
=============

Building this documentation locally requires `Sphinx <http://www.sphinx-doc.org/en/stable/>`_, `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_, and typing the command ``make html`` inside the ``docs`` directory. The environment variable ``PFNET_DOCS`` must be defined and should point to the output directory for the html files. The source code documentation will be placed in the directory ``PFNET_DOCS/c``.

.. _start_example:

Example
=======

Coming soon. 

.. _Graphviz: http://www.graphviz.org/
