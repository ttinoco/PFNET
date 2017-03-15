.. _start:

***************
Getting Started
***************

.. _start_requirements:

Dependencies
============

The PFNET C library has the following optional dependencies:

* `Graphviz`_ (libgvc >= 2.38)
* Raw parser (libraw_parser >= 1.2.2)

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

If raw parser is available, the environment variable ``RAW_PARSER`` should point to the root directory of that library. 

.. _start_docs:

Documentation
=============

Building this documentation locally requires `Sphinx <http://www.sphinx-doc.org/en/stable/>`_, `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_, and typing the command ``make html`` inside the ``docs`` directory. The environment variable ``PFNET_DOCS`` must be defined and should point to the output directory for the html files. The source code documentation will be placed in the directory ``PFNET_DOCS/c``.

.. _start_example:

Example
=======

Coming soon. 

.. _Graphviz: http://www.graphviz.org/
