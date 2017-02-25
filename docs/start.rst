.. _start:

***************
Getting Started
***************

.. _start_requirements:

Dependencies
============

PFNET for Python has the following dependencies:

* `Graphviz`_ (libgvc >= 2.38) (optional) 
* Raw parser (libraw_parser >= 1.2.1) (optional)

.. _start_download:

Download
========

The latest version of PFNET can be obtained from `<https://github.com/ttinoco/PFNET>`_.

.. _start_linux:

Build Instructions (Linux)
==========================

Building PFNET in Linux requires typing the command ``make`` in the root directory of the library, say ``$PFNET``. This creates the shared library ``libpfnet.so`` in the ``$PFNET/lib`` directory.

To build the library without visualization capabilities (no ``libgvc`` dependency), ``make`` should be passed the argument ``NO_GRAPHVIZ=1``. Otherwise, the compiler needs to find the header ``graphviz/gvc.h``, and the linker needs to find ``libgvc.so`` (to build the tests). If `Graphviz`_ has been installed in some non-standard directory, then ``$GRAPHVIZ`` should be set to that directory so that the compiler can find the above header in ``$GRAPHVIZ/include``. The directory ``$GRAPHVIZ/lib`` can be added to ``$LD_LIBRARY_PATH`` so that the linker can find the above shared library.

To build the library without raw parsing capabilities (no ``libraw_parser`` dependency), ``make`` should be passed the argument ``NO_RAW_PARSER=1``.

.. _start_linux_docs:

Documentation
-------------

Building this documentation for PFNET requires `Sphinx <http://www.sphinx-doc.org/en/stable/>`_ and `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_. The environment variable ``$PFNET_DOCS`` must be set to the location where the documentation will be moved to once it is built. After this environment variable has been set, run the command ``make docs`` from the roof directory of PFNET.

.. _start_mac:

Build Instructions (Mac OS X)
=============================

In Mac OS X, PFNET can be easily installed using `Homebrew <http://brew.sh>`_ with::

  brew install --HEAD https://raw.githubusercontent.com/ttinoco/PFNET/master/pfnet.rb

Since Homebrew does currently not support external dependencies, PFNET cannot be built with raw parsing capabilities with this method. Homebrew will install `Graphviz`_ as a dependency unless the option ``--without-graphviz`` is used. 

If building from source directly, the instructions for Linux should suffice but substitute ``$DYLD_FALLBACK_LIBRARY_PATH`` for ``$LD_LIBRARY_PATH``.

.. _Graphviz: http://www.graphviz.org/
