.. _parsers:

************
Data Parsers
************

This section describes the different data parsers available in PFNET and the supported file types. 

.. _parsers_overview:

Overview
========

Parsers in PFNET are subclasses of the :class:`ParserBase <pfnet.ParserBase>` class. They can be used to read a power network data file and create a :class:`Power Network <pfnet.Network>`. For convenience, a format-specific parser can be instantiated from the :class:`Parser <pfnet.Parser>` class by specifying the file extension or by specifying a sample file name::

  >>> import pfnet
  >>> parser = pfnet.Parser("mat")
  >>> network = parser.parse("ieee14.mat")

In this example, is it assumed that the Python interpreter was started from the ``data`` directory of the PFNET library, where the sample case ``ieee14.mat`` is located.

.. _parsers_mat:

MATPOWER Data Files
===================

`MATPOWER`_ is a `MATLAB`_ package for solving power flow and optimal power flow problems. It contains several power flow and optimal power flow cases defined in `MATLAB`_ files. These "m" files can be converted to CSV files using the script :download:`mpc2mat.m <../../tools/mpc2mat.m>`. These MATPOWER-converted CSV files have extension ``.mat`` and can then be used to create power networks in PFNET. A parser for these data files can be constructed from the class :class:`ParserMAT <pfnet.ParserMAT>`.

.. _parser_artere:

ARTERE Data Files
=================

PFNET can construct networks from data files used by `ARTERE`_, which is a software for performing power flow computations using the Newton-Raphson method. These files should have extension ``.art``. Details about these data files can be found in the document `"ARTERE: description of data files" <http://www.montefiore.ulg.ac.be/~vct/software/ARTERE_data.pdf>`_. A parser for these data files can be constructed from the class :class:`ParserART <pfnet.ParserART>`.

Currently, there is limited support for these files. More specifically:

* Components with open breakers are ignored.
* For LTC-V devices, tap positions are treated as continuous and the optional fields are ignored. 
* The SWITCH, TRFO, PSHIFT-P, TURLIM, SVC, LFRESV, BUSPART and BRAPART records are not supported.
* Computation control parameters are ignored.

.. _parser_raw:

RAW Data Files
==============

.. include:: <isonum.txt> 

If built with raw parsing capabilities, PFNET can construct power networks from files with extension ``.raw``. These files are used by the software PSS |reg| E, which is widely used by North American power system operators. A parser for these data files can be constructed from the class :class:`ParserRAW <pfnet.ParserRAW>`.

.. _ARTERE: http://www.montefiore.ulg.ac.be/~vct/software.html
.. _MATPOWER: http://www.pserc.cornell.edu//matpower/
.. _MATLAB: http://www.mathworks.com/products/matlab/
