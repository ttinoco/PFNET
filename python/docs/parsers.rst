.. include:: defs.hrst

.. _parsers:

************
Data Parsers
************

This section describes the different data parsers available in PFNET. 

.. _parsers_overview:

Overview
========

Parsers in PFNET are subclasses of the |ParserBase| class. They can be used to read a power network data file and create a |Network|. For convenience, a format-specific parser can be instantiated from the |Parser| class by specifying the file extension or a sample file name::

  >>> import pfnet
  >>> parser = pfnet.Parser("mat")
  >>> network = parser.parse("ieee14.mat")

For this and subsequent examples, is it assumed that the Python interpreter was started in a directory where the sample case |ieee14| can be found.

.. _parsers_json:

JSON Data Parser
================

PFNET networks can be constructed from data files in the widely-popular lightweight data-interchange format |JSON|. These network data files have extension ``.json`` and parsers for them can be instantiated from the class |ParserJSON|. These JSON parsers also allow writing a given network to a file, as the example below shows::

  >>> pfnet.ParserJSON().write(network,"new_network.json")

For creating, visualizing, or modifying these JSON network files, online editors such as the following may be used for convenience:

* `<http://jsoneditoronline.org>`_
* `<http://www.cleancss.com/json-editor>`_

In the top-level object of the JSON data, *i.e.*, the network, the field ``version`` indicates the PFNET version associated with the data.    
  
.. _parsers_mat:

MATPOWER Data Parser
=====================

|MATPOWER| is a popular |MATLAB| package for solving power flow and optimal power flow problems. It contains several power flow and optimal power flow cases defined in |MATLAB| files. These "m" files can be converted to CSV files using the script :download:`mpc2mat.m <../../tools/mpc2mat.m>`. These MATPOWER-converted CSV files have extension ``.mat`` and can be used to create power networks in PFNET. A parser for these data files can be constructed from the class |ParserMAT|.

.. _parser_artere:

ARTERE Data Parser
==================

PFNET can construct networks from data files used by |ARTERE|, which is a software for performing power flow computations using the Newton-Raphson method. These files should have extension ``.art``. Details about these data files can be found in the document `"ARTERE: description of data files" <http://www.montefiore.ulg.ac.be/~vct/software/ARTERE_data.pdf>`_. A parser for these data files can be constructed from the class |ParserART|.

Currently, there is limited support for these files. More specifically:

* Components with open breakers are ignored.
* For LTC-V devices, tap positions are treated as continuous and the optional fields are ignored. 
* The SWITCH, TRFO, PSHIFT-P, TURLIM, SVC, LFRESV, BUSPART and BRAPART records are not supported.
* Computation control parameters are ignored.

.. _parser_raw:

RAW Data Parser
===============

.. include:: <isonum.txt> 

If built with "raw" parsing capabilities, PFNET can construct power networks from files with extension ``.raw``. These files are used by the software PSS |reg| E, which is widely used by North American power system operators. A parser for these data files can be constructed from the class |ParserRAW|.
