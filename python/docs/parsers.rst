.. _parsers:

************
Data Parsers
************

This section describes the different data parsers available in PFNET and the supported file types. 

.. _parsers_mat:

MATPOWER case files
-------------------

`MATPOWER <http://www.pserc.cornell.edu//matpower/>`_ is a `MATLAB <http://www.mathworks.com/products/matlab/>`_ package for solving power flow and optimal power flow problems. It contains several power flow and optimal power flow cases defined in `MATLAB <http://www.mathworks.com/products/matlab/>`_ files. These "M" files can be converted to CSV files using the script :download:`mpc2mat.m <../../tools/mpc2mat.m>`. These MATPOWER-converted CSV files have extension ``.mat`` and can be used to load power networks in PFNET.

.. _parser_raw:

RAW case files
--------------

.. include:: <isonum.txt> 

If built with raw parsing capabilities, which requires linking PFNET with ``libraw_parser``, PFNET can load power networks from files with extension ``.raw``. These files are used by the software PSS |reg| E and are widely used by North American power system operators.
