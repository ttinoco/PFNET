.. _parsers:

************
Data Parsers
************

This section describes the different data parsers available in PFNET and the supported file types. 

.. _parsers_mat:

MATPOWER case files
-------------------

`MATPOWER <http://www.pserc.cornell.edu//matpower/>`_ is a `MATLAB <http://www.mathworks.com/products/matlab/>`_ package for solving power flow and optimal power flow problems. It contains several power flow and optimal power flow cases defined in `MATLAB <http://www.mathworks.com/products/matlab/>`_ files. These "M" files can be converted to CSV files using the script :download:`mpc2mat.m <../../tools/mpc2mat.m>`. These MATPOWER-converted CSV files have extension ``.mat`` and can be used to load power networks in PFNET.

.. _parser_artere:

ARTERE case files
-----------------

PFNET can load networks from case files used by `ARTERE <http://www.montefiore.ulg.ac.be/~vct/software.html>`_, which is a software for performing power flow computations using the Newton-Raphson method. These files should have extension ``.art``. Details about these data files can be found in the document `"ARTERE: description of data files" <http://www.montefiore.ulg.ac.be/~vct/software/ARTERE_data.pdf>`_.

Currently, PFNET has limited support of these files. More specifically:

* Components with open breakers are ignored.
* For LTC-V devices, tap positions are treated as continuous and the optional fields are ignored. 
* The SWITCH, TRFO, PSHIFT-P, TURLIM, SVC, LFRESV, BUSPART and BRAPART records are not supported.
* Computation control parameters are ignored.

.. _parser_raw:

RAW case files
--------------

.. include:: <isonum.txt> 

If built with raw parsing capabilities, which requires linking PFNET with ``libraw_parser``, PFNET can load power networks from files with extension ``.raw``. These files are used by the software PSS |reg| E and are widely used by North American power system operators.
