.. _ext:

**********
Extensions
**********

This section describes how to add custom parsers, functions, and constraints to PFNET.

.. _ext_parser:

Adding a Parser
===============

To add a new data parser to PFNET, one should create a subclass of the :class:`CustomParser <pfnet.CustomParser>` class and provide five methods:

* :func:`init(self) <pfnet.ParserBase.init>`: Initializes parser data. This is called every time the :func:`parse() <pfnet.ParserBase.parse>` method is called.
* :func:`parse(self,filename,num_periods=1) <pfnet.ParserBase.parse>`: Reads power network data from ``filename`` and constructs and returns a :class:`Network <pfnet.Network>` with the given number of time periods. The custom parser only needs to allocate the :ref:`network components <net_components>` and store their information for time period zero. PFNET allocates flag and other utility arrays, copies the data to all time periods, and computes the initial network properties automatically. 
* :func:`set(self,key,value) <pfnet.ParserBase.set>`: Sets the parser's ``key`` parameter or option to ``value``, if it exists. 
* :func:`show(self) <pfnet.ParserBase.show>`: Shows information about the parser's internal data. 
* :func:`write(self,net,filename) <pfnet.ParserBase.write>`: Writes the given :class:`Network <pfnet.Network>` to ``filename`` using the format assoicated with the custom parser. 

.. _ext_func:

Adding a Function
=================



.. _ext_constr:

Adding a Constraint
===================
