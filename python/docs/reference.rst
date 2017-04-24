.. _reference:

*************
API Reference
*************

.. _ref_vec:

Vector
======

.. class:: numpy.ndarray

           See `numpy documentation <http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html>`_.

.. _ref_mat:

Matrix
======

.. class:: scipy.sparse.coo_matrix

	   See `scipy documentation <http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.coo_matrix.html>`_.

.. _ref_parser:

Parser
======

.. autoclass:: pfnet.ParserBase
   :members: 

.. autoclass:: pfnet.Parser
.. autoclass:: pfnet.ParserMAT
.. autoclass:: pfnet.ParserART
.. autoclass:: pfnet.ParserRAW
.. autoclass:: pfnet.CustomParser

.. _ref_bus:

Bus
===

.. _ref_bus_prop:

Bus Properties
--------------

================================ ========
================================ ========
``"any"``
``"slack"``
``"regulated by generator"``
``"regulated by transformer"``
``"regulated by shunt"``
``"not slack"``
``"not regulated by generator"``
================================ ========

.. _ref_bus_q:

Bus Quantities
--------------

================================= ========
================================= ========
``"all"``
``"voltage angle"``
``"voltage magnitude"``
================================= ========

.. _ref_bus_sens:

Bus Sensitivities
-----------------

.. data:: pfnet.BUS_SENS_LARGEST

          Largest objective function sensitivity with respect to constraints involving this bus.
	  
.. data:: pfnet.BUS_SENS_P_BALANCE

	  Objective function sensitivity with respect to active power balance.

.. data:: pfnet.BUS_SENS_Q_BALANCE

	  Objective function sensitivity with respect to reactive power balance.

.. data:: pfnet.BUS_SENS_V_MAG_U_BOUND

	  Objective function sensitivity with respect to voltage magnitude upper bound.

.. data:: pfnet.BUS_SENS_V_MAG_L_BOUND

	  Objective function sensitivity with respect to voltage magnitude lower bound.

.. data:: pfnet.BUS_SENS_V_ANG_U_BOUND

	  Objective function sensitivity with respect to voltage angle upper bound.

.. data:: pfnet.BUS_SENS_V_ANG_L_BOUND

	  Objective function sensitivity with respect to voltage angle lower bound.

.. data:: pfnet.BUS_SENS_V_REG_BY_GEN

	  Objective function sensitivity with respect to voltage magnitude regulation by generators.

.. data:: pfnet.BUS_SENS_V_REG_BY_TRAN

	  Objective function sensitivity with respect to voltage magnitude regulation by tap-changing transformers.

.. data:: pfnet.BUS_SENS_V_REG_BY_SHUNT

	  Objective function sensitivity with respect to voltage magnitude regulation by switched shunt devices.

.. _ref_bus_mis:

Bus Power Mismatches
--------------------

.. data:: pfnet.BUS_MIS_LARGEST

	  Largest bus power mismatch.

.. data:: pfnet.BUS_MIS_ACTIVE

	  Bus active power mismatch.

.. data:: pfnet.BUS_MIS_REACTIVE

	  Bus reactive power mismatch.

.. _ref_bus_class:

Bus Class
---------

.. autoclass:: pfnet.Bus
   :members:

.. _ref_branch:

Branch
======

.. _ref_branch_prop:

Branch Properties
-----------------

===================== ============================
===================== ============================
``"any"``
``"tap changer"``
``"tap changer - v"`` Controls voltage magnitude
``"tap changer - Q"`` Controls reactive flow
``"phase shifter"``
``"not on outage"``
===================== ============================

.. _ref_branch_q:

Branch Quantities
-----------------

========================= =======
========================= =======
``"all"``
``"phase shift"``
``"tap ratio"``
========================= =======

.. _ref_branch_class:

Branch Class
------------

.. autoclass:: pfnet.Branch
   :members:

.. _ref_gen:

Generator
=========

.. _ref_gen_prop:

Generator Properties
--------------------

============================= ===========================
============================= ===========================
``"any"``
``"slack"``
``"regulator"``
``"not slack"``
``"not regulator"``
``"not on outage"``
``"adjustable active power"`` :math:`P_{\min} < P_{\max}`
============================= ===========================

.. _ref_gen_q:

Generator Quantities
--------------------

==================== =======
==================== =======
``"all"``
``"active power"``
``"reactive power"``
==================== =======

.. _ref_gen_class:

Generator Class
---------------

.. autoclass:: pfnet.Generator
   :members:

.. _ref_shunt:

Shunt
=====

.. _ref_shunt_prop:

Shunt Properties
----------------

=================== ============================
=================== ============================
``"any"``
``"switching - v"`` Controls voltage magnitude
=================== ============================	
  
.. _ref_shunt_q:

Shunt Quantities
----------------

=========================== =======
=========================== =======
``"all"``
``"susceptance"``
=========================== =======

.. _ref_shunt_class:

Shunt Class
-----------

.. autoclass:: pfnet.Shunt
   :members:

.. _ref_load:

Load
====

.. _ref_load_prop:

Load Properties
---------------

============================= ===========================
============================= ===========================
``"any"``
``"adjustable active power"`` :math:`P_{\min} < P_{\max}`
============================= ===========================

.. _ref_load_q:

Load Quantities
---------------

==================== =======
==================== =======
``"all"``
``"active power"``
``"reactive power"``
==================== =======

.. _ref_load_class:

Load Class
----------

.. autoclass:: pfnet.Load
   :members:

.. _ref_vargen:

Variable Generator
==================

.. _ref_vargen_prop:

Variable Generator Properties
-----------------------------

========= =======
========= =======
``"any"``
========= =======

.. _ref_vargen_q:

Variable Generator Quantities
-----------------------------

==================== =======
==================== =======
``"all"``
``"active power"``
``"reactive power"``
==================== =======

.. _ref_vargen_class:

Variable Generator Class
------------------------

.. autoclass:: pfnet.VarGenerator
   :members:

.. _ref_bat:

Battery
=======

.. _ref_bat_prop:

Battery Properties
------------------

========= =======
========= =======
``"any"``
========= =======

.. _ref_bat_q:

Battery Quantities
------------------

==================== =======
==================== =======
``"all"``
``"charging power"``
``"energy level"``
==================== =======

.. _ref_bat_class:

Battery Class
-------------

.. autoclass:: pfnet.Battery
   :members:

.. _ref_net:

Network
=======

.. _ref_net_obj:

Component Names
---------------

======================== =======
======================== =======
``"all"``
``"bus"``
``"battery"``
``"branch"``
``"generator"``
``"load"``
``"shunt"``
``"variable generator"``
======================== =======

.. _ref_net_flag:

Flags
-----

============== ==============================================
============== ==============================================
``"variable"`` For selecting quantities to be variables
``"fixed"``    For selecting variables to be fixed
``"bounded"``  For selecting variables to be bounded.
``"sparse"``   For selecting control adjustments to be sparse
============== ==============================================

.. _ref_var_values:

Variable Value Options
----------------------

================== =======
================== =======
``"current"``
``"upper limits"``
``"lower limits"``
================== =======

.. _ref_net_class:

Network Class
-------------

.. autoclass:: pfnet.Network
   :members:

.. _ref_cont:

Contingency
===========

.. autoclass:: pfnet.Contingency
   :members:

.. _ref_graph:

Graph
=====

.. autoclass:: pfnet.Graph
   :members:

.. _ref_func:

Function
========

.. _ref_func_names:

Function Names
--------------

====================================== =======
====================================== =======
``"consumption utility"`` 
``"generation cost"``
``"generator powers regularization"``
``"net consumption cost"``
``"phase shift regularization"``
``"susceptance regularization"``
``"tap ratio regularization"``
``"soft voltage magnitude limits"``
``"sparse controls penalty"``
``"voltage magnitude regularization"`` 
``"voltage angle regularization"``
====================================== =======
 
.. _ref_func_class:

Function Classes
----------------

.. autoclass:: pfnet.FunctionBase
   :members:

.. autoclass:: pfnet.Function
.. autoclass:: pfnet.CustomFunction
   :members:

.. _ref_constr:

Constraint
==========

.. _ref_constr_names:

Constraint Names
----------------

============================================ =======
============================================ =======
``"AC power balance"``
``"DC power balance"``
``"linearized AC power balance"``
``"variable fixing"``
``"variable bounds"``
``"variable nonlinear bounds"``
``"generator active power participation"``
``"generator reactive power participation"``
``"generator ramp limits"``
``"voltage regulation by generators"``
``"voltage regulation by transformers"``
``"voltage regulation by shunts"``
``"AC branch flow limits"``
``"DC branch flow limits"``
``"battery dynamics"``
============================================ =======

.. _ref_constr_class:

Constraint Classes
------------------

.. autoclass:: pfnet.ConstraintBase
   :members:

.. autoclass:: pfnet.Constraint
.. autoclass:: pfnet.CustomConstraint
   :members:

.. _ref_problem:

Optimization Problem
====================

.. _ref_problem_class:

Problem Class
-------------

.. autoclass:: pfnet.Problem
   :members:
   :exclude-members: add_heuristic, apply_heuristics

.. _ref_references:

References
==========

.. [TTR2015] T\. Tinoco De Rubira, *Numerical Optimization and Modeling Techniques for Power System Operations and Planning*. PhD thesis, Stanford University, March 2015.

