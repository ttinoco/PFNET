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

.. _ref_bus:

Bus
===

.. _ref_bus_prop:

Bus Properties
--------------

.. data:: 'any'
.. data:: 'slack'
.. data:: 'regulated by generator'
.. data:: 'regulated by transformer'
.. data:: 'regulated by shunt'
.. data:: 'not slack'
.. data:: 'not regulated by generator'

.. _ref_bus_q:

Bus Quantities
--------------

.. data:: 'all'
.. data:: 'voltage angle'
.. data:: 'voltage magnitude'
.. data:: 'voltage magnitude deviation' 
.. data:: 'voltage magnitude violation'

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

.. data:: 'any'
.. data:: 'tap changer'
.. data:: 'tap changer - v' (controls voltage magnitude)
.. data:: 'tap changer - Q' (controls reactive flow)
.. data:: 'phase shifter'
.. data:: 'not on outage'

.. _ref_branch_q:

Branch Quantities
-----------------

.. data:: 'all'
.. data:: 'phase shift'
.. data:: 'tap ratio'
.. data:: 'tap ratio deviation'

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

.. data:: 'any'
.. data:: 'slack'
.. data:: 'regulator'
.. data:: 'not slack' 
.. data:: 'not regulator'
.. data:: 'not on outage'
.. data:: 'adjustable active power'

.. _ref_gen_q:

Generator Quantities
--------------------

.. data:: 'all'
.. data:: 'active power'
.. data:: 'reactive power' 

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

.. data:: 'any'
.. data:: 'switching - v' (controls voltage magnitude)

.. _ref_shunt_q:

Shunt Quantities
----------------

.. data:: 'all'
.. data:: 'susceptance'
.. data:: 'susceptance deviation'

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

.. data:: 'any'
.. data:: 'adjustable active power'

.. _ref_load_q:

Load Quantities
---------------

.. data:: 'all'
.. data:: 'active power'

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

.. data:: 'any' 

.. _ref_vargen_q:

Variable Generator Quantities
-----------------------------

.. data:: 'all' 
.. data:: 'active power'
.. data:: 'reactive power'

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

.. data:: 'any'

.. _ref_bat_q:

Battery Quantities
------------------

.. data:: 'all' 
.. data:: 'charging power'
.. data:: 'energy level'

.. _ref_bat_class:

Battery Class
-------------

.. autoclass:: pfnet.Battery
   :members:

.. _ref_net:

Network
=======

.. _ref_net_obj:

Component Types
---------------

.. data:: 'all'
.. data:: 'bus'
.. data:: 'generator'
.. data:: 'branch'
.. data:: 'shunt' 
.. data:: 'load'
.. data:: 'variable generator'
.. data:: 'battery'
.. data:: 'unknown' 

.. _ref_net_flag:

Flag Types
----------

.. data:: 'variable'

	  For selecting quantities to be variables.

.. data:: 'fixed'

	  For selecting variables to be fixed.

.. data:: 'bounded'

	  For selecting variables to be bounded.

.. data:: 'sparse'

	  For selecting control adjustments to be sparse.

.. _ref_var_values:

Variable Value Options
----------------------

.. data:: 'current'
.. data:: 'upper limits'
.. data:: 'lower limits'

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

.. _ref_func_type:

Function Types
--------------

.. data:: 'voltage magnitude regularization' 
.. data:: 'voltage angle regularization'
.. data:: 'generator powers regularization'
.. data:: 'tap ratio regularization'
.. data:: 'phase shift regularization'
.. data:: 'susceptance regularization'
.. data:: 'generation cost'
.. data:: 'sparse controls penalty'
.. data:: 'soft voltage magnitude limits'
.. data:: 'consumption utility' 
.. data:: 'net consumption cost'
 
.. _ref_func_class:

Function Class
--------------

.. autoclass:: pfnet.Function
   :members:

.. _ref_constr:

Constraint
==========

.. _ref_constr_type:

Constraint Types
----------------

.. data:: 'AC power balance'
.. data:: 'DC power balance'
.. data:: 'linearized AC power balance'
.. data:: 'variable fixing' 
.. data:: 'variable nonlinear bounds'
.. data:: 'generator active power participation'
.. data:: 'generator reactive power participation'
.. data:: 'voltage regulation by generators'
.. data:: 'voltage regulation by transformers'
.. data:: 'voltage regulation by shunts'
.. data:: 'DC branch flow limits'
.. data:: 'variable bounds'
.. data:: 'generator ramp limits'

.. _ref_constr_class:

Constraint Class
----------------

.. autoclass:: pfnet.Constraint
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

