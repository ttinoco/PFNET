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

.. data:: pfnet.BUS_PROP_ANY 
 
	  Any bus.
           
.. data:: pfnet.BUS_PROP_SLACK

	  Slack bus.
 
.. data:: pfnet.BUS_PROP_REG_BY_GEN

          Bus with voltage magnitude regulated by one or more generators.

.. data:: pfnet.BUS_PROP_REG_BY_TRAN

	  Bus with voltage magnitude regulated by one or more tap-changing transformers.

.. data:: pfnet.BUS_PROP_REG_BY_SHUNT

	  Bus with voltage magnitude regulated by one or more switched shunt devices.
	  
.. data:: pfnet.BUS_PROP_NOT_REG_BY_GEN

	  Bus with voltage magnitude that is not regulated by generators.

.. data:: pfnet.BUS_PROP_NOT_SLACK

	  Bus that is not a slack bus. 

.. _ref_bus_var:

Bus Variables
-------------

.. data:: pfnet.BUS_VAR_VMAG

          Bus voltage magnitude.

.. data:: pfnet.BUS_VAR_VANG

          Bus voltage angle.

.. data:: pfnet.BUS_VAR_VDEV

          Bus voltage magnitude positive and negative set-point deviations.

.. data:: pfnet.BUS_VAR_VVIO

          Bus voltage magnitude upper and lower bound violations.

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

.. data:: pfnet.BRANCH_PROP_ANY

	  Any branch.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER

	  Branch that is tap-changing transformer.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER_V

	  Branch that is tap-changing transformer regulating a bus voltage magnitude.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER_Q 

	  Branch that is tap-changing transformer regulating reactive power flow.

.. data:: pfnet.BRANCH_PROP_PHASE_SHIFTER

	  Branch that is phase-shifting transformer regulating active power flow.

.. data:: pfnet.BRANCH_PROP_NOT_OUT

	  Branch that is not on outage.

.. _ref_branch_var:

Branch Variables
----------------

.. data:: pfnet.BRANCH_VAR_RATIO

	  Transformer tap ratio.

.. data:: pfnet.BRANCH_VAR_RATIO_DEV

	  Transformer tap ratio deviations from current value.

.. data:: pfnet.BRANCH_VAR_PHASE

	  Transformer phase shift.

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

.. data:: pfnet.GEN_PROP_ANY

	  Any generator.

.. data:: pfnet.GEN_PROP_SLACK

	  Slack generator.

.. data:: pfnet.GEN_PROP_REG

	  Generator that regulates a bus voltage magnitude.

.. data:: pfnet.GEN_PROP_NOT_REG

	  Generator that does not regulate a bus voltage magnitude.

.. data:: pfnet.GEN_PROP_NOT_SLACK

	  Generator that is not a slack generator.

.. data:: pfnet.GEN_PROP_NOT_OUT

	  Generator that is not on outage.

.. data:: pfnet.GEN_PROP_P_ADJUST

	  Generator that can adjust its active power, e.g., :math:`P_{\min} < P_{\max}`.

.. _ref_gen_var:

Generator Variables
-------------------

.. data:: pfnet.GEN_VAR_P

	  Generator active power.

.. data:: pfnet.GEN_VAR_Q

	  Generator reactive power.

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

.. data:: pfnet.SHUNT_PROP_ANY

	  Any shunt.

.. data:: pfnet.SHUNT_PROP_SWITCHED_V

	  Switched shunt devices that regulates a bus voltage magnitude.

.. _ref_shunt_var:

Shunt Variables
---------------

.. data:: pfnet.SHUNT_VAR_SUSC

	  Switched shunt susceptance.

.. data:: pfnet.SHUNT_VAR_SUSC_DEV

	  Switched shunt susceptance deviations from current point.

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

.. data:: pfnet.LOAD_PROP_ANY

	  Any load.

.. data:: pfnet.LOAD_PROP_P_ADJUST

	  Load that can adjust its active power, e.g., :math:`P_{\min} < P_{\max}`.

.. _ref_load_var:

Load Variables
--------------

.. data:: pfnet.LOAD_VAR_P

	  Load active power.

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

.. data:: pfnet.VARGEN_PROP_ANY

	  Any variable generator.

.. _ref_vargen_var:

Variable Generator Variables
----------------------------

.. data:: pfnet.VARGEN_VAR_P

	  Variable generator active power.

.. data:: pfnet.VARGEN_VAR_Q

	  Variable generator reactive power.

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

.. data:: pfnet.BAT_PROP_ANY

	  Any battery.

.. _ref_bat_var:

Battery Variables
-----------------

.. data:: pfnet.BAT_VAR_P

	  Battery charging/discharging power.

.. data:: pfnet.BAT_VAR_E

	  Battery energy level.

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

.. data:: pfnet.OBJ_BUS

	  Bus.

.. data:: pfnet.OBJ_GEN

	  Generator.

.. data:: pfnet.OBJ_BRANCH

	  Branch.

.. data:: pfnet.OBJ_SHUNT

	  Shunt device.

.. data:: pfnet.OBJ_LOAD

	  Load.

.. data:: pfnet.OBJ_VARGEN

	  Variable generator (solar, wind, etc).

.. data:: pfnet.OBJ_BAT

	  Battery.

.. data:: pfnet.OBJ_UNKNOWN

	  Unknown network component.

.. _ref_net_flag:

Flag Types
----------

.. data:: pfnet.FLAG_VARS

	  For specifying quantities as variables.

.. data:: pfnet.FLAG_FIXED

	  For specifying variables that should be fixed.

.. data:: pfnet.FLAG_BOUNDED

	  For specifying variables that should be bounded.

.. data:: pfnet.FLAG_SPARSE

	  For specifying control adjustments that should be sparse.

.. _ref_var_values:

Variable Value Codes
--------------------

.. data:: pfnet.CURRENT

	  Current variable value.

.. data:: pfnet.UPPER_LIMIT

	  Upper limit of variable.

.. data:: pfnet.LOWER_LIMIT

	  Lower limit of variable.

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

.. data:: pfnet.FUNC_TYPE_UNKNOWN

	  Unknown function.

.. data:: pfnet.FUNC_TYPE_REG_VMAG
	  
	  Bus voltage magnitude regularization.

.. data:: pfnet.FUNC_TYPE_SLIM_VMAG
	  
	  Bus voltage magnitude soft limits penalty.

.. data:: pfnet.FUNC_TYPE_REG_VANG

	  Bus voltage angle regularization.

.. data:: pfnet.FUNC_TYPE_REG_PQ

	  Generator active and reactive power regularization.

.. data:: pfnet.FUNC_TYPE_GEN_COST

          Active power generation cost.

.. data:: pfnet.FUNC_TYPE_NETCON_COST

          Net power consumption cost.

.. data:: pfnet.FUNC_TYPE_LOAD_UTIL

          Active power consumption utility.

.. data:: pfnet.FUNC_TYPE_REG_RATIO

	  Transformer tap ratio regularization.

.. data:: pfnet.FUNC_TYPE_REG_PHASE

	  Transformer phase shift regularization.

.. data:: pfnet.FUNC_TYPE_REG_SUSC

	  Switched shunt susceptance regularization.

.. data:: pfnet.FUNC_TYPE_SP_CONTROLS

	  Sparsity-inducing penalty for control adjustments.

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

.. data:: pfnet.CONSTR_TYPE_PF

	  Constraint for enforcing AC power balance at every bus of the network. 

.. data:: pfnet.CONSTR_TYPE_DCPF

	  Constraint for enforcing DC power balance at every bus of the network. 

.. data:: pfnet.CONSTR_TYPE_LINPF

	  Constraint for enforcing linearized power balance at every bus of the network. 

.. data:: pfnet.CONSTR_TYPE_FIX

	  Constraint for fixing a subset of variables to their current value.

.. data:: pfnet.CONSTR_TYPE_BOUND

	  Constraint for forcing a subset of variables to be within their bounds (nonlinear).

.. data:: pfnet.CONSTR_TYPE_LBOUND

	  Constraint for forcing a subset of variables to be within their bounds (linear).

.. data:: pfnet.CONSTR_TYPE_PAR_GEN_P

	  Constraint for enforcing generator active power participations.

.. data:: pfnet.CONSTR_TYPE_PAR_GEN_Q

	  Constraint for enforcing generator reactive power participations.

.. data:: pfnet.CONSTR_TYPE_REG_GEN

	  Constraint for enforcing voltage set point regulation by generators.

.. data:: pfnet.CONSTR_TYPE_REG_TRAN
	  
	  Constraint for enforcing voltage band regulation by tap-changing transformers.

.. data:: pfnet.CONSTR_TYPE_REG_SHUNT

	  Constraint for enforcing voltage band regulation by switched shunt devices.

.. data:: pfnet.CONSTR_TYPE_DC_FLOW_LIM

	  Constraint for enforcing DC power flow limits on every branch

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

