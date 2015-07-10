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

Bus Property Masks
------------------

.. data:: pfnet.BUS_PROP_ANY 
 
	  Any bus.
           
.. data:: pfnet.BUS_PROP_SLACK

	  Slack bus.
 
.. data:: pfnet.BUS_PROP_REG_BY_GEN

          Bus voltage magnitude is regulated by generators.

.. data:: pfnet.BUS_PROP_REG_BY_TRAN

	  Bus voltage magnitude is regulated by tap-changing transformers.

.. data:: pfnet.BUS_PROP_REG_BY_SHUNT

	  Bus voltage magnitude is regulated by switched shunt devices.
	  
.. data:: pfnet.BUS_PROP_NOT_REG_BY_GEN

	  Bus voltage magnitude is not regulated by generators.

.. data:: pfnet.BUS_PROP_NOT_SLACK

	  Bus is not slack. 

.. _ref_bus_var:

Bus Variable Masks
------------------

.. data:: pfnet.BUS_VAR_VMAG

          Voltage magnitude.

.. data:: pfnet.BUS_VAR_VANG

          Voltage angle.

.. data:: pfnet.BUS_VAR_VDEV

          Voltage magnitude positive and negative set point deviations.

.. data:: pfnet.BUS_VAR_VVIO

          Voltage magnitude upper and lower bound violations.

.. _ref_bus_sens:

Bus Sensitivities
-----------------

.. data:: pfnet.BUS_SENS_LARGEST

          Largest objective function sensitivity with respect to nonlinear equality constraints involving this bus.
	  
.. data:: pfnet.BUS_SENS_P_BALANCE

	  Objective function sensitivity with respect to bus active power balance.

.. data:: pfnet.BUS_SENS_Q_BALANCE

	  Objective function sensitivity with respect to bus reactive power balance.

.. data:: pfnet.BUS_SENS_V_MAG_U_BOUND

	  Objective function sensitivity with respect to bus upper voltage bound.

.. data:: pfnet.BUS_SENS_V_MAG_L_BOUND

	  Objective function sensitivity with respect to bus lower voltage bound.

.. data:: pfnet.BUS_SENS_V_REG_BY_GEN

	  Objective function sensitivity with respect to bus voltage magnitude regulation by generators.

.. data:: pfnet.BUS_SENS_V_REG_BY_TRAN

	  Objective function sensitivity with respect to bus voltage magnitude regulation by tap-changing transformers.

.. data:: pfnet.BUS_SENS_V_REG_BY_SHUNT

	  Objective function sensitivity with respect to bus voltage magnitude regulation by switched shunt devices.

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

Branch Property Masks
---------------------

.. data:: pfnet.BRANCH_PROP_ANY

	  Any branch.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER

	  Branch is tap-changing transformer.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER_V

	  Branch is tap-changing transformer regulating bus voltage magnitude.

.. data:: pfnet.BRANCH_PROP_TAP_CHANGER_Q 

	  Branch is tap-changing transformer regulating reactive power flow.

.. data:: pfnet.BRANCH_PROP_PHASE_SHIFTER

	  Branch is phase-shifting transformer regulating active power flow.

.. _ref_branch_var:

Branch Variable Masks
---------------------

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

Generator Property Masks
------------------------

.. data:: pfnet.GEN_PROP_ANY

	  Any generator.

.. data:: pfnet.GEN_PROP_SLACK

	  Slack generator.

.. data:: pfnet.GEN_PROP_REG

	  Generator that regulates bus voltage magnitude.

.. data:: pfnet.GEN_PROP_NOT_REG

	  Generator that does not regulate bus voltage magnitude.

.. data:: pfnet.GEN_PROP_NOT_SLACK

	  Generator that is not slack.

.. _ref_gen_var:

Generator Variable Masks
------------------------

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

Shunt Property Masks
--------------------

.. data:: pfnet.SHUNT_PROP_ANY

	  Any shunt.

.. data:: pfnet.SHUNT_PROP_SWITCHED_V

	  Switched shunt devices that regulates bus voltage magnitude.

.. _ref_shunt_var:

Shunt Variable Masks
--------------------

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

.. _ref_load_class:

Load Class
----------

.. autoclass:: pfnet.Load
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

.. _ref_net_flag:

Flag Masks
----------

.. data:: pfnet.FLAG_VARS

	  For specifying quantities as variable.

.. data:: pfnet.FLAG_FIXED

	  For specifying variables that should be fixed.

.. data:: pfnet.FLAG_BOUNDED

	  For specifying variables that should be bounded.

.. data:: pfnet.FLAG_SPARSE

	  For specifying control adjustments that should be sparse.

.. _ref_net_class:

Network Class
-------------

.. autoclass:: pfnet.Network
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

	  Constraint for enforcing power balance at every bus of the network. 

.. data:: pfnet.CONSTR_TYPE_FIX

	  Constraint for fixing a subset of variables to their current value.

.. data:: pfnet.CONSTR_TYPE_BOUND

	  Constraint for forcing a subset of variables to be within their bounds.

.. data:: pfnet.CONSTR_TYPE_PAR_GEN

	  Constraint for enforcing generator participations.

.. data:: pfnet.CONSTR_TYPE_REG_GEN

	  Constraint for enforcing voltage set point regulation by generators.

.. data:: pfnet.CONSTR_TYPE_REG_TRAN
	  
	  Constraint for enforcing voltage band regulation by tap-changing transformers.

.. data:: pfnet.CONSTR_TYPE_REG_SHUNT

	  Constraint for enforcing voltage band regulation by switched shunt devices.

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

