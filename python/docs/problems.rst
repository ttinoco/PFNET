.. _prob:

*********************
Optimization Problems
*********************

This section describes how to formulate power network optimization problems using PFNET.

.. _prob_func:

Objective Function
==================

The objective function :math:`\phi` for a network optimization problem created using PFNET is of the form

.. math::

   \varphi(x) = \sum_i w_i \varphi_i (x),

where :math:`w_i` are weights, :math:`\varphi_i` are general linear or nonlinear functions, and :math:`x` is a vector of values of network quantities that have been set as variables. Each weight-function pair in the summation is represented by an object of type :class:`Function <pfnet.Function>`. To instantiate an object of this type, the function type and weight need to be specified as well as the :class:`Network <pfnet.Network>` object that is to be associated with the function. The following example sets all bus voltage magnitudes as variables and constructs a function that penalizes voltage magnitude deviations from ideal values::

  >>> import pfnet as pf
  
  >>> net = pf.Network()
  >>> net.load('ieee14.mat')

  >>> net.set_flags(pf.OBJ_BUS,
  ...               pf.FLAG_VARS,
  ...               pf.BUS_PROP_ANY,
  ...               pf.BUS_VAR_VMAG)

  >>> func = pf.Function(pf.FUNC_TYPE_REG_VMAG,0.3,net)

  >>> print func.type == pf.FUNC_TYPE_REG_VMAG
  True

  >>> print func.weight
  0.3

After a :class:`Function <pfnet.Function>` object is created, its value, gradient and Hessian are zero, an empty vector, and an empty matrix, respectively. Before evaluating the function at a specific vector of values, it must be analyzed using the :class:`Function <pfnet.Function>` class method :func:`analyze() <pfnet.Function.analyze>`. This routine analyzes the function and allocated the required vectors and matrices for storing its gradient and Hessian. After this, the function can be evaluated using the method :func:`eval() <pfnet.Function.eval>`::

  >>> x = net.get_var_values()

  >>> func.analyze()

  >>> func.eval(x + 0.01)
  >>> func.eval(x)

The value :math:`\varphi_i(x)`, gradient :math:`\nabla \varphi_i(x)` and Hessian :math:`\nabla^2 \varphi_i(x)` of a function can then be extracted from the :data:`phi <pfnet.Function.phi>`, :data:`gphi <pfnet.Function.gphi>` and :data:`Hphi <pfnet.Function.Hphi>` attributes, respectively::

  >>> print x.shape
  (14,)

  >>> print func.phi
  0.255

  >>> print type(func.gphi), func.gphi.shape
  <type 'numpy.ndarray'> (14,)

  >>> print type(func.Hphi), func.Hphi.shape
  <class 'scipy.sparse.coo.coo_matrix'> (14, 14)

For the Hessian matrix, only the lower triangular part is stored.

Details about each of the different function types available in PFNET are provided below.

.. _prob_func_REG_VMAG:

Voltage magnitude regularization
--------------------------------

This function is of type :data:`FUNC_TYPE_REG_VMAG <pfnet.FUNC_TYPE_REG_VMAG>`. It penalizes deviations of bus voltage magnitudes from ideal values. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{v_k - v^t_k}{\Delta v} \Bigg)^2 + 
                 \frac{1}{2} \sum_k \Bigg( \frac{v^y_k}{\Delta v} \Bigg)^2 +
	         \frac{1}{2} \sum_k \Bigg( \frac{v^z_k}{\Delta v} \Bigg)^2 + 
                 \frac{1}{2} \sum_k \Bigg( \frac{v^h_k}{\Delta v} \Bigg)^2 +
	         \frac{1}{2} \sum_k \Bigg( \frac{v^l_k}{\Delta v} \Bigg)^2,

where :math:`v` are bus voltage magnitudes, :math:`v^t` are voltage magnitude set points (one for buses not regulated by generators), :math:`v^y` and :math:`v^z` are positive and negative deviations of :math:`v` from :math:`v^t`, :math:`v^h` and :math:`v^l` are voltage band upper and lower limit violations, and :math:`\Delta v` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_SLIM_VMAG:

Voltage magnitude soft limit penalty
------------------------------------

This function is of type :data:`FUNC_TYPE_SLIM_VMAG <pfnet.FUNC_TYPE_SLIM_VMAG>`. It reduces voltage (soft) limit violations by penalizing deviations of bus voltage magnitudes from the mid point of their ranges. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{v_k - \bar{v}_k}{\Delta v} \Bigg)^2,

where :math:`v` are bus voltage magnitudes, :math:`\bar{v}` are the mid points of their ranges, and :math:`\Delta v` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_ANG:

Voltage angle regularization
----------------------------

This function is of type :data:`FUNC_TYPE_REG_VANG <pfnet.FUNC_TYPE_REG_VANG>`. It penalizes large bus voltage angles and voltage angle differences across branches. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{\theta_k}{\Delta \theta} \Bigg)^2 + 
                 \frac{1}{2} \sum_{(k,m)} \Bigg( \frac{\theta_k - \theta_m - \phi_{km}}{\Delta \theta} \Bigg)^2,

where :math:`\theta` are bus voltage angles, :math:`\phi` are branch phase shifts, and :math:`\Delta \theta` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_PQ:

Generator powers regularization
-------------------------------

This function is of type :data:`FUNC_TYPE_REG_PQ <pfnet.FUNC_TYPE_REG_PQ>`. It penalizes deviations of generator powers from the midpoint of their ranges. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{P^g_k - \bar{P}_k}{\Delta P} \Bigg)^2 + 
                 \frac{1}{2} \sum_k \Bigg( \frac{Q^g_k - \bar{Q}_k}{\Delta Q} \Bigg)^2,

where :math:`P^g` and :math:`Q^g` are generator active and reactive powers, :math:`\bar{P}` and :math:`\bar{Q}` are midpoints of generator active and reactive power ranges, and :math:`\Delta P = \Delta Q` are normalization factors. Only terms that include optimization variables are included in the summation.

.. _prob_func_GEN_COST:

Active power generation cost
----------------------------

This function is of type :data:`FUNC_TYPE_GEN_COST <pfnet.FUNC_TYPE_GEN_COST>`. It measures active power generation cost by the expression

.. math::

   \varphi(x) := \sum_k q_{k0} + q_{k1} P_k + q_{k2} P_k^2,

where :math:`P_k` are generator active powers in per unit base system power, and :math:`q_{k0}`, :math:`q_{k1}`, and :math:`q_{k2}` are constant coefficients. These coefficients correspond to the attributes :data:`cost_coeff_Q0 <pfnet.Generator.cost_coeff_Q0>`, :data:`cost_coeff_Q1 <pfnet.Generator.cost_coeff_Q1>` and :data:`cost_coeff_Q2 <pfnet.Generator.cost_coeff_Q2>` of each :class:`Generator <pfnet.Generator>` object. 

.. _prob_func_LOAD_UTIL:

Active power consumption utility
--------------------------------

This function is of type :data:`FUNC_TYPE_LOAD_UTIL <pfnet.FUNC_TYPE_LOAD_UTIL>`. It measures active power consumption utility by the expression

.. math::

   \varphi(x) := \sum_k q_{k0} + q_{k1} P_k + q_{k2} P_k^2,

where :math:`P_k` are load active powers in per unit base system power, and :math:`q_{k0}`, :math:`q_{k1}`, and :math:`q_{k2}` are constant coefficients. These coefficients correspond to the attributes :data:`util_coeff_Q0 <pfnet.Load.util_coeff_Q0>`, :data:`util_coeff_Q1 <pfnet.Load.util_coeff_Q1>` and :data:`util_coeff_Q2 <pfnet.Load.util_coeff_Q2>` of each :class:`Load <pfnet.Load>` object. 

.. _prob_func_REG_RATIO:

Transformer tap ratio regularization
------------------------------------

This function is of type :data:`FUNC_TYPE_REG_RATIO <pfnet.FUNC_TYPE_REG_RATIO>`. It penalizes deviations of tap ratios of tap-changing transformers from their initial value. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{t_k - t^0_k}{\Delta t} \Bigg)^2 + 
                 \frac{1}{2} \sum_k \Bigg( \frac{t^y_k}{\Delta t} \Bigg)^2 + 
	         \frac{1}{2} \sum_k \Bigg( \frac{t^z_k}{\Delta t} \Bigg)^2,

where :math:`t` are tap ratios of tap-changing transformers, :math:`t^0` are their initial values, :math:`t^y` and :math:`t^z` are positive and negative deviations of :math:`t` from :math:`t^0`, and :math:`\Delta t` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_PHASE:

Transformer phase shift regularization
--------------------------------------

This function is of type :data:`FUNC_TYPE_REG_PHASE <pfnet.FUNC_TYPE_REG_PHASE>`. It penalizes deviations of phase shifts of phase shifting transformers from their initial value. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{\phi_k - \phi^0_k}{\Delta \phi} \Bigg)^2

where :math:`\phi` are phase shifts of phase-shifting transformers, :math:`\phi^0` are their initial values, and :math:`\Delta \phi` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_SUSC:

Switched shunt susceptance regularization
-----------------------------------------

This function is of type :data:`FUNC_TYPE_REG_SUSC <pfnet.FUNC_TYPE_REG_SUSC>`. It penalizes deviations of susceptances of switched shunt devices from their initial value. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_k \Bigg( \frac{b_k - b^0_k}{\Delta b} \Bigg)^2 + 
                 \frac{1}{2} \sum_k \Bigg( \frac{b^y_k}{\Delta b} \Bigg)^2 + 
	         \frac{1}{2} \sum_k \Bigg( \frac{b^z_k}{\Delta b} \Bigg)^2,

where :math:`b` are susceptances of switched shunt devices, :math:`b^0` are their initial values, :math:`b^y` and :math:`b^z` are positive and negative deviations of :math:`b` from :math:`b^0`, and :math:`\Delta b` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_SP_CONTROLS:

Sparsity inducing penalty for controls
--------------------------------------

This function is of type :data:`FUNC_TYPE_SP_CONTROLS <pfnet.FUNC_TYPE_SP_CONTROLS>`. It encourages sparse control adjustments with the expression

.. math::

   \varphi(x) := \sum_k \sqrt{ \Bigg( \frac{u_k - u_k^0}{\Delta u_k} \Bigg)^2 + \epsilon },

where :math:`u` are control quantities, :math:`u^0` are their current values, and :math:`\epsilon` is a small positive scalar. The normalization factors :math:`\Delta u_k` are given by

.. math::

   \Delta u_k := \max\{u^{\max}_k-u^{\min}_k, \delta\},

where :math:`u^{\max}` and :math:`u^{\min}` are control limits, and :math:`\delta` is a small positive scalar. The control quantities that are considered by this function are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` using the flag type :data:`FLAG_SPARSE <pfnet.FLAG_SPARSE>`.

.. _prob_constr:

Constraints
===========

Constraints in PFNET are of the form

.. math::
   
   & A x = b \\
   & f(x) = 0 \\
   & l \le G x \le u,

where :math:`A` and :math:`G`  are sparse matrices, :math:`b`, :math:`l` and :math:`u`  are vectors, :math:`f` is a vector-valued nonlinear function, and :math:`x` is a vector of values of network quantities that have been set as variables. They are represented by objects of type :class:`Constraint <pfnet.Constraint>`. To create an object of this type, the constraint type and the network to be associated with the constraint need to be specified. The following example sets all bus voltage magnitudes and angles as variables and constructs the power flow constraints::

  >>> import pfnet as pf

  >>> net = pf.Network()
  >>> net.load('ieee14.mat')

  >>> net.set_flags(pf.OBJ_BUS,
  ...               pf.FLAG_VARS,
  ...               pf.BUS_PROP_ANY,
  ...               pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

  >>> print net.num_vars == 2*net.num_buses
  True

  >>> constr = pf.Constraint(pf.CONSTR_TYPE_PF,net)

  >>> print constr.type == pf.CONSTR_TYPE_PF
  True

Before a :class:`Constraint <pfnet.Constraint>` object can be used, it must be initialized using the :class:`Constraint <pfnet.Constraint>` class method :func:`analyze() <pfnet.Constraint.analyze>`. This routine analyzes the constraint and allocates the required vectors and matrices. After this, the constraint can be evaluated using the method :func:`eval() <pfnet.Constraint.eval>`::

  >>> x = net.get_var_values()

  >>> constr.analyze()

  >>> constr.eval(x + 0.01)
  >>> constr.eval(x)

The matrices and vectors associated with the linear constraints can be extracted from the :data:`A <pfnet.Constraint.A>`, :data:`G <pfnet.Constraint.G>`, :data:`b <pfnet.Constraint.b>`, :data:`l <pfnet.Constraint.l>` and :data:`u <pfnet.Constraint.u>` attributes of the :class:`Constraint <pfnet.Constraint>` object. The vector of violations and Jacobian matrix of the nonlinear constraints can be extracted from the attributes :data:`f <pfnet.Constraint.f>` and :data:`J <pfnet.Constraint.J>`, respectively. Also, the Hessian matrix of any individual nonlinear constraint :math:`f_i(x) = 0` can be extracted using the class method :func:`get_H_single() <pfnet.Constraint.get_H_single>`. The following example shows how to extract the largest power flow mismatch in per unit :data:`system base power <pfnet.Network.base_power>` and the Hessian matrix corresponding to the active power balance constraint of a bus::

  >>> import numpy as np

  >>> f = constr.f

  >>> print type(f), f.shape
  <type 'numpy.ndarray'>  (28,)

  >>> print np.linalg.norm(f,np.inf)
  0.042
  
  >>> bus = net.get_bus(5)
  >>> Hi = constr.get_H_single(bus.index_P)

  >>> print type(Hi), Hi.shape, Hi.nnz
  <class 'scipy.sparse.coo.coo_matrix'> (28, 28) 27

As before, all Hessian matrices have stored only the lower triangular part. In addition to being possible to extract Hessian matrices of individual nonlinear constraints, it is also possible to construct any linear combination of these individual Hessian matrices. This can be done using the :class:`Constraint <pfnet.Constraint>` class method :func:`combine_H() <pfnet.Constraint.combine_H>`. After this, the resulting matrix can be extracted from the :data:`H_combined <pfnet.Constraint.H_combined>` attribute::

  >>> coefficients = np.random.randn(f.size)

  >>> constr.combine_H(coefficients)
  >>> H = constr.H_combined

  >>> print type(H), H.shape, H.nnz
  <class 'scipy.sparse.coo.coo_matrix'> (28, 28) 564

Lastly, Lagrange multiplier estimates of the linear and nonlinear constraints can be used to store sensitivity information in the network components associated with the constraints. This is done using the class method :func:`store_sensitivities() <pfnet.Constraint.store_sensitivities>`. Component-specific attributes that store sensitivity information are described in the :ref:`reference` section.

Details about each of the different constraint types available in PFNET are provided below.

.. _prob_constr_ACPF:

AC Power balance
----------------

This constraint is of type :data:`CONSTR_TYPE_PF <pfnet.CONSTR_TYPE_PF>`. It enforces active and reactive power balance at every bus of the network. It is given by

.. math:: 
   
   (P^g_k + j Q^g_k) - (P^l_k + j Q^l_k) - S_k^{sh} - \sum_{m \in [n]} S_{km} = 0, \ \forall \ k \in [n],

where :math:`P^g` and :math:`Q^g` are generator active and reactive powers, :math:`P^l` and :math:`Q^l` are load active and reactive powers, :math:`S^{sh}` are apparent powers flowing out of buses through shunt devices, :math:`S` are apparent powers flowing out of buses through branches, :math:`n` is the number of buses, and :math:`[n] := \{1,\ldots,n\}`. 

.. _prob_constr_DCPF:

DC Power balance
----------------

This constraint is of type :data:`CONSTR_TYPE_DCPF <pfnet.CONSTR_TYPE_DCPF>`. It enforces "DC" active power balance at every bus of the network. It is given by

.. math:: 
   
   P^g_k - P^l_k + \sum_{m \in [n]} b_{km} \left( \theta_k - \theta_m - \phi_{km} \right) = 0, \ \forall \ k \in [n],

where :math:`P^g` are generator active powers, :math:`P^l` are load active powers, :math:`b_{km}` are branch susceptances, :math:`\theta_k` are bus voltage angles, :math:`\phi_{km}` are phase shifts of phase-shifting transformers, :math:`n` is the number of buses, and :math:`[n] := \{1,\ldots,n\}`.

.. _prob_constr_DC_FLOW_LIM:

Branch DC power flow limits
---------------------------

This constraint is of type :data:`CONSTR_TYPE_DC_FLOW_LIM <pfnet.CONSTR_TYPE_DC_FLOW_LIM>`. It enforces branch "DC" power flow limits due to thermal ratings. It is given by

.. math:: 

   -P^{\max}_{km} \le -b_{km} \left( \theta_k - \theta_m - \phi_{km} \right) \le P^{\max}_{km},

for each branch :math:`(k,m)`, where :math:`b_{km}` are branch susceptances, :math:`\theta_k` are bus voltage angles, :math:`\phi_{km}` are phase shifts of phase-shifting transformers, and :math:`P^{\max}_{km}` are branch power flow limits. 

.. _prob_constr_FIX:

Variable fixing
---------------

This constraint is of type :data:`CONSTR_TYPE_FIX <pfnet.CONSTR_TYPE_FIX>`. It constrains specific variables to be fixed at their current value. The variables to be fixed are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` with the flag type :data:`FLAG_FIXED <pfnet.FLAG_FIXED>`.

.. _prob_constr_BOUND:

Variable bounding
-----------------

This constraint is of type :data:`CONSTR_TYPE_BOUND <pfnet.CONSTR_TYPE_BOUND>`. It constrains specific variables to be inside their bounds. The variables to be bounded are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` with the flag type :data:`FLAG_BOUNDED <pfnet.FLAG_BOUNDED>`. These constraints are expressed as nonlinear equality constraints using the techniques described in Section 4.3.3 of [TTR2015]_.

For conventional linear bounds, the constraint type :data:`CONSTR_TYPE_LBOUND <pfnet.CONSTR_TYPE_LBOUND>` can be used.

.. _prob_constr_PAR_GEN:

Generator participation
-----------------------

This constraint is of type :data:`CONSTR_TYPE_PAR_GEN <pfnet.CONSTR_TYPE_PAR_GEN>`. It enforces specific active power participations among slack generators, and reactive power participations among generators regulating the same bus voltage magnitude. For slack generators, all participate with equal active powers. For voltage regulating generators, each one participates with the same fraction of its total resources. More specifically, this constraint enforces

.. math:: 

   P^g_k = P^g_m,

for all slack generators :math:`k` and :math:`m` connected to the same bus, and

.. math::

   \frac{Q^g_k - Q^{\min}_k}{Q^{\max}_k - Q^{\min}_k} = \frac{Q^g_m - Q^{\min}_m}{Q^{\max}_m - Q^{\min}_m},

for all generators :math:`k` and :math:`m` regulating the same bus voltage magnitude, where :math:`Q^{\min}` and :math:`Q^{\max}` are generator reactive power limits.

.. _prob_constr_REG_GEN:

Voltage set-point regulation by generators
------------------------------------------

This constraint is of type :data:`CONSTR_TYPE_REG_GEN <pfnet.CONSTR_TYPE_REG_GEN>`. It enforces voltage set-point regulation by generators. It approximates the constraints

.. math:: 
   
   v_k & = v_k^t + v^y_k - v^z_k \\
   0 & \le (Q_k - Q^{\min}_k) \perp v^y_k \ge 0 \\
   0 & \le (Q^{\max}_k - Q_k) \perp v^z_k \ge 0,

for each bus :math:`k` whose voltage is regulated by generators, where :math:`v` are bus voltage magnitudes, :math:`v^t` are their set points, :math:`v^y` and :math:`v^z` are positive and negative deviations of :math:`v` from :math:`v^t`, and :math:`Q`, :math:`Q^{\max}` and :math:`Q^{\min}` are aggregate reactive powers and limits of the generators regulating the same bus voltage magnitude.

.. _prob_constr_REG_TRAN:

Voltage band regulation by transformers
---------------------------------------

This constraint is of type :data:`CONSTR_TYPE_REG_TRAN <pfnet.CONSTR_TYPE_REG_TRAN>`. It enforces voltage band regulation by tap-changing transformers. It approximates the constraints

.. math:: 
   
   t_k & = t_k^0 + t^y_k - t^z_k \\
   0 & \le (v_k + v^l_k - v^{\min}_k) \perp t^y_k \ge 0 \\
   0 & \le (v^{\max}_k - v_k + v^h_k) \perp t^z_k \ge 0 \\
   0 & \le (t^{\max}_k - t_k) \perp v^l_k \ge 0 \\
   0 & \le (t_k - t^{\min}_k) \perp v^h_k \ge 0,

for each bus :math:`k` whose voltage is regulated by tap-changing transformers, where :math:`v` are bus voltage magnitudes, :math:`v^{\max}` and :math:`v^{\min}` are their band limits, :math:`v^l` and :math:`v^h` are voltage violations of band lower and upper limits, :math:`t` are transformer tap ratios, :math:`t^0`, :math:`t^{\max}` and :math:`t^{\min}` are their current values and limits, and :math:`t^y` and :math:`t^z` are positive and negative deviations of :math:`t` from :math:`t^0`. The above equations assume that the sensitivity between voltage magnitude and transformer tap ratio is positive. If it is negative, :math:`t^y` and :math:`t^z` are interchanged in the first two complementarity constraints, and :math:`v^l` and :math:`v^h` are interchanged in the bottom two complementarity constraints. 

.. _prob_constr_REG_SHUNT:

Voltage band regulation by switched shunts
------------------------------------------

This constraint is of type :data:`CONSTR_TYPE_REG_SHUNT <pfnet.CONSTR_TYPE_REG_SHUNT>`. It enforces voltage band regulation by switched shunt devices. It approximates the constraints

.. math:: 
   
   b_k & = b_k^0 + b^y_k - b^z_k \\
   0 & \le (v_k + v^l_k - v^{\min}_k) \perp b^y_k \ge 0 \\
   0 & \le (v^{\max}_k - v_k + v^h_k) \perp b^z_k \ge 0 \\
   0 & \le (b^{\max}_k - b_k) \perp v^l_k \ge 0 \\
   0 & \le (b_k - b^{\min}_k) \perp v^h_k \ge 0,

for each bus :math:`k` whose voltage is regulated by switched shunt devices, where :math:`v` are bus voltage magnitudes, :math:`v^{\max}` and :math:`v^{\min}` are their band limits, :math:`v^l` and :math:`v^h` are voltage violations of band lower and upper limits, :math:`b` are switched shunt susceptances, :math:`b^0`, :math:`b^{\max}` and :math:`b^{\min}` are their current values and limits, and :math:`b^y` and :math:`b^z` are positive and negative deviations of :math:`b` from :math:`b^0` .

.. _prob_prob:

Problems
========

Optimization problems constructed with PFNET are of the form

.. math:: 
   :nowrap:

   \begin{alignat*}{2}
   & \mbox{minimize}   \quad && \varphi(x) \\
   & \mbox{subject to} \quad && Ax = b \\
   &                   \quad && f(x) = 0 \\
   &                   \quad && l \le Gx \le u,
   \end{alignat*}

As already noted, the objective function :math:`\varphi` is a weighted sum of functions :math:`\varphi_i`. The linear and nonlinear constraints :math:`Ax = b`, :math:`l \le Gx \le u`, and :math:`f(x) = 0` correspond to one or more of the constraints described above. An optimization problem in PFNET is represented by an object of type :class:`Problem <pfnet.Problem>`. 

After instantiation, a :class:`Problem <pfnet.Problem>` is empty and one needs to specify the :class:`Network <pfnet.Network>` that is to be associated with the problem, the :class:`Constraints <pfnet.Constraint>` to include, and the :class:`Functions <pfnet.Function>` that form the objective function. This can be done using the :class:`Problem <pfnet.Problem>` class methods :func:`set_network() <pfnet.Problem.set_network>`, :func:`add_constraint() <pfnet.Problem.add_constraint>`, and :func:`add_function() <pfnet.Problem.add_function>`. The following example shows how to construct a simple power flow problem and solve it using the Newton-Raphson method:

.. literalinclude:: ../examples/ex6.py

The above routine can then be used as follows::

   >>> net = Network()
   >>> net.load('case3012wp.mat')

   >>> print net.bus_P_mis, net.bus_Q_mis
   2.79e+0 1.56e+1
   
   >>> NRsolve(net)

   >>> print net.bus_P_mis, net.bus_Q_mis
   2.37e-6 3.58e-6

As shown in the example, the :class:`Problem <pfnet.Problem>` class method :func:`analyze() <pfnet.Problem.analyze>` needs to be called before the vectors and matrices associated with the problem constraints and functions can be used. The method :func:`eval() <pfnet.Problem.eval>` can then be used for evaluating the problem objective and constraint functions at different points. As is the case for :class:`Constraints <pfnet.Constraint>`, a :class:`Problem <pfnet.Problem>` has a method :func:`combine_H() <pfnet.Problem.combine_H>` for forming linear combinations of individual constraint Hessians, and a method :func:`store_sensitivities() <pfnet.Problem.store_sensitivities>` for storing sensitivity information in the network components associated with the constraints.
