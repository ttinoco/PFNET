.. include:: defs.hrst

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

where :math:`w_i` are weights, :math:`\varphi_i` are general linear or nonlinear functions, and :math:`x` is a vector of network variables. Each weight-function pair in the summation is represented by an object of type |FunctionBase|. An instance of this type can be constructed from the :class:`Function <pfnet.Function>` class, which requires specifying the function name, weight, and the |Network| object to be associated with the function. The following example sets all bus voltage magnitudes as variables and constructs a function that penalizes voltage magnitude deviations from ideal values::

  >>> import pfnet
  
  >>> net = net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> net.set_flags('bus',
  ...               'variable',
  ...               'any',
  ...               'voltage magnitude')

  >>> func = pfnet.Function('voltage magnitude regularization', 0.3, net)

  >>> print func.name == 'voltage magnitude regularization'
  True

  >>> print func.weight
  0.3

After a :class:`FunctionBase <pfnet.FunctionBase>` object is created, its value, gradient, and Hessian are zero, an empty vector, and an empty matrix, respectively. Before evaluating the function at a specific vector of values, it must be analyzed using the :class:`FunctionBase <pfnet.FunctionBase>` class method :func:`analyze() <pfnet.FunctionBase.analyze>`. This routine analyzes the function and allocates the required vectors and matrices for storing its gradient and Hessian. After this, the function can be evaluated using the method :func:`eval() <pfnet.FunctionBase.eval>`::

  >>> x = net.get_var_values()

  >>> func.analyze()

  >>> func.eval(x)

The value :math:`\varphi_i(x)`, gradient :math:`\nabla \varphi_i(x)`, and Hessian :math:`\nabla^2 \varphi_i(x)` of a function can then be extracted from the :data:`phi <pfnet.FunctionBase.phi>`, :data:`gphi <pfnet.FunctionBase.gphi>` and :data:`Hphi <pfnet.FunctionBase.Hphi>` attributes, respectively::

  >>> print x.shape
  (14,)

  >>> print func.phi
  0.255

  >>> print type(func.gphi), func.gphi.shape
  <type 'numpy.ndarray'> (14,)

  >>> print type(func.Hphi), func.Hphi.shape
  <class 'scipy.sparse.coo.coo_matrix'> (14, 14)

For the Hessian matrix, only the lower triangular part is stored.

Details about each of the different functions available in PFNET are provided below.

.. _prob_func_GEN_COST:

Active power generation cost
----------------------------

This function is associated with the string ``"generation cost"``. It measures active power generation cost by the expression

.. math::

   \varphi(x) := \sum_t \sum_k q_{k0} + q_{k1} P_k(t) + q_{k2} P_k(t)^2,

where :math:`P_k(t)` are generator active powers in per unit system base power, :math:`t` are time periods, and :math:`q_{k0}`, :math:`q_{k1}`, and :math:`q_{k2}` are constant coefficients. These coefficients correspond to the attributes :data:`cost_coeff_Q0 <pfnet.Generator.cost_coeff_Q0>`, :data:`cost_coeff_Q1 <pfnet.Generator.cost_coeff_Q1>`, and :data:`cost_coeff_Q2 <pfnet.Generator.cost_coeff_Q2>` of each :class:`Generator <pfnet.Generator>` object. 

.. _prob_func_LOAD_UTIL:

Active power consumption utility
--------------------------------

This function is associated with the string ``"consumption utility"``. It measures active power consumption utility by the expression

.. math::

   \varphi(x) := \sum_t \sum_k q_{k0} + q_{k1} P_k(t) + q_{k2} P_k(t)^2,

where :math:`P_k(t)` are load active powers in per unit system base power, :math:`t` are time periods, and :math:`q_{k0}`, :math:`q_{k1}`, and :math:`q_{k2}` are constant coefficients. These coefficients correspond to the attributes :data:`util_coeff_Q0 <pfnet.Load.util_coeff_Q0>`, :data:`util_coeff_Q1 <pfnet.Load.util_coeff_Q1>`, and :data:`util_coeff_Q2 <pfnet.Load.util_coeff_Q2>` of each :class:`Load <pfnet.Load>` object.

.. _prob_func_NETCON_COST:

Net Active Power Consumption Cost
---------------------------------

This function is associated with the string ``"net consumption cost"``. It measures the total cost of net active power consumption over the time periods using the price defined by the :data:`price <pfnet.Bus.price>` attribute of each :class:`Bus <pfnet.Bus>` object.

.. _prob_func_REG_VMAG:

Voltage magnitude regularization
--------------------------------

This function is associated with the string ``"voltage magnitude regularization"``. It penalizes deviations of bus voltage magnitudes from ideal values. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{v_k(t) - v^s_k(t)}{\Delta v} \Bigg)^2,

where :math:`t` are time periods, :math:`v` are bus voltage magnitudes, :math:`v^s` are voltage magnitude set points (one for buses not regulated by generators), and :math:`\Delta v` is a normalization factor. 

.. _prob_func_REG_ANG:

Voltage angle regularization
----------------------------

This function is associated with the string ``"voltage angle regularization"``. It penalizes large bus voltage angles and voltage angle differences across branches. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{\theta_k(t)}{\Delta \theta} \Bigg)^2 + 
                 \frac{1}{2} \sum_t \sum_{(k,m)} \Bigg( \frac{\theta_k(t) - \theta_m(t) - \phi_{km}(t)}{\Delta \theta} \Bigg)^2,

where :math:`t` are time periods, :math:`\theta` are bus voltage angles, :math:`\phi` are branch phase shifts, and :math:`\Delta \theta` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_PQ:

Generator powers regularization
-------------------------------

This function is associated with the string ``"generator powers regularization"``. It penalizes deviations of generator powers from the midpoint of their ranges. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{P^g_k(t) - \bar{P}_k}{\Delta P} \Bigg)^2 + 
                 \frac{1}{2} \sum_t \sum_k \Bigg( \frac{Q^g_k(t) - \bar{Q}_k}{\Delta Q} \Bigg)^2,

where :math:`t` are time periods, :math:`P^g` and :math:`Q^g` are generator active and reactive powers, :math:`\bar{P}` and :math:`\bar{Q}` are midpoints of generator active and reactive power ranges, and :math:`\Delta P = \Delta Q` are normalization factors. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_RATIO:

Transformer tap ratio regularization
------------------------------------

This function is associated with the string ``"tap ratio regularization"``. It penalizes deviations of tap ratios of tap-changing transformers from their initial values. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_{\tau} \sum_k \Bigg( \frac{t_k(\tau) - t^0_k(\tau)}{\Delta t} \Bigg)^2,

where :math:`\tau` are time periods, :math:`t` are tap ratios of tap-changing transformers, :math:`t^0` are their initial values, and :math:`\Delta t` is a normalization factor.

.. _prob_func_REG_PHASE:

Transformer phase shift regularization
--------------------------------------

This function is associated with the string ``"phase shift regularization"``. It penalizes deviations of phase shifts of phase shifting transformers from their initial value. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{\phi_k(t) - \phi^0_k(t)}{\Delta \phi} \Bigg)^2

where :math:`t` are time periods, :math:`\phi` are phase shifts of phase-shifting transformers, :math:`\phi^0` are their initial values, and :math:`\Delta \phi` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_REG_SUSC:

Switched shunt susceptance regularization
-----------------------------------------

This function is associated with the string ``"susceptance regularization"``. It penalizes deviations of susceptances of switched shunt devices from their initial values. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{b_k(t) - b^0_k(t)}{\Delta b} \Bigg)^2,

where :math:`t` are time periods, :math:`b` are susceptances of switched shunt devices, :math:`b^0` are their initial values, and :math:`\Delta b` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_SLIM_VMAG:

Voltage magnitude soft limit penalty
------------------------------------

This function is associated with the string ``"soft voltage magnitude limits"``. It penalizes deviations of bus voltage magnitudes from the mid point of their ranges. It is defined by the expression

.. math::

   \varphi(x) := \frac{1}{2} \sum_t \sum_k \Bigg( \frac{v_k(t) - \bar{v}_k}{\Delta v} \Bigg)^2,

where :math:`t` are time periods, :math:`v` are bus voltage magnitudes, :math:`\bar{v}` are the mid points of their ranges, and :math:`\Delta v` is a normalization factor. Only terms that include optimization variables are included in the summation.

.. _prob_func_SP_CONTROLS:

Sparsity inducing penalty for controls
--------------------------------------

This function is associated with the string ``"sparse controls penalty"``. It encourages sparse control adjustments with the expression

.. math::

   \varphi(x) := \sum_t \sum_k \sqrt{ \Bigg( \frac{u_k(t) - u_k^0(t)}{\Delta u_k} \Bigg)^2 + \epsilon },

where :math:`t` are time periods, :math:`u` are control quantities, :math:`u^0` are their current values, and :math:`\epsilon` is a small positive scalar. The normalization factors :math:`\Delta u_k` are given by

.. math::

   \Delta u_k := \max\{u^{\max}_k-u^{\min}_k, \delta\},

where :math:`u^{\max}` and :math:`u^{\min}` are control limits, and :math:`\delta` is a small positive scalar. The control quantities that are considered by this function are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` using the flag ``"sparse"``.

.. _prob_constr:

Constraints
===========

Constraints in PFNET are of the form

.. math::
   
   A \left[ \begin{array}{c} x \\ y \end{array} \right] = b, \quad
   f(x,y) = 0, \quad
   l \le G \left[ \begin{array}{c} x \\ y \end{array} \right] \le u,

where :math:`A` and :math:`G`  are sparse matrices, :math:`b`, :math:`l` and :math:`u`  are vectors, :math:`f` is a vector-valued nonlinear function, and :math:`x` and :math:`y` are vectors of network variables and extra or auxiliary constraint variables, respectively. They are represented by objects of type :class:`ConstraintBase <pfnet.ConstraintBase>`. An instance of this type can be constructed from the class :class:`Constraint <pfnet.Constraint>`, whose constructor requires specifying the constraint name and the |Network| to be associated with the constraint. The following example sets all bus voltage magnitudes and angles as variables and constructs the AC power balance constraints::

  >>> import pfnet

  >>> pfnet.ParserMAT().parse('ieee14.mat')

  >>> net.set_flags('bus',
  ...               'variable',
  ...               'any',
  ...               ['voltage magnitude','voltage angle'])

  >>> print net.num_vars == 2*net.num_buses
  True

  >>> constr = pfnet.Constraint('AC power balance',net)

  >>> print constr.name == 'AC power balance'
  True

Before a :class:`ConstraintBase <pfnet.ConstraintBase>` object can be used, it must be initialized using the :class:`ConstraintBase <pfnet.ConstraintBase>` class method :func:`analyze() <pfnet.ConstraintBase.analyze>`. This routine analyzes the constraint and allocates the required vectors and matrices. After this, the number of extra variables can be obtained from the attribute :data:`num_extra_vars <pfnet.ConstraintBase.num_extra_vars>`, and the constraint can be evaluated using the method :func:`eval() <pfnet.ConstraintBase.eval>`::

  >>> x = net.get_var_values()

  >>> constr.analyze()

  >>> print constr.num_extra_vars
  0

  >>> constr.eval(x)

The matrices and vectors associated with the linear constraints can be extracted from the :data:`A <pfnet.ConstraintBase.A>`, :data:`G <pfnet.ConstraintBase.G>`, :data:`b <pfnet.ConstraintBase.b>`, :data:`l <pfnet.ConstraintBase.l>` and :data:`u <pfnet.ConstraintBase.u>` attributes of the :class:`ConstraintBase <pfnet.ConstraintBase>` object. The vector of violations and Jacobian matrix of the nonlinear constraints can be extracted from the attributes :data:`f <pfnet.ConstraintBase.f>` and :data:`J <pfnet.ConstraintBase.J>`, respectively. Also, the Hessian matrix of any individual nonlinear constraint function :math:`f_i(x)` can be extracted using the class method :func:`get_H_single() <pfnet.ConstraintBase.get_H_single>`. The following example shows how to extract the largest power flow mismatch in per unit :data:`system base power <pfnet.Network.base_power>` and the Hessian matrix corresponding to the active power balance constraint of a bus::

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

As before, all Hessian matrices have stored only the lower triangular part. In addition to being possible to extract Hessian matrices of individual nonlinear constraints, it is also possible to construct any linear combination of these individual Hessian matrices. This can be done using the :class:`ConstraintBase <pfnet.ConstraintBase>` class method :func:`combine_H() <pfnet.ConstraintBase.combine_H>`. After this, the resulting matrix can be extracted from the :data:`H_combined <pfnet.ConstraintBase.H_combined>` attribute::

  >>> coefficients = np.random.randn(f.size)

  >>> constr.combine_H(coefficients)
  >>> H = constr.H_combined

  >>> print type(H), H.shape, H.nnz
  <class 'scipy.sparse.coo.coo_matrix'> (28, 28) 564

Lagrange multiplier estimates of the linear and nonlinear constraints can be used to store sensitivity information in the network components associated with the constraints. This is done using the class method :func:`store_sensitivities() <pfnet.ConstraintBase.store_sensitivities>`. Component-specific attributes that store sensitivity information are described in the :ref:`reference` section.

Lastly, the |ConstraintBase| class provides methods for extracting information about the rows of each linear and nonlinear constraint. These methods are :func:`get_A_row_info_string() <pfnet.ConstraintBase.get_A_row_info_string>`, :func:`get_J_row_info_string() <pfnet.ConstraintBase.get_J_row_info_string>`, and :func:`get_G_row_info_string() <pfnet.ConstraintBase.get_G_row_info_string>`.

Details about each of the different constraints available in PFNET are provided below.

.. _prob_constr_ACPF:

AC Power balance
----------------

This constraint is associated with the string ``"AC power balance"``. It enforces "AC" active and reactive power balance at every bus of the network. It is given by

.. math:: 
   
   P^g_k(t) + j Q^g_k(t) - P^l_k(t) - j Q^l_k(t) - S_k^{sh}(t) - \sum_{m \in [n]} S_{km}(t) = 0, \ \forall \ k \in [n], \ t \in [T],

where :math:`t` are time periods, :math:`P^g` and :math:`Q^g` are generator active and reactive powers, :math:`P^l` and :math:`Q^l` are load active and reactive powers, :math:`S^{sh}` are apparent powers flowing out of buses through shunt devices, :math:`S` are apparent powers flowing out of buses through branches, :math:`n` is the number of buses, :math:`T` is the number of time periods, and :math:`[n] := \{1,\ldots,n\}`. 

.. _prob_constr_DCPF:

DC Power balance
----------------

This constraint is associated with the string ``"DC power balance"``. It enforces "DC" active power balance at every bus of the network. It is given by

.. math:: 
   
   P^g_k(t) - P^l_k(t) + \sum_{m \in [n]} b_{km} \left( \theta_k(t) - \theta_m(t) - \phi_{km}(t) \right) = 0, \ \forall \ k \in [n], \ t \in [T], 

where :math:`t` are time periods, :math:`P^g` are generator active powers, :math:`P^l` are load active powers, :math:`b_{km}` are branch susceptances, :math:`\theta_k` are bus voltage angles, :math:`\phi_{km}` are phase shifts of phase-shifting transformers, :math:`n` is the number of buses, :math:`T` is the number of time periods, and :math:`[n] := \{1,\ldots,n\}`.

.. _prob_constr_LINPF:

Linearized AC Power balance
---------------------------

This constraint is associated with the string ``"linearized AC power balance"``. It enforces active and reactive power balance at every bus of the network using a first-order Taylor expansion of the AC power balance constraints. It is given by

.. math:: 
   
   J(x_0) x = J(x_0) x_0 - f(x_0),

where :math:`x_0` is the vector of current variable values, :math:`f(x_0)` is the vector of AC bus power mismatches, and :math:`J(x_0)` is the Jacobian of :math:`f` at :math:`x_0`.

.. _prob_constr_DC_FLOW_LIM:

DC branch flow limits
---------------------

This constraint is associated with the string ``"DC branch flow limits"``. It enforces branch "DC" power flow limits due to thermal ratings. It is given by

.. math:: 

   -P^{\max}_{km} \le -b_{km} \left( \theta_k(t) - \theta_m(t) - \phi_{km}(t) \right) \le P^{\max}_{km},

for each branch :math:`(k,m)` and time period :math:`t`, where :math:`b_{km}` are branch susceptances, :math:`\theta_k` are bus voltage angles, :math:`\phi_{km}` are phase shifts of phase-shifting transformers, and :math:`P^{\max}_{km}` are branch power flow limits. 

.. _prob_constr_AC_FLOW_LIM:

AC branch flow limits
---------------------

This constraint is associated with the string ``"AC branch flow limits"``. It enforces branch "AC" power flow limits due to thermal ratings based on current magnitudes. It utilizes auxiliary variables (slacks). It is given by

.. _prob_constr_AC_LIN_FLOW_LIM:

Linearized AC branch flow limits
--------------------------------

This constraint is associated with the string ``"linearized AC branch flow limits"``. It enforces branch "AC" power flow limits due to thermal ratings based on current magnitudes using conservative linear constraints constructed with the external library ``Line Flow`` [DS2017]_. This external library is used with its default parameters.

.. _prob_constr_FIX:

Variable fixing
---------------

This constraint is associated with the string ``"variable fixing"``. It constrains specific variables to be fixed at their current value. The variables to be fixed are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` with the flag ``"fixed"``.

.. _prob_constr_BOUND:

Variable bounds
---------------

This constraint is associated with the string ``"variable bounds"``. It constrains specific variables to be inside their bounds. The variables to be bounded are specified using the :class:`Network <pfnet.Network>` class methods :func:`set_flags() <pfnet.Network.set_flags>` or :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` with the flag ``"bounded"``.

.. _prob_constr_PAR_GEN:

Generator participation
-----------------------

This constraint is associated with the string ``"generator active power participation"`` and ``"generator reactive power participation"``. It enforces specific active power participations among slack generators connected to the same bus, or reactive power participations among generators regulating the same bus voltage magnitude. For slack generators, all participate with equal active powers. For voltage-regulating generators, each one participates with the same fraction of its total reactive resources. More specifically, this constraint enforces

.. math:: 

   P^g_k(t) = P^g_m(t),

for all slack generators :math:`k` and :math:`m` connected to the same bus and time period :math:`t`, or

.. math::

   \frac{Q^g_k(t) - Q^{\min}_k}{Q^{\max}_k - Q^{\min}_k} = \frac{Q^g_m(t) - Q^{\min}_m}{Q^{\max}_m - Q^{\min}_m},

for all generators :math:`k` and :math:`m` regulating the same bus voltage magnitude and time period :math:`t`, where :math:`Q^{\min}` and :math:`Q^{\max}` are generator reactive power limits.

.. _prob_constr_REG_GEN:

Voltage set-point regulation by generators
------------------------------------------

This constraint is associated with the string ``"voltage regulation by generators"``. It enforces voltage set-point regulation by generators. It approximates the constraints

.. math:: 
   
   v_k(t) & = v_k^s(t) + v^y_k(t) - v^z_k(t) \\
   0 & \le (Q_k(t) - Q^{\min}_k) \perp v^y_k(t) \ge 0 \\
   0 & \le (Q^{\max}_k - Q_k(t)) \perp v^z_k(t) \ge 0,

for each time period :math:`t` and bus :math:`k` whose voltage is regulated by generators, where :math:`v` are bus voltage magnitudes, :math:`v^s` are their set points, :math:`v^y` and :math:`v^z` are positive and negative deviations of :math:`v` from :math:`v^s` (auxiliary variables of the constraint), and :math:`Q`, :math:`Q^{\max}` and :math:`Q^{\min}` are aggregate reactive powers and limits of the generators regulating the same bus voltage magnitude.

.. _prob_constr_REG_TRAN:

Voltage band regulation by transformers
---------------------------------------

This constraint is associated with the string ``"voltage regulation by transformers"``. It enforces voltage band regulation by tap-changing transformers. It approximates the constraints

.. math:: 
   
   t_k(\tau) & = t_k^0(\tau) + t^y_k(\tau) - t^z_k(\tau) \\
   0 & \le (v_k(\tau) + v^l_k(\tau) - v^{\min}_k) \perp t^y_k(\tau) \ge 0 \\
   0 & \le (v^{\max}_k - v_k(\tau) + v^h_k(\tau)) \perp t^z_k(\tau) \ge 0 \\
   0 & \le (t^{\max}_k - t_k(\tau)) \perp v^l_k(\tau) \ge 0 \\
   0 & \le (t_k(\tau) - t^{\min}_k) \perp v^h_k(\tau) \ge 0,

for each time period :math:`\tau` and bus :math:`k` whose voltage is regulated by tap-changing transformers, where :math:`v` are bus voltage magnitudes, :math:`v^{\max}` and :math:`v^{\min}` are their band limits, :math:`v^l` and :math:`v^h` are voltage violations of band lower and upper limits (auxiliary variables), :math:`t` are transformer tap ratios, :math:`t^0`, :math:`t^{\max}` and :math:`t^{\min}` are their current values and limits, and :math:`t^y` and :math:`t^z` are positive and negative deviations of :math:`t` from :math:`t^0` (auxiliary variables). The above equations assume that the sensitivity between voltage magnitude and transformer tap ratio is positive. If it is negative, :math:`t^y` and :math:`t^z` are interchanged in the first two complementarity constraints, and :math:`v^l` and :math:`v^h` are interchanged in the bottom two complementarity constraints.

.. _prob_constr_REG_SHUNT:

Voltage band regulation by switched shunts
------------------------------------------

This constraint is associated with the string ``"voltage regulation by shunts"``. It enforces voltage band regulation by switched shunt devices. It approximates the constraints

.. math:: 
   
   b_k(t) & = b_k^0(t) + b^y_k(t) - b^z_k(t) \\
   0 & \le (v_k(t) + v^l_k(t) - v^{\min}_k) \perp b^y_k(t) \ge 0 \\
   0 & \le (v^{\max}_k - v_k(t) + v^h_k(t)) \perp b^z_k(t) \ge 0 \\
   0 & \le (b^{\max}_k - b_k(t)) \perp v^l_k(t) \ge 0 \\
   0 & \le (b_k(t) - b^{\min}_k) \perp v^h_k(t) \ge 0,

for each time period :math:`t` and bus :math:`k` whose voltage is regulated by switched shunt devices, where :math:`v` are bus voltage magnitudes, :math:`v^{\max}` and :math:`v^{\min}` are their band limits, :math:`v^l` and :math:`v^h` are voltage violations of band lower and upper limits (auxiliary variables), :math:`b` are switched shunt susceptances, :math:`b^0`, :math:`b^{\max}` and :math:`b^{\min}` are their current values and limits, and :math:`b^y` and :math:`b^z` are positive and negative deviations of :math:`b` from :math:`b^0` (auxiliary variables).

.. _prob_constr_GEN_RAMP:

Generator active power ramp limits
----------------------------------

This constraint is associated with the string ``"generator ramp limits"``. It enforces generator active power ramping limits. It is given by

.. math:: 
   
    -\delta P^{\max}_k \le P^g_k(t) - P^g_k(t-1) \le \delta P^{\max}_k

for each generator :math:`k` and time period :math:`t \in \{1,\ldots,T\}`, where :math:`P^g` are generator active powers, and :math:`\delta P^{\max}` are generator ramping limits. The ramping limits are defined by the :data:`dP_max <pfnet.Generator.dP_max>` attribute of each :class:`Generator <pfnet.Generator>` object. For :math:`t = 1`, :math:`P^g_k(t-1)` is the :data:`P_prev <pfnet.Generator.P_prev>` attribute of a :class:`Generator <pfnet.Generator>`.

.. _prob_constr_BAT_DYN:

Battery dynamics
----------------

This constraint is associated with the string ``"battery dynamics"``. It enforces the dynamic equations of the batteries' energy levels. It is given by

.. math::
   :nowrap:

   \begin{align*}
   E_k(1) &= E_k^i \\
   E_k(T+1) &= E_k^f \\
   E_k(t+1) &= E_k(t) + \eta_k^c P_k^c(t) - (\eta_k^d)^{-1} P_k^d(t), \ \forall t \in \{1,\ldots,T\}, 
   \end{align*}

for each battery :math:`k`, where :math:`E`, :math:`P^c`, and :math:`P^d` denote battery energy levels and charging and discharging powers, respectively, and :math:`E^i`, :math:`E^f`, :math:`\eta^c`, and :math:`\eta^d` correspond to the attributes :data:`E_init <pfnet.Battery.E_init>` , :data:`E_final <pfnet.Battery.E_final>`, :data:`eta_c <pfnet.Battery.eta_c>`, and :data:`eta_d <pfnet.Battery.eta_d>` of a :class:`Battery <pfnet.Battery>`, respectively. It is noted here that the units of the charging/discharging powers are p.u. system base power, and the units of the energy levels are p.u. system base power times the duration of a time period.

.. _prob_constr_LOAD_PF:

Load Constant Power Factor
--------------------------

This constraint is associated with the string ``"load constant power factor"``. It forces loads to have a constant power factor at each time. It is given by

.. math:: 
   
    Q_k^l(t) - P_k^l(t) \frac{\sqrt{1-\gamma^2}}{|\gamma|} = 0

for each load :math:`k` and time period :math:`t \in \{1,\ldots,T\}`, where :math:`P^l` are active powers, :math:`Q^l` are reactive powers, and :math:`\gamma` is the load's :data:`target_power_factor <pfnet.Load.target_power_factor>`.

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

After instantiation, a :class:`Problem <pfnet.Problem>` is empty and one needs to specify the :class:`Constraints <pfnet.ConstraintBase>` to include, and the :class:`Functions <pfnet.FunctionBase>` that form the objective function. This can be done using the :class:`Problem <pfnet.Problem>` class methods :func:`add_constraint() <pfnet.Problem.add_constraint>`, and :func:`add_function() <pfnet.Problem.add_function>`, respectively. The following example shows how to construct a simple power flow problem and solve it using the Newton-Raphson method:

.. literalinclude:: ../examples/power_flow.py

The above routine can then be used as follows::

   >>> pfnet.ParserMAT().parse('case3012wp.mat')

   >>> print net.bus_P_mis, net.bus_Q_mis
   2.79e+0 1.56e+1
   
   >>> NRsolve(net)

   >>> print net.bus_P_mis, net.bus_Q_mis
   2.37e-6 3.58e-6

As shown in the example, the :class:`Problem <pfnet.Problem>` class method :func:`analyze() <pfnet.Problem.analyze>` needs to be called before the vectors and matrices associated with the problem constraints and functions can be used. The method :func:`eval() <pfnet.Problem.eval>` can then be used for evaluating the problem objective and constraint functions at different points. As is the case for :class:`Constraints <pfnet.ConstraintBase>`, a :class:`Problem <pfnet.Problem>` has a method :func:`combine_H() <pfnet.Problem.combine_H>` for forming linear combinations of individual constraint Hessians, and a method :func:`store_sensitivities() <pfnet.Problem.store_sensitivities>` for storing sensitivity information in the network components associated with the constraints.
