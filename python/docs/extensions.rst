.. _ext:

**********
Extensions
**********

This section describes how to add custom parsers, functions, and constraints to PFNET.

.. _ext_parser:

Adding a Parser
===============

To add a new data parser to PFNET, one should create a subclass of the :class:`CustomParser <pfnet.CustomParser>` class and provide five methods:

* :func:`init(self) <pfnet.ParserBase.init>`
  
  * Initializes parser data. This is called every time the :func:`parse() <pfnet.ParserBase.parse>` method is called.

* :func:`parse(self,filename,num_periods=1) <pfnet.ParserBase.parse>`

  * Reads power network data from ``filename`` and constructs and returns a :class:`Network <pfnet.Network>` with the given number of time periods. The custom parser only needs to allocate the :ref:`network components <net_components>` and store their information for time period zero. PFNET allocates flag and other utility arrays, copies the data to all time periods, and computes the initial network properties automatically. 

* :func:`set(self,key,value) <pfnet.ParserBase.set>`

  * Sets the parser's ``key`` parameter or option to ``value``, if it exists. 

* :func:`show(self) <pfnet.ParserBase.show>`

  * Shows information about the parser's internal data. 

* :func:`write(self,net,filename) <pfnet.ParserBase.write>`

  * Writes the given :class:`Network <pfnet.Network>` to ``filename`` using the format assoicated with the custom parser. 

A template for creating a custom parser is provided below:

.. literalinclude:: ../examples/custom_parser_template.py

.. _ext_func:

Adding a Function
=================

To add a new function to PFNET, one should create a subclass of the :class:`CustomFunction <pfnet.CustomFunction>` class and provide six methods:

* :func:`init(self) <pfnet.CustomFunction.init>`

  * This method initializes any custom function data. 

* :func:`count_step(self,branch,t) <pfnet.CustomFunction.count_step>`

  * This method is called for every :class:`branch <pfnet.Branch>` and time period, and is responsible for updating the counter :data:`Hphi_nnz <pfnet.FunctionBase.Hphi_nnz>` of nonzero entries of the Hessian matrix of the function (only lower triangular part). 

* :func:`clear(self) <pfnet.CustomFunction.clear>`

  * This method resets the value of any attribute that is updated by the other methods, such as the value of the function :data:`phi <pfnet.FunctionBase.phi>`, the counter :data:`Hphi_nnz <pfnet.FunctionBase.Hphi_nnz>`, etc. If used, the array of flags :class:`bus_counted <pfnet.FunctionBase.bus_counted>`, which has size equal to the number of buses times the number of time periods and can be used to keep track of which buses have already been processed, should also be reset here. 

* :func:`allocate(self) <pfnet.CustomFunction.allocate>`

  * This method allocates the gradient vector :data:`gphi <pfnet.FunctionBase.gphi>` and Hessian matrix :data:`Hphi <pfnet.FunctionBase.Hphi>` (lower triangular part) using the methods :func:`set_gphi() <pfnet.FunctionBase.set_gphi>` and :func:`set_Hphi() <pfnet.FunctionBase.set_Hphi>`, respectively.

* :func:`analyze_step(self,branch,t) <pfnet.CustomFunction.analyze_step>`

  * This method is called for every :class:`branch <pfnet.Branch>` and time period, and is responsible for storing the structural or constant information of the Hessian matrix :data:`Hphi <pfnet.FunctionBase.Hphi>` (lower triangular part). 

* :func:`eval_step(self,branch,t,x) <pfnet.CustomFunction.eval_step>`

  * This method is called for every :class:`branch <pfnet.Branch>` and time period, and is responsible for updating the values of :data:`phi <pfnet.FunctionBase.phi>`, :data:`gphi <pfnet.FunctionBase.gphi>`, and :data:`Hphi <pfnet.FunctionBase.Hphi>` using the given vector of variable values. 

A template for creating a custom function is provided below:

.. literalinclude:: ../examples/custom_function_template.py

An example of a custom function that computes the quadratic active power generation cost can be found in `here <https://github.com/ttinoco/PFNET/blob/master/python/pfnet/functions/dummy_function.py>`_. 

.. _ext_constr:

Adding a Constraint
===================

To add a new constraint to PFNET, one should create a subclass of the :class:`CustomConstraint <pfnet.CustomConstraint>` class. The subclass needs to define the following seven methods:

* :func:`init(self) <pfnet.CustomConstraint.init>`

  * This method initializes any custom constraint data. If the constraint has nonlinear equality constraints, the array :data:`H_nnz <pfnet.ConstraintBase.H_nnz>` of counters of nonzero elements of each constraint Hessian needs to be allocated using the method :func:`set_H_nnz <pfnet.ConstraintBase.set_H_nnz()>`. For off-diagonal pairs of elements, only one should be counted.

* :func:`count_step(self,branch,t) <pfnet.CustomConstraint.count_step>`

  * This method is called for every :class:`branch <pfnet.Branch>` and time period, and is responsible for updating the counters :data:`A_row <pfnet.ConstraintBase.A_row>`, :data:`A_nnz <pfnet.ConstraintBase.A_nnz>`, :data:`G_row <pfnet.ConstraintBase.G_row>`, :data:`G_nnz <pfnet.ConstraintBase.G_nnz>`, :data:`J_row <pfnet.ConstraintBase.J_row>`, and :data:`J_nnz <pfnet.ConstraintBase.J_nnz>`, which count the number of rows and nonzero elements of the matrices :data:`A <pfnet.ConstraintBase.A>`, :data:`G <pfnet.ConstraintBase.G>`, and Jacobian :data:`J <pfnet.ConstraintBase.J>`, respectively, and for updating the entries of the array :data:`H_nnz <pfnet.ConstraintBase.H_nnz>`. 

* :func:`clear(self) <pfnet.CustomConstraint.clear>`

  * This method resets the value of any attribute that is updated by the other methods, such as the counters :data:`A_row <pfnet.ConstraintBase.A_row>`, :data:`A_nnz <pfnet.ConstraintBase.A_nnz>`, :data:`G_row <pfnet.ConstraintBase.G_row>`, :data:`G_nnz <pfnet.ConstraintBase.G_nnz>`, :data:`J_row <pfnet.ConstraintBase.J_row>`, :data:`J_nnz <pfnet.ConstraintBase.J_nnz>`, and :data:`H_nnz <pfnet.ConstraintBase.H_nnz>`. If used, the array of flags :class:`bus_counted <pfnet.ConstraintBase.bus_counted>`, which has size equal to the number of buses times the number of time periods and can be used to keep track of which buses have already being processed, should also be reset here. 

* :func:`allocate(self) <pfnet.CustomConstraint.allocate>`

  * This method allocates the vectors :data:`b <pfnet.ConstraintBase.b>`, :data:`f <pfnet.ConstraintBase.f>`, :data:`l <pfnet.ConstraintBase.l>`, and :data:`u <pfnet.ConstraintBase.u>`, and matrices :data:`A <pfnet.ConstraintBase.A>`, :data:`G <pfnet.ConstraintBase.G>`, and Jacobian :data:`J <pfnet.ConstraintBase.J>` using the methods :func:`set_b() <pfnet.ConstraintBase.set_b>`, :func:`set_f() <pfnet.ConstraintBase.set_f>`, :func:`set_l() <pfnet.ConstraintBase.set_l>`, :func:`set_u() <pfnet.ConstraintBase.set_u>`, :func:`set_A() <pfnet.ConstraintBase.set_A>`, :func:`set_J() <pfnet.ConstraintBase.set_J>`, and :func:`set_G() <pfnet.ConstraintBase.set_G>`, respectively. For constraints with nonlinear equality constraint, this method should also allocate the array of constraint Hessians using the method :func:`allocate_H_array() <pfnet.ConstraintBase.allocate_H_array>` as well as the individual Hessians using the method :func:`set_H_single() <pfnet.ConstraintBase.set_H_single>` and information available in :data:`H_nnz <pfnet.ConstraintBase.H_nnz>`.

* :func:`analyze_step(self,branch,t) <pfnet.CustomConstraint.analyze_step>`

  * This method is called for every :class:`branch <pfnet.Branch>` and time period, and is responsible for storing the structural or constant information of the matrices :data:`A <pfnet.ConstraintBase.A>`, :data:`G <pfnet.ConstraintBase.G>`, Jacobian :data:`J <pfnet.ConstraintBase.J>`, and the Hessians of the nonlinear equality constraint functions. The latter can be extracted using the method :func:`get_H_single() <pfnet.ConstraintBase.get_H_single>`. For these Hessian matrices, only one element of each off-diagonal pairs should be stored. After all branches and time periods have been processed, PFNET automatically makes these constraint Hessians lower triangular by swapping elements as necessary.

* :func:`eval_step(self,branch,t,x) <pfnet.CustomConstraint.eval_step>`

  * This method is used for updating the values of the nonlinear constraint functions :data:`f <pfnet.ConstraintBase.f>`, their Jacobian :data:`J <pfnet.ConstraintBase.J>`, and Hessians.

* :func:`store_sens_step(self,branch,t,sA,sf,sGu,sGl) <pfnet.CustomConstraint.store_sens_step>`

  * This method is used for storing constraint sensitivity information in the :ref:`network components <net_components>` and will be documented in the near future. It can be ignored for now.

A template for creating a custom constraint is provided below:

.. literalinclude:: ../examples/custom_constraint_template.py

An example of a custom constraint that constructs the DC power balance equations can be found in `here <https://github.com/ttinoco/PFNET/blob/master/python/pfnet/constraints/dummy_constraint.py>`_.
