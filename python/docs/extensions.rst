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

  * This method is called for every :class:`network branch <pfnet.Branch>` and time period combination, and is responsible for updating the counter :data:`Hphi_nnz <pfnet.FunctionBase.Hphi_nnz>` of nonzero entries of the Hessian matrix of the function. 

* :func:`clear(self) <pfnet.CustomFunction.clear>`

  * This method resets the value of any attribute that is updated by the other methods, such as the value of the function :data:`phi <pfnet.FunctionBase.phi>`, the counter :data:`Hphi_nnz <pfnet.FunctionBase.Hphi_nnz>`, etc. If used, the array of flags :class:`bus_counted <pfnet.FunctionBase.bus_counted>`, which has size equal to the number of buses times the number of time periods and can be used to keep track of which buses have already being processed, should also be reset here. 

* :func:`allocate(self) <pfnet.CustomFunction.allocate>`

  * This method allocates the gradient vector :data:`gphi <pfnet.FunctionBase.gphi>` and Hessian matrix :data:`Hphi <pfnet.FunctionBase.Hphi>` (lower triangular part) using the methods :func:`set_gphi() <pfnet.FunctionBase.set_gphi>` and :func:`set_Hphi() <pfnet.FunctionBase.set_Hphi>`, respectively.

* :func:`analyze_step(self,branch,t) <pfnet.CustomFunction.analyze_step>`

  * This method is called for every :class:`network branch <pfnet.Branch>` and time period combination, and is responsible for storing the structural or constant information of the Hessian matrix :data:`Hphi <pfnet.FunctionBase.Hphi>` (lower triangular part). 

* :func:`eval_step(self,branch,t,x) <pfnet.CustomFunction.eval_step>`

  * This method is called for every :class:`network branch <pfnet.Branch>` and time period combination, and is responsible for updating the values of :data:`phi <pfnet.FunctionBase.phi>`, :data:`gphi <pfnet.FunctionBase.gphi>`, and :data:`Hphi <pfnet.FunctionBase.Hphi>` using the given vector of variable values. 

A template for creating a custom function is provided below:

.. literalinclude:: ../examples/custom_function_template.py

An example of a custom function that computes the quadratic active power generation cost can be found in `here <https://github.com/ttinoco/PFNET/blob/master/python/pfnet/functions/dummy_function.py>`_. 

.. _ext_constr:

Adding a Constraint
===================

Similarly, to add a new constraint to PFNET, one should create a subclass of the :class:`CustomConstraint <pfnet.CustomConstraint>` class. Currently, this documentation covers linear constraints only. Adding nonlinear constraints and constraints that use auxiliary variables will be documented in the near future. The subclass needs to define the following seven methods:

* :func:`init(self) <pfnet.CustomConstraint.init>`

  * This method initializes any custom constraint data. 

* :func:`count_step(self,branch,t) <pfnet.CustomConstraint.count_step>`

  * This method is called for every :class:`network branch <pfnet.Branch>` and time period combination, and is responsible for updating the counters :data:`A_row <pfnet.ConstraintBase.A_row>`, :data:`A_nnz <pfnet.ConstraintBase.A_nnz>`, :data:`G_row <pfnet.ConstraintBase.G_row>`, and :data:`G_nnz <pfnet.ConstraintBase.G_nnz>`, which count the number of rows and nonzero elements of the matrices :data:`A <pfnet.ConstraintBase.A>` and :data:`G <pfnet.ConstraintBase.G>`, respectively. 

* :func:`clear(self) <pfnet.CustomConstraint.clear>`

  * This method resets the value of any attribute that is updated by the other methods, such as the counters :data:`A_row <pfnet.ConstraintBase.A_row>`, :data:`A_nnz <pfnet.ConstraintBase.A_nnz>`, :data:`G_row <pfnet.ConstraintBase.G_row>`, and :data:`G_nnz <pfnet.ConstraintBase.G_nnz>`. If used, the array of flags :class:`bus_counted <pfnet.ConstraintBase.bus_counted>`, which has size equal to the number of buses times the number of time periods and can be used to keep track of which buses have already being processed, should also be reset here. 

* :func:`allocate(self) <pfnet.CustomConstraint.allocate>`

  * This method allocates the vectors :data:`b <pfnet.ConstraintBase.b>`, :data:`l <pfnet.ConstraintBase.l>`, and :data:`u <pfnet.ConstraintBase.u>`, and matrices :data:`A <pfnet.ConstraintBase.A>` and :data:`G <pfnet.ConstraintBase.G>` using the methods :func:`set_b() <pfnet.ConstraintBase.set_b>`, :func:`set_l() <pfnet.ConstraintBase.set_l>`, :func:`set_u() <pfnet.ConstraintBase.set_u>`, :func:`set_A() <pfnet.ConstraintBase.set_A>`, and :func:`set_G() <pfnet.ConstraintBase.set_G>`, respectively.

* :func:`analyze_step(self,branch,t) <pfnet.CustomConstraint.analyze_step>`

  * This method is called for every :class:`network branch <pfnet.Branch>` and time period combination, and is responsible for storing the structural or constant information of the matrices :data:`A <pfnet.ConstraintBase.A>` and :data:`G <pfnet.ConstraintBase.G>`. 

* :func:`eval_step(self,branch,t,x) <pfnet.CustomConstraint.eval_step>`

  * This method is used for evaluating nonlinear constraints and will be documented in the near future.

* :func:`store_sens_step(self,branch,t,sA,sf,sGu,sGl) <pfnet.CustomConstraint.store_sens_step>`

  * This method is used for storing constraint sensitivity information in the :ref:`network components <net_components>` and will be documented in the near future.

A template for creating a custom constraint is provided below:

.. literalinclude:: ../examples/custom_constraint_template.py

An example of a custom constraint that constructs the DC power balance equations can be found in `here <https://github.com/ttinoco/PFNET/blob/master/python/pfnet/constraints/dummy_constraint.py>`_. 
