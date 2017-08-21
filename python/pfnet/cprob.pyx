#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cprob

class ProblemError(Exception):
    """
    Problem error exception.
    """

    pass

cdef class Problem:
    """
    Optimization problem class.
    """

    cdef cprob.Prob* _c_prob
    cdef bint alloc
    cdef list _functions
    cdef list _constraints

    def __init__(self, Network net):
        """
        Optimization problem class.

        Parameters
        ----------
        net : :class:`Network <pfnet.Network>`
        """

        pass

    def __cinit__(self, Network net):

        self._c_prob = cprob.PROB_new(net._c_net)
        self.alloc = True
        self._functions = []
        self._constraints = []

    def __dealloc__(self):
        """
        Frees problem C data structure.
        """

        if self.alloc:
            cprob.PROB_del(self._c_prob)
            self._c_prob = NULL

    def add_constraint(self,ConstraintBase constr):
        """
        Adds constraint to optimization problem.

        Parameters
        ----------
        constr : ConstraintBase
        """

        # Prevent __dealloc__ of function
        constr._alloc = False
        
        # Add constraint to list
        # Prevents python object from being garbage collected
        # This is needed for CustomConstraints
        self._constraints.append(constr)
        
        # Add constraint to problem    
        cprob.PROB_add_constr(self._c_prob,constr._c_constr)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def add_function(self,FunctionBase func):
        """
        Adds function to optimization problem objective.

        Parameters
        ----------
        func : FunctionBase
        """
        
        # Prevent __dealloc__ of function
        func._alloc = False
        
        # Add function to list
        # Prevents python object from being garbage collected
        # This is needed for CustomFunctions
        self._functions.append(func)

        # Add function to problem
        cprob.PROB_add_func(self._c_prob,func._c_func)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def add_heuristic(self,htype):

        cprob.PROB_add_heur(self._c_prob,htype)

    def analyze(self):
        """
        Analyzes function and constraint structures and allocates
        required vectors and matrices.
        """

        cprob.PROB_analyze(self._c_prob)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def apply_heuristics(self,var_values):
        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cprob.REAL*>(x.data),x.size)
        cprob.PROB_apply_heuristics(self._c_prob,v)
        free(v)

    def clear(self):
        """
        Resets optimization problem data.
        """

        self._functions = []
        cprob.PROB_clear(self._c_prob)

    def clear_error(self):
        """
        Clears error flag and string.
        """

        cprob.PROB_clear_error(self._c_prob)

    def combine_H(self,coeff,ensure_psd=False):
        """
        Forms and saves a linear combination of the individual constraint Hessians.

        Parameters
        ----------
        coeff : :class:`ndarray <numpy.ndarray>`
        ensure_psd : {``True``, ``False``}
        """

        cdef np.ndarray[double,mode='c'] x = coeff
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cprob.REAL*>(x.data),x.size)
        cprob.PROB_combine_H(self._c_prob,v,ensure_psd)
        free(v)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def has_error(self):
        """
        Indicates whether the problem has the error flag set due to an
        invalid operation.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cprob.PROB_has_error(self._c_prob)

    def eval(self,var_values):
        """
        Evaluates objective function and constraints as well as their first and
        second derivatives using the given variable values.

        Parameters
        ----------
        var_values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cprob.REAL*>(x.data),x.size)
        cprob.PROB_eval(self._c_prob,v)
        free(v)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def store_sensitivities(self,sA,sf,sGu,sGl):
        """
        Stores Lagrange multiplier estimates of the constraints in
        the power network components.

        Parameters
        ----------
        sA : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear equality constraints (:math:`Ax = b`)
        sf : :class:`ndarray <numpy.ndarray>`
             sensitivities for nonlinear equality constraints (:math:`f(x) = 0`)
        sGu : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`Gx \le u`)
        sGl : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`l \le Gx`)
        """

        cdef np.ndarray[double,mode='c'] xA = sA
        cdef np.ndarray[double,mode='c'] xf = sf
        cdef np.ndarray[double,mode='c'] xGu = sGu
        cdef np.ndarray[double,mode='c'] xGl = sGl
        cdef cvec.Vec* vA = cvec.VEC_new_from_array(<cprob.REAL*>(xA.data),xA.size) if sA is not None else NULL
        cdef cvec.Vec* vf = cvec.VEC_new_from_array(<cprob.REAL*>(xf.data),xf.size) if sf is not None else NULL
        cdef cvec.Vec* vGu = cvec.VEC_new_from_array(<cprob.REAL*>(xGu.data),xGu.size) if sGu is not None else NULL
        cdef cvec.Vec* vGl = cvec.VEC_new_from_array(<cprob.REAL*>(xGl.data),xGl.size) if sGl is not None else NULL
        cprob.PROB_store_sens(self._c_prob,vA,vf,vGu,vGl)
        if vA != NULL:
            free(vA)
        if vf != NULL:
            free(vf)
        if vGu != NULL:
            free(vGu)
        if vGl != NULL:
            free(vGl)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def find_constraint(self,name):
        """
        Finds constraint of give type among the constraints of this optimization problem.

        Parameters
        ----------
        name : string

        Returns
        -------
        constr : Constraint
        """

        name = name.encode('UTF-8')
        c = cprob.PROB_find_constr(self._c_prob,name)
        if c is not NULL:
            return new_Constraint(c)
        else:
            raise ProblemError('constraint not found')

    def get_init_point(self):
        """
        Gets initial solution estimate from the current values of the network variables.

        Returns
        -------
        point : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_init_point(self._c_prob),owndata=True)

    def get_upper_limits(self):
        """
        Gets vector of upper limits for the network and extra variables.

        Returns
        -------
        limits : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_upper_limits(self._c_prob),owndata=True)

    def get_lower_limits(self):
        """
        Gets vector of lower limits for the network and extra variables.

        Returns
        -------
        limits : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_lower_limits(self._c_prob),owndata=True)

    def get_network(self):
        """
        Gets the power network associated with this optimization problem.
        """

        return new_Network(cprob.PROB_get_network(self._c_prob))

    def show(self):
        """
        Shows information about this optimization problem.
        """

        print(cprob.PROB_get_show_str(self._c_prob).decode('UTF-8'))

    def update_lin(self):
        """
        Updates linear equality constraints.
        """

        cprob.PROB_update_lin(self._c_prob)

    def get_num_primal_variables(self):
        """ 
        Gets number of primal variables. 

        Returns
        -------
        num : int
        """
        
        return self.num_primal_variables

    def get_num_linear_equality_constraints(self):    
        """ 
        Gets number of linear equality constraints.

        Returns
        -------
        num : int
        """
        
        return self.num_linear_equality_constraints

    def get_num_nonlinear_equality_constraints(self):
        """ 
        Number of nonlinear equality constraints.

        Returns
        -------
        num : int
        """
        
        return self.num_nonlinear_equality_constraints

    property network:
        """ Power network associated with this optimization problem (:class:`Network <pfnet.Network>`). """
        def __get__(self): return new_Network(cprob.PROB_get_network(self._c_prob))

    property constraints:
        """ List of :class:`constraints <pfnet.Constraint>` of this optimization problem (list). """
        def __get__(self):
            clist = []
            cdef cconstr.Constr* c = cprob.PROB_get_constr(self._c_prob)
            while c is not NULL:
                clist.append(new_Constraint(c))
                c = cconstr.CONSTR_get_next(c)
            return clist

    property functions:
        """ List of :class:`functions <pfnet.Function>` that form the objective function of this optimization problem (list). """
        def __get__(self):
            flist = []
            cdef cfunc.Func* f = cprob.PROB_get_func(self._c_prob)
            while f is not NULL:
                flist.append(new_Function(f))
                f = cfunc.FUNC_get_next(f)
            return flist

    property A:
        """ Constraint matrix of linear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_A(self._c_prob))

    property b:
        """ Right hand side vectors of the linear equality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_b(self._c_prob))

    property G:
        """ Constraint matrix of linear inequality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_G(self._c_prob))

    property l:
        """ Lower bound for linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_l(self._c_prob))

    property u:
        """ Upper bound for linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_u(self._c_prob))

    property J:
        """ Jacobian matrix of the nonlinear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_J(self._c_prob))

    property f:
        """ Vector of nonlinear equality constraints violations (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_f(self._c_prob))

    property phi:
        """ Objective function value (float). """
        def __get__(self): return cprob.PROB_get_phi(self._c_prob)

    property gphi:
        """ Objective function gradient vector (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_gphi(self._c_prob))

    property Hphi:
        """ Objective function Hessian matrix (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_Hphi(self._c_prob))

    property H_combined:
        """ Linear combination of Hessian matrices of individual nonlinear equality constraints (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_H_combined(self._c_prob))

    property x:
        """ Initial primal point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return self.get_init_point()

    property lam:
        """ Initial dual point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return None

    property nu:
        """ Initial dual point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return None

    property num_primal_variables:
        """ Number of primal variables (int). """
        def __get__(self): return cprob.PROB_get_num_primal_variables(self._c_prob)

    property num_linear_equality_constraints:    
        """ Number of linear equality constraints (int). """
        def __get__(self): return cprob.PROB_get_num_linear_equality_constraints(self._c_prob)

    property num_nonlinear_equality_constraints:
        """ Number of nonlinear equality constraints (int). """
        def __get__(self): return cprob.PROB_get_num_nonlinear_equality_constraints(self._c_prob)

    property num_extra_vars:
        """ Number of extra varaibles (set during analyze) (int). """
        def __get__(self): return cprob.PROB_get_num_extra_vars(self._c_prob)
