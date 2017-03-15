#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cconstr

class ConstraintError(Exception):
    """
    Constraint error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Constraint:
    """
    Constraint class.
    """

    cdef cconstr.Constr* _c_constr
    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self,ctype,Network net,alloc=True):
        """
        Contraint class.

        Parameters
        ----------
        ctype : string (:ref:`ref_constr_type`)
        net : :class:`Network <pfnet.Network>`
        alloc : {``True``, ``False``}
        """

        pass

    def __cinit__(self,ctype,Network net,alloc=True):

        self._c_net = net._c_net
        if alloc:
            self._c_constr = cconstr.CONSTR_new(str2constr[ctype],net._c_net)
        else:
            self._c_constr = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees constraint C data structure.
        """

        if self.alloc:
            cconstr.CONSTR_del(self._c_constr)
            self._c_constr = NULL

    def del_matvec(self):
        """
        Deletes matrices and vectors associated with
        this constraint.
        """

        cconstr.CONSTR_del_matvec(self._c_constr)

    def update_network(self):
        """
        Updates internal arrays to be compatible
        with any network changes.
        """

        cconstr.CONSTR_update_network(self._c_constr)

    def clear_error(self):
        """
        Clears error flag and string.
        """

        cconstr.CONSTR_clear_error(self._c_constr)

    def analyze(self):
        """
        Analyzes constraint and allocates required vectors and matrices.
        """

        cconstr.CONSTR_del_matvec(self._c_constr)
        cconstr.CONSTR_count(self._c_constr)
        cconstr.CONSTR_allocate(self._c_constr)
        cconstr.CONSTR_analyze(self._c_constr)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def combine_H(self,coeff,ensure_psd=False):
        """
        Forms and saves a linear combination of the individual constraint Hessians.

        Parameters
        ----------
        coeff : :class:`ndarray <numpy.ndarray>`
        ensure_psd : {``True``, ``False``}
        """

        cdef np.ndarray[double,mode='c'] x = coeff
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if coeff.size else NULL
        cconstr.CONSTR_combine_H(self._c_constr,v,ensure_psd)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def eval(self,values):
        """
        Evaluates constraint violations, Jacobian, and individual Hessian matrices.

        Parameters
        ----------
        values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),x.size) if values.size else NULL
        cconstr.CONSTR_eval(self._c_constr,v)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

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
        cdef cvec.Vec* vA = cvec.VEC_new_from_array(&(xA[0]),len(xA)) if (sA is not None and sA.size) else NULL
        cdef cvec.Vec* vf = cvec.VEC_new_from_array(&(xf[0]),len(xf)) if (sf is not None and sf.size) else NULL
        cdef cvec.Vec* vGu = cvec.VEC_new_from_array(&(xGu[0]),len(xGu)) if (sGu is not None and sGu.size) else NULL
        cdef cvec.Vec* vGl = cvec.VEC_new_from_array(&(xGl[0]),len(xGl)) if (sGl is not None and sGl.size) else NULL
        cconstr.CONSTR_store_sens(self._c_constr,vA,vf,vGu,vGl)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def get_H_single(self,i):
        """
        Gets the Hessian matrix (only lower triangular part) of an individual constraint.

        Parameters
        ----------
        i : int

        Returns
        -------
        H : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """
        return Matrix(cconstr.CONSTR_get_H_single(self._c_constr,i))

    property type:
        """ Constraint type (string) (:ref:`ref_constr_type`). """
        def __get__(self): return constr2str[cconstr.CONSTR_get_type(self._c_constr)]

    property A_nnz:
        """ Number of nonzero entries in the matrix of linear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_A_nnz(self._c_constr)

    property G_nnz:
        """ Number of nonzero entries in the matrix of linear inequality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_G_nnz(self._c_constr)

    property J_nnz:
        """ Number of nonzero entries in the Jacobian matrix of the nonlinear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_J_nnz(self._c_constr)

    property A_row:
        """ Number of linear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_A_row(self._c_constr)

    property G_row:
        """ Number of linear ineqquality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_G_row(self._c_constr)

    property J_row:
        """ Number of nonlinear equality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_J_row(self._c_constr)

    property f:
        """ Vector of nonlinear equality constraint violations (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_f(self._c_constr))

    property J:
        """ Jacobian matrix of nonlinear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_J(self._c_constr))

    property Jbar:
        """ Jacobian matrix of nonlinear equality constraints wrt. extra variables (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_Jbar(self._c_constr))

    property b:
        """ Right-hand side vector of linear equality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_b(self._c_constr))

    property A:
        """ Matrix for linear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_A(self._c_constr))

    property l:
        """ Lower bound vector of linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_l(self._c_constr))

    property u:
        """ Upper bound vector of linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_u(self._c_constr))

    property G:
        """ Matrix for linear inequality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_G(self._c_constr))

    property Gbar:
        """ Matrix for linear inequality constraints, for extra variables (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_Gbar(self._c_constr))

    property H_combined:
        """ Linear combination of Hessian matrices of individual nonlinear equality constraints (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_H_combined(self._c_constr))

    property num_extra_vars:
            """ Number of extra variables (set during count) (int). """
            def __get__(self): return cconstr.CONSTR_get_num_extra_vars(self._c_constr)

cdef new_Constraint(cconstr.Constr* c, cnet.Net* n):
    if c is not NULL and n is not NULL:
        constr = Constraint(0,new_Network(n),alloc=False)
        constr._c_constr = c
        return constr
    else:
        raise ConstraintError('invalid constraint data')
