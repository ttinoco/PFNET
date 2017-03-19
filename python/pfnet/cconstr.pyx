#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cvec
cimport cbranch
cimport cconstr

class ConstraintError(Exception):
    """
    Constraint error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class ConstraintBase:
    """
    Base constraint class.
    """

    cdef cconstr.Constr* _c_constr
    cdef bint _alloc

    def __init__(self):
        """
        Base constraint class.
        """

        pass

    def __cinit__(self):

        self._c_constr = NULL
        self._alloc = False

    def __dealloc__(self):
        """
        Frees constraint C data structure.
        """

        if self._alloc:
            cconstr.CONSTR_del(self._c_constr)
            self._c_constr = NULL
            self._alloc = False

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
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(x.data),x.size)
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
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(x.data),x.size)
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
        cdef cvec.Vec* vA = cvec.VEC_new_from_array(<cconstr.REAL*>(xA.data),xA.size) if sA is not None else NULL
        cdef cvec.Vec* vf = cvec.VEC_new_from_array(<cconstr.REAL*>(xf.data),xf.size) if sf is not None else NULL
        cdef cvec.Vec* vGu = cvec.VEC_new_from_array(<cconstr.REAL*>(xGu.data),xGu.size) if sGu is not None else NULL
        cdef cvec.Vec* vGl = cvec.VEC_new_from_array(<cconstr.REAL*>(xGl.data),xGl.size) if sGl is not None else NULL
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

    def set_b(self,b):
        """
        Sets b vector.

        Parameters
        ----------
        b : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] bb = b
        PyArray_CLEARFLAGS(bb,np.NPY_OWNDATA)
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(bb.data),bb.size)
        cconstr.CONSTR_set_b(self._c_constr,v)

    def set_A(self,A):
        """
        Sets A matrix.

        Parameters
        ----------
        A : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """
        
        cdef np.ndarray[int,mode='c'] row = A.row
        cdef np.ndarray[int,mode='c'] col = A.col
        cdef np.ndarray[double,mode='c'] data = A.data
        PyArray_CLEARFLAGS(row,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(col,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(data,np.NPY_OWNDATA)
        cdef cmat.Mat* m = cmat.MAT_new_from_arrays(A.shape[0],A.shape[1],A.nnz, 
                                                    <int*>(row.data),
                                                    <int*>(col.data), 
                                                    <cconstr.REAL*>(data.data))
        cconstr.CONSTR_set_A(self._c_constr,m)

    def set_l(self,l):
        """
        Sets l vector.

        Parameters
        ----------
        l : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] ll = l
        PyArray_CLEARFLAGS(ll,np.NPY_OWNDATA)
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(ll.data),ll.size)
        cconstr.CONSTR_set_l(self._c_constr,v)  

    def set_u(self,u):
        """
        Sets u vector.

        Parameters
        ----------
        u : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] uu = u
        PyArray_CLEARFLAGS(uu,np.NPY_OWNDATA)
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(uu.data),uu.size)
        cconstr.CONSTR_set_u(self._c_constr,v)  

    def set_G(self,G):
        """
        Sets G matrix.

        Parameters
        ----------
        G : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """
        
        cdef np.ndarray[int,mode='c'] row = G.row
        cdef np.ndarray[int,mode='c'] col = G.col
        cdef np.ndarray[double,mode='c'] data = G.data
        PyArray_CLEARFLAGS(row,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(col,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(data,np.NPY_OWNDATA)
        cdef cmat.Mat* m = cmat.MAT_new_from_arrays(G.shape[0],G.shape[1],G.nnz, 
                                                    <int*>(row.data),
                                                    <int*>(col.data), 
                                                    <cconstr.REAL*>(data.data))
        cconstr.CONSTR_set_G(self._c_constr,m)

    def set_f(self,f):
        """
        Sets f vector.

        Parameters
        ----------
        f : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] ff = f
        PyArray_CLEARFLAGS(ff,np.NPY_OWNDATA)
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cconstr.REAL*>(ff.data),ff.size)
        cconstr.CONSTR_set_f(self._c_constr,v)  

    def set_J(self,J):
        """
        Sets J matrix.

        Parameters
        ----------
        J : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """
        
        cdef np.ndarray[int,mode='c'] row = J.row
        cdef np.ndarray[int,mode='c'] col = J.col
        cdef np.ndarray[double,mode='c'] data = J.data
        PyArray_CLEARFLAGS(row,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(col,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(data,np.NPY_OWNDATA)
        cdef cmat.Mat* m = cmat.MAT_new_from_arrays(J.shape[0],J.shape[1],J.nnz, 
                                                    <int*>(row.data),
                                                    <int*>(col.data), 
                                                    <cconstr.REAL*>(data.data))
        cconstr.CONSTR_set_J(self._c_constr,m)

    property name:
        """ Constraint name (string). """
        def __get__(self): return cconstr.CONSTR_get_name(self._c_constr)
        def __set__(self,name): 
            name = name.encode('UTF-8')
            cconstr.CONSTR_set_name(self._c_constr,name)

    property A_nnz:
        """ Number of nonzero entries in the matrix of linear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_A_nnz(self._c_constr)
        def __set__(self,nnz): cconstr.CONSTR_set_A_nnz(self._c_constr,nnz)
        
    property G_nnz:
        """ Number of nonzero entries in the matrix of linear inequality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_G_nnz(self._c_constr)
        def __set__(self,nnz): cconstr.CONSTR_set_G_nnz(self._c_constr,nnz)

    property J_nnz:
        """ Number of nonzero entries in the Jacobian matrix of the nonlinear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_J_nnz(self._c_constr)
        def __set__(self,nnz): cconstr.CONSTR_set_J_nnz(self._c_constr,nnz)

    property A_row:
        """ Number of linear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_A_row(self._c_constr)
        def __set__(self,row): cconstr.CONSTR_set_A_row(self._c_constr,row)

    property G_row:
        """ Number of linear ineqquality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_G_row(self._c_constr)
        def __set__(self,row): cconstr.CONSTR_set_G_row(self._c_constr,row)

    property J_row:
        """ Number of nonlinear equality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_J_row(self._c_constr)
        def __set__(self,row): cconstr.CONSTR_set_J_row(self._c_constr,row)

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

    property network:
        """ Network associated with constraint. """
        def __get__(self): return new_Network(cconstr.CONSTR_get_network(self._c_constr))

    property bus_counted:
        """ Boolean array of flags for processing buses during count/analyze/eval, etc. """
        def __get__(self): return BoolArray(cconstr.CONSTR_get_bus_counted(self._c_constr),
                                            cconstr.CONSTR_get_bus_counted_size(self._c_constr))

cdef new_Constraint(cconstr.Constr* c):
    if c is not NULL:
        constr = ConstraintBase()
        constr._c_constr = c
        return constr
    else:
        raise ConstraintError('invalid constraint data')

cdef class Constraint(ConstraintBase):
    
    def __init__(self,name,Network net):
        """
        Function class.
        
        Parameters
        ----------
        name : string
        net : :class:`Network <pfnet.Network>`
        """
        
        pass
    
    def __cinit__(self,name,Network net):

        if name == "AC power balance":
            self._c_constr = cconstr.CONSTR_ACPF_new(net._c_net)
        elif name == "DC power balance":
            self._c_constr = cconstr.CONSTR_DCPF_new(net._c_net)
        elif name == "linearized AC power balance":
            self._c_constr = cconstr.CONSTR_LINPF_new(net._c_net)
        elif name == "variable fixing":
            self._c_constr = cconstr.CONSTR_FIX_new(net._c_net)
        elif name == "variable nonlinear bounds":
            self._c_constr = cconstr.CONSTR_NBOUND_new(net._c_net)
        elif name == "variable bounds":
            self._c_constr = cconstr.CONSTR_LBOUND_new(net._c_net)
        elif name == "generator active power participation":
            self._c_constr = cconstr.CONSTR_PAR_GEN_P_new(net._c_net)
        elif name == "generator reactive power participation":
            self._c_constr = cconstr.CONSTR_PAR_GEN_Q_new(net._c_net)
        elif name == "voltage regulation by generators":
            self._c_constr = cconstr.CONSTR_REG_GEN_new(net._c_net)
        elif name == "voltage regulation by transformers":
            self._c_constr = cconstr.CONSTR_REG_TRAN_new(net._c_net)
        elif name == "voltage regulation by shunts":
            self._c_constr = cconstr.CONSTR_REG_SHUNT_new(net._c_net)
        elif name == "DC branch flow limits":
            self._c_constr = cconstr.CONSTR_DC_FLOW_LIM_new(net._c_net)
        elif name == "AC branch flow limits":
            self._c_constr = cconstr.CONSTR_AC_FLOW_LIM_new(net._c_net)
        elif name == "generator ramp limits":
            self._c_constr = cconstr.CONSTR_GEN_RAMP_new(net._c_net)
        else:
            ConstraintError('invalid constraint name')
            
        self._alloc = True

cdef class CustomConstraint(ConstraintBase):

    def __init__(self,Network net):
        """
        Custom constraint class.
        
        Parameters
        ----------
        net : :class:`Network <pfnet.Network>`
        """

        pass

    def __cinit__(self,Network net):

        self._c_constr = cconstr.CONSTR_new(net._c_net)
        cconstr.CONSTR_set_data(self._c_constr,<void*>self)
        cconstr.CONSTR_set_func_init(self._c_constr,constr_init)
        cconstr.CONSTR_set_func_count_step(self._c_constr,constr_count_step)
        cconstr.CONSTR_set_func_allocate(self._c_constr,constr_allocate)
        cconstr.CONSTR_set_func_clear(self._c_constr,constr_clear)
        cconstr.CONSTR_set_func_analyze_step(self._c_constr,constr_analyze_step)
        cconstr.CONSTR_set_func_eval_step(self._c_constr,constr_eval_step)
        cconstr.CONSTR_set_func_store_sens_step(self._c_constr,constr_store_sens_step)
        self._alloc = True
        self.init()
    
    def init(self):
        """"
        Performs constraint initialization.
        """

        pass
        
    def count_step(self,branch,t):
        """
        Performs count step.

        Parameters
        ----------
        branch : Branch
        t : time period (int)
        """
        
        pass
        
    def allocate(self):
        """
        Allocates matrices and vectors.
        """

        pass

    def clear(self):
        """
        Clears counters and values.
        """

        pass

    def analyze_step(self,branch,t):
        """
        Performs analyze step.
       
        Parameters
        ----------
        branch : Branch
        t : time period (int)
        """
        
        pass

    def eval_step(self,branch,t,x):
        """
        Performs eval step.
       
        Parameters
        ----------
        branch : Branch
        t : time period (int)
        x : ndarray
        """
 
        pass

    def store_sens_step(self,branch,t,sA,sf,sGu,sGl):
        """
        Performs step for storing sensitivities.
       
        Parameters
        ----------
        branch : Branch
        t : time period (int)
        sA : ndarray
        sf : ndarray
        sGu : ndarray
        sGl : ndarray
        """
 
        pass

cdef void constr_init(cconstr.Constr* c):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.init()

cdef void constr_count_step(cconstr.Constr* c, cbranch.Branch* br, int t):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.count_step(new_Branch(br),t)

cdef void constr_allocate(cconstr.Constr* c):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.allocate()
        
cdef void constr_clear(cconstr.Constr* c):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.clear()

cdef void constr_analyze_step(cconstr.Constr* c, cbranch.Branch* br, int t):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.analyze_step(new_Branch(br),t)

cdef void constr_eval_step(cconstr.Constr* c, cbranch.Branch* br, int t, cvec.Vec* v):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.eval_step(new_Branch(br),t,Vector(v))

cdef void constr_store_sens_step(cconstr.Constr* c, cbranch.Branch* br, int t, cvec.Vec* sA, cvec.Vec* sf, cvec.Vec* sGu, cvec.Vec* sGl):
    cdef CustomConstraint cc = <CustomConstraint>cconstr.CONSTR_get_data(c)
    cc.eval_step(new_Branch(br),t,Vector(sA),Vector(sf),Vector(sGu),Vector(sGl))
