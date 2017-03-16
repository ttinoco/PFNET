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
cimport cfunc

class FunctionError(Exception):
    """
    Function error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Function:
    """
    Function class.
    """

    cdef cfunc.Func* _c_func
    cdef bint alloc

    def __init__(self,weight,Network net,alloc=True):
        """
        Function class.

        Parameters
        ----------
        weight : float
        net : :class:`Network <pfnet.Network>`
        alloc : {``True``, ``False``}
        """

        pass

    def __cinit__(self,weight,Network net,alloc=True):

        if alloc:
            self._c_func = cfunc.FUNC_new(weight,net._c_net)
        else:
            self._c_func = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees function C data structure.
        """

        if self.alloc:
            cfunc.FUNC_del(self._c_func)
            self._c_func = NULL

    def del_matvec(self):
        """
        Deletes matrices and vectors associated with
        this function.
        """

        cfunc.FUNC_del_matvec(self._c_func)

    def update_network(self):
        """
        Updates internal arrays to be compatible
        with any network changes.
        """

        cfunc.FUNC_update_network(self._c_func)

    def clear_error(self):
        """
        Clears error flag and string.
        """

        cfunc.FUNC_clear_error(self._c_func)

    def analyze(self):
        """
        Analyzes function and allocates required vectors and matrices.
        """

        cfunc.FUNC_del_matvec(self._c_func)
        cfunc.FUNC_count(self._c_func)
        cfunc.FUNC_allocate(self._c_func)
        cfunc.FUNC_analyze(self._c_func)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    def eval(self,var_values):
        """
        Evaluates function value, gradient, and Hessian using
        the given variable values.

        Parameters
        ----------
        var_values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if var_values.size else NULL
        cfunc.FUNC_eval(self._c_func,v)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    property name:
        """ Function name (string). """
        def __get__(self): return cfunc.FUNC_get_name(self._c_func).decode('UTF-8')

    property Hcounter:
        """ Number of nonzero entries in Hessian matrix (int). """
        def __get__(self): return cfunc.FUNC_get_Hcounter(self._c_func)

    property phi:
        """ Function value (float). """
        def __get__(self): return cfunc.FUNC_get_phi(self._c_func)

    property gphi:
        """ Function gradient vector (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cfunc.FUNC_get_gphi(self._c_func))

    property Hphi:
        """ Function Hessian matrix (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cfunc.FUNC_get_Hphi(self._c_func))

    property weight:
        """ Function weight (float). """
        def __get__(self): return cfunc.FUNC_get_weight(self._c_func)

cdef new_Function(cfunc.Func* f, cnet.Net* n):
    if f is not NULL and n is not NULL:
        func = Function(0,new_Network(n),alloc=False)
        func._c_func = f
        return func
    else:
        raise FunctionError('invalid function data')

cdef class CustomFunction(Function):
    
    def count_step(self,branch,t):
        pass
        
    def allocate(self):
        pass

    def clear(self):
        pass

    def analyze_step(self,branch,t):
        pass

cdef void func_init(cfunc.Func* f):
    pass # nothing

cdef void func_count_step(cfunc.Func* f, cbranch.Branch* br, int t):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.count_step(new_Branch(br),t)
    # HERE: update fc._c_func

cdef void func_allocate(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.allocate()
    # HERE: update fc._c_func
        
cdef void func_clear(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.clear()
    # HERE: update fc._c_func

cdef void func_analyze_step(cfunc.Func* f, cbranch.Branch* br, int t):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.analyze_step(new_Branch(br),t)
    # HERE: update fc._c_func

cdef void func_eval_step(cfunc.Func* f, cbranch.Branch* br, int t, cvec.Vec* v):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.eval_step(new_Branch(br),t,Vector(v))
    # HERE: update fc._c_func

cdef void func_free(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.free()
    # HERE: update fc._c_func


