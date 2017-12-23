#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
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

    pass

cdef class FunctionBase:
    """
    Base function class.
    """

    cdef cfunc.Func* _c_func
    cdef bint _alloc

    def __init__(self):
        """
        Base function class.
        """

        pass

    def __cinit__(self):

        self._c_func = NULL
        self._alloc = False

    def __dealloc__(self):
        """
        Frees function C data structure.
        """

        if self._alloc:
            cfunc.FUNC_del(self._c_func)
            self._c_func = NULL
            self._alloc = False
            
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
        Analyzes the structure of the function and allocates required vectors and matrices.
        """

        cfunc.FUNC_del_matvec(self._c_func)
        cfunc.FUNC_count(self._c_func)
        cfunc.FUNC_allocate(self._c_func)
        cfunc.FUNC_analyze(self._c_func)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    def eval(self, values):
        """
        Evaluates function value, gradient, and Hessian using
        the given variable values.

        Parameters
        ----------
        values : |Array|
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cfunc.REAL*>(x.data),x.size)
        cfunc.FUNC_eval(self._c_func,v)
        free(v)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    def set_gphi(self, gphi):
        """
        Sets gradient vector.

        Parameters
        ----------
        gphi : |Array|
        """
        
        cdef np.ndarray[double,mode='c'] g = gphi
        PyArray_CLEARFLAGS(g,np.NPY_OWNDATA)
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cfunc.REAL*>(g.data),g.size)
        cfunc.FUNC_set_gphi(self._c_func,v)

    def set_Hphi(self, Hphi):
        """
        Sets Hessian matrix.

        Parameters
        ----------
        Hphi : |CooMatrix| (lower triangular)
        """
        
        cdef np.ndarray[int,mode='c'] row = Hphi.row
        cdef np.ndarray[int,mode='c'] col = Hphi.col
        cdef np.ndarray[double,mode='c'] data = Hphi.data
        PyArray_CLEARFLAGS(row,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(col,np.NPY_OWNDATA)
        PyArray_CLEARFLAGS(data,np.NPY_OWNDATA)
        cdef cmat.Mat* m = cmat.MAT_new_from_arrays(Hphi.shape[0],Hphi.shape[1],Hphi.nnz, 
                                                    <int*>(row.data),
                                                    <int*>(col.data), 
                                                    <cfunc.REAL*>(data.data))
        cfunc.FUNC_set_Hphi(self._c_func,m)

    def set_parameter(self, key, value):
        """
        Sets function parameter.

        Parameters
        ----------
        key : string
        value : object
        """

        cdef void* cvalue = NULL
        cdef np.ndarray[double, mode='c'] a
        key = key.encode('UTF-8')

        # int

        # float

        # ndarray
        if issubclass(type(value), np.ndarray):
            a = value
            cvalue = <void*>a.data

        # Unknown
        if cvalue == NULL:
            raise FunctionError('invalid value type')

        # Set parameter
        cfunc.FUNC_set_parameter(self._c_func, key, cvalue)
        
        # Error
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    property name:
        """ Function name (string). """
        def __get__(self): return cfunc.FUNC_get_name(self._c_func).decode('UTF-8')
        def __set__(self,name): 
            name = name.encode('UTF-8')
            cfunc.FUNC_set_name(self._c_func,name)

    property phi:
        """ Function value (float). """
        def __get__(self): return cfunc.FUNC_get_phi(self._c_func)
        def __set__(self,phi): cfunc.FUNC_set_phi(self._c_func,phi)

    property gphi:
        """ Function gradient vector (|Array|). """
        def __get__(self): return Vector(cfunc.FUNC_get_gphi(self._c_func))

    property Hphi:
        """ Function Hessian matrix (only the lower triangular part) (|CooMatrix|). """
        def __get__(self): return Matrix(cfunc.FUNC_get_Hphi(self._c_func))

    property Hphi_nnz:
        """ Counter of nonzero entries in Hessian matrix (int). """
        def __get__(self): return cfunc.FUNC_get_Hphi_nnz(self._c_func)
        def __set__(self,nnz): cfunc.FUNC_set_Hphi_nnz(self._c_func,nnz);

    property weight:
        """ Function weight (float). """
        def __get__(self): return cfunc.FUNC_get_weight(self._c_func)

    property network:
        """ Network associated with function (|Network|). """
        def __get__(self): return new_Network(cfunc.FUNC_get_network(self._c_func))

    property bus_counted:
        """ Boolean array of flags for processing buses during count/analyze/eval, etc (|Array|). """
        def __get__(self): return BoolArray(cfunc.FUNC_get_bus_counted(self._c_func),
                                            cfunc.FUNC_get_bus_counted_size(self._c_func))

cdef new_Function(cfunc.Func* f):
    if f is not NULL:
        func = FunctionBase()
        func._c_func = f
        return func
    else:
        raise FunctionError('invalid function data')

cdef class Function(FunctionBase):
    
    def __init__(self, name, weight, Network net):
        """
        Function class.
        
        Parameters
        ----------
        name : string
        weight : float
        net : |Network|
        """
        
        pass
    
    def __cinit__(self, name, weight, Network net):
                
        if name == "generation cost":
            self._c_func = cfunc.FUNC_GEN_COST_new(weight, net._c_net)
        elif name == "consumption utility":
            self._c_func = cfunc.FUNC_LOAD_UTIL_new(weight, net._c_net)
        elif name == "net consumption cost":
            self._c_func = cfunc.FUNC_NETCON_COST_new(weight, net._c_net)
        elif name == "phase shift regularization":
            self._c_func = cfunc.FUNC_REG_PHASE_new(weight, net._c_net)
        elif name == "generator powers regularization":
            self._c_func = cfunc.FUNC_REG_PQ_new(weight, net._c_net)
        elif name == "tap ratio regularization":
            self._c_func = cfunc.FUNC_REG_RATIO_new(weight, net._c_net)
        elif name == "susceptance regularization":
            self._c_func = cfunc.FUNC_REG_SUSC_new(weight, net._c_net)
        elif name == "voltage angle regularization":
            self._c_func = cfunc.FUNC_REG_VANG_new(weight, net._c_net)
        elif name == "voltage magnitude regularization":
            self._c_func = cfunc.FUNC_REG_VMAG_new(weight, net._c_net)
        elif name == "variable regularization":
            self._c_func = cfunc.FUNC_REG_VAR_new(weight, net._c_net)
        elif name == "soft voltage magnitude limits":
            self._c_func = cfunc.FUNC_SLIM_VMAG_new(weight, net._c_net)
        elif name == "sparse controls penalty":
            self._c_func = cfunc.FUNC_SP_CONTROLS_new(weight, net._c_net)
        else:
            raise FunctionError('invalid function name')
            
        self._alloc = True
    
cdef class CustomFunction(FunctionBase):
    """
    Custom function class.
    """

    def __init__(self, weight, Network net):
        """
        Custom function class.
        
        Parameters
        ----------
        weight : float
        net : |Network|
        """
        
        pass

    def __cinit__(self, weight, Network net):
        
        self._c_func = cfunc.FUNC_new(weight,net._c_net)
        cfunc.FUNC_set_data(self._c_func,<void*>self)
        cfunc.FUNC_set_func_init(self._c_func,func_init)
        cfunc.FUNC_set_func_count_step(self._c_func,func_count_step)
        cfunc.FUNC_set_func_allocate(self._c_func,func_allocate)
        cfunc.FUNC_set_func_clear(self._c_func,func_clear)
        cfunc.FUNC_set_func_analyze_step(self._c_func,func_analyze_step)
        cfunc.FUNC_set_func_eval_step(self._c_func,func_eval_step)
        cfunc.FUNC_init(self._c_func)
        self._alloc = True

    def init(self):
        """
        Performs function initialization.
        """

        pass
        
    def count_step(self,branch,t):
        """
        Performs count step.

        Parameters
        ----------
        branch : |Branch|
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

    def analyze_step(self, branch, t):
        """
        Performs analyze step.
       
        Parameters
        ----------
        branch : |Branch|
        t : time period (int)
        """
        
        pass

    def eval_step(self, branch, t, x):
        """
        Performs eval step.
       
        Parameters
        ----------
        branch : |Branch|
        t : time period (int)
        x : |Array|
        """
 
        pass

cdef void func_init(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.init()

cdef void func_count_step(cfunc.Func* f, cbranch.Branch* br, int t):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.count_step(new_Branch(br),t)

cdef void func_allocate(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.allocate()
        
cdef void func_clear(cfunc.Func* f):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.clear()

cdef void func_analyze_step(cfunc.Func* f, cbranch.Branch* br, int t):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.analyze_step(new_Branch(br),t)

cdef void func_eval_step(cfunc.Func* f, cbranch.Branch* br, int t, cvec.Vec* v):
    cdef CustomFunction fc = <CustomFunction>cfunc.FUNC_get_data(f)
    fc.eval_step(new_Branch(br),t,Vector(v))


