#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cload

# Infinity
LOAD_INF_P = cload.LOAD_INF_P
LOAD_INF_Q = cload.LOAD_INF_Q

# Others
LOAD_MIN_TARGET_PF = cload.LOAD_MIN_TARGET_PF

class LoadError(Exception):
    """
    Load error exception.
    """

    pass

cdef class Load:
    """
    Load class.
    """

    cdef cload.Load* _c_ptr
    cdef bint alloc

    def __init__(self, num_periods=1, alloc=True):
        """
        Load class.

        Parameters
        ----------
        num_periods : int
        alloc : |TrueFalse|
        """

        pass

    def __cinit__(self, num_periods=1, alloc=True):

        if alloc:
            self._c_ptr = cload.LOAD_new(num_periods)
        else:
            self._c_ptr = NULL
        self.alloc = alloc

    def __dealloc__(self):

        if self.alloc:
            cload.LOAD_array_del(self._c_ptr,1)
            self._c_ptr = NULL    

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_equal(self, other):
        """
        Determines whether load is equal to given load.

        Parameters
        ----------
        other : |Load|
        """

        cdef Load l_other

        if not isinstance(other,Load):
            return False

        l_other = other

        return cload.LOAD_is_equal(self._c_ptr, l_other._c_ptr)

    def is_P_adjustable(self):
        """
        Determines whether the load has adjustable active power.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cload.LOAD_is_P_adjustable(self._c_ptr)

    def has_flags(self, flag_type, q):
        """
        Determines whether the load has the flags associated with
        certain quantities set.

        Parameters
        ----------
        flag_type : string (|RefFlags|)
        q : string or list of strings (|RefLoadQuantities|)

        Returns
        -------
        flag : |TrueFalse|
        """

        q = q if isinstance(q,list) else [q]

        return cload.LOAD_has_flags(self._c_ptr,
                                    str2flag[flag_type],
                                    reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0))

    def get_var_info_string(self, index):
        """
        Gets info string of variable associated with index.

        Parameters
        ----------
        index : int

        Returns
        -------
        info : string
        """

        cdef char* info_string = cload.LOAD_get_var_info_string(self._c_ptr, index)
        if info_string:
            s = info_string.decode('UTF-8')
            free(info_string)
            return s
        else:
            raise LoadError('index does not correspond to any variable')

    def set_P(self, P, t=0):
        """
        Sets active power.

        Parameters
        ----------
        P : float
        t : int
        """

        cload.LOAD_set_P(self._c_ptr,P,t)

    def set_P_max(self, P, t=0):
        """
        Sets active power upper limit.

        Parameters
        ----------
        P : float
        t : int
        """

        cload.LOAD_set_P_max(self._c_ptr,P,t)

    def set_P_min(self, P, t=0):
        """
        Sets active power lower limit.

        Parameters
        ----------
        P : float
        t : int
        """

        cload.LOAD_set_P_min(self._c_ptr,P,t)

    def set_Q(self, Q, t=0):
        """
        Sets reactive power.

        Parameters
        ----------
        Q : float
        t : int
        """

        cload.LOAD_set_Q(self._c_ptr,Q,t)

    property name:
        """ Load name (string). """
        def __get__(self):
            return cload.LOAD_get_name(self._c_ptr).decode('UTF-8')
        def __set__(self,name):
            name = name.encode('UTF-8')
            cload.LOAD_set_name(self._c_ptr,name)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cload.LOAD_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (string). """
        def __get__(self): return obj2str[cload.LOAD_get_obj_type(self._c_ptr)]

    property index:
        """ Load index (int). """
        def __get__(self): return cload.LOAD_get_index(self._c_ptr)

    property index_P:
        """ Index of load active power variable (int or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_Q:
        """ Index of load reactive power variable (int or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_index_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property bus:
        """ |Bus| to which load is connected. """
        def __get__(self): 
            return new_Bus(cload.LOAD_get_bus(self._c_ptr))
        def __set__(self,bus):
            cdef Bus cbus
            if not isinstance(bus,Bus) and bus is not None:
                raise LoadError('Not a Bus type object')
            cbus = bus
            cload.LOAD_set_bus(self._c_ptr,cbus._c_ptr if bus is not None else NULL)

    property P:
        """ Load active power (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_P)
        def __set__(self,P):
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cload.LOAD_set_P(self._c_ptr,Par[t],t)

    property P_max:
        """ Load active power upper limit (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_P_max(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_P_max)
        def __set__(self,P):
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cload.LOAD_set_P_max(self._c_ptr,Par[t],t)

    property P_min:
        """ Load active power lower limit (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_P_min(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_P_min)
        def __set__(self,P):
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cload.LOAD_set_P_min(self._c_ptr,Par[t],t)

    property Q:
        """ Load reactive power (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_Q)
        def __set__(self,Q):
            cdef int t
            cdef np.ndarray Qar = np.array(Q).flatten()
            for t in range(np.minimum(Qar.size,self.num_periods)):
                cload.LOAD_set_Q(self._c_ptr,Qar[t],t)

    property P_util:
        """ Active power load utility ($/hr) (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_P_util(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property power_factor:
        """ Load power factor (float or |Array|). """
        def __get__(self):
            r = [cload.LOAD_get_power_factor(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property target_power_factor:
        """ Target load power factor in (0,1] (float). """
        def __get__(self): return cload.LOAD_get_target_power_factor(self._c_ptr)
        def __set__(self,pf): cload.LOAD_set_target_power_factor(self._c_ptr,pf)

    property util_coeff_Q0:
        """ Coefficient for consumption utility function (constant term, units of $/hr) (float). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q0(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q0(self._c_ptr,c)

    property util_coeff_Q1:
        """ Coefficient for consumption utility function (linear term, units of $/(hr p.u.)) (float). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q1(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q1(self._c_ptr,c)

    property util_coeff_Q2:
        """ Coefficient for consumption utility function (quadratic term, units of $/(hr p.u.^2)) (float). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q2(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q2(self._c_ptr,c)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = cload.LOAD_get_json_string(self._c_ptr, NULL)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power upper bound (float or |Array|). """
        def __get__(self): return DoubleArray(cload.LOAD_get_sens_P_u_bound_array(self._c_ptr),
                                              cload.LOAD_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_P_u_bound[:] = x

    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power lower bound (float or |Array|). """
        def __get__(self): return DoubleArray(cload.LOAD_get_sens_P_l_bound_array(self._c_ptr),
                                              cload.LOAD_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_P_l_bound[:] = x

    property flags_vars:
        """ Flags associated with variable quantities (byte). """
        def __get__(self): return cload.LOAD_get_flags_vars(self._c_ptr)

    property flags_fixed:
        """ Flags associated with fixed quantities (byte). """
        def __get__(self): return cload.LOAD_get_flags_fixed(self._c_ptr)

    property flags_bounded:
        """ Flags associated with bounded quantities (byte). """
        def __get__(self): return cload.LOAD_get_flags_bounded(self._c_ptr)

    property flags_sparse:
        """ Flags associated with sparse quantities (byte). """
        def __get__(self): return cload.LOAD_get_flags_sparse(self._c_ptr)
            
cdef new_Load(cload.Load* l):
    if l is not NULL:
        load = Load(alloc=False)
        load._c_ptr = l
        return load
    else:
        raise LoadError('no load data')
