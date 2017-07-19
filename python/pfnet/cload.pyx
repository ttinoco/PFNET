#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
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

    def __init__(self,num_periods=1,alloc=True):
        """
        Load class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cload.LOAD_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_P_adjustable(self):
        """
        Determines whether the load has adjustable active power.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cload.LOAD_is_P_adjustable(self._c_ptr)

    def has_flags(self,flag_type,q):
        """
        Determines whether the load has the flags associated with
        certain quantities set.

        Parameters
        ----------
        flag_type : string (:ref:`ref_net_flag`)
        q : string or list of strings (:ref:`ref_load_q`)

        Returns
        -------
        flag : {``True``, ``False``}
        """

        q = q if isinstance(q,list) else [q]

        return cload.LOAD_has_flags(self._c_ptr,
                                    str2flag[flag_type],
                                    reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0))

    def set_P(self,P,t=0):
        """"
        Sets active power.

        Parameters
        ----------
        P : float
        t = int
        """

        cload.LOAD_set_P(self._c_ptr,P,t)

    def set_P_max(self,P,t=0):
        """"
        Sets active power upper limit.

        Parameters
        ----------
        P : float
        t = int
        """

        cload.LOAD_set_P_max(self._c_ptr,P,t)

    def set_P_min(self,P,t=0):
        """"
        Sets active power lower limit.

        Parameters
        ----------
        P : float
        t = int
        """

        cload.LOAD_set_P_min(self._c_ptr,P,t)

    def set_Q(self,Q,t=0):
        """"
        Sets reactive power.

        Parameters
        ----------
        Q : float
        t = int
        """

        cload.LOAD_set_Q(self._c_ptr,Q,t)

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
        """ Index of load active power variable (int or array). """
        def __get__(self):
            r = [cload.LOAD_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_Q:
        """ Index of load reactive power variable (int or array). """
        def __get__(self):
            r = [cload.LOAD_get_index_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which load is connected. """
        def __get__(self): return new_Bus(cload.LOAD_get_bus(self._c_ptr))

    property P:
        """ Load active power (p.u. system base MVA) (float or array). """
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
        """ Load active power upper limit (p.u. system base MVA) (float or array). """
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
        """ Load active power lower limit (p.u. system base MVA) (float or array). """
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
        """ Load reactive power (p.u. system base MVA) (float or array). """
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
        """ Active power load utility ($/hr) (float or array). """
        def __get__(self):
            r = [cload.LOAD_get_P_util(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property power_factor:
        """ Load power factor (float or array). """
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
        """ Coefficient for consumption utility function (constant term, units of $/hr). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q0(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q0(self._c_ptr,c)

    property util_coeff_Q1:
        """ Coefficient for consumption utility function (linear term, units of $/(hr p.u.)). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q1(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q1(self._c_ptr,c)

    property util_coeff_Q2:
        """ Coefficient for consumption utility function (quadratic term, units of $/(hr p.u.^2)). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q2(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q2(self._c_ptr,c)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = cload.LOAD_get_json_string(self._c_ptr)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power upper bound (float or array). """
        def __get__(self):
            r = [cload.LOAD_get_sens_P_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power lower bound (float or array). """
        def __get__(self):
            r = [cload.LOAD_get_sens_P_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

cdef new_Load(cload.Load* l):
    if l is not NULL:
        load = Load(alloc=False)
        load._c_ptr = l
        return load
    else:
        raise LoadError('no load data')
