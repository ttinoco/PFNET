#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cgen

# Infinity
GEN_INF_P = cgen.GEN_INF_P
GEN_INF_Q = cgen.GEN_INF_Q

class GeneratorError(Exception):
    """
    Generator error exception.
    """

    pass

cdef class Generator:
    """
    Generator class.
    """

    cdef cgen.Gen* _c_ptr

    def __init__(self, num_periods=1, alloc=True):
        """
        Generator class.

        Parameters
        ----------
        num_periods : int
        alloc : |TrueFalse|
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cgen.GEN_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_equal(self, other):
        """
        Determines whether generator is equal to given generator.

        Parameters
        ----------
        other : |Generator|
        """

        cdef Generator g_other

        if not isinstance(other,Generator):
            return False

        g_other = other

        return cgen.GEN_is_equal(self._c_ptr,g_other._c_ptr)

    def __richcmp__(self, other, op):
        """
        Compares two generators.

        Parameters
        ----------
        other : |Generator|
        op : comparison type

        Returns
        -------
        flag : {``True``, ``False``}
        """

        if op == 2:
            return self.is_equal(other)
        elif op == 3:
            return not self.is_equal(other)
        else:
            return False

    def is_on_outage(self):
        """
        Determines whether generator in on outage.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cgen.GEN_is_on_outage(self._c_ptr)

    def is_slack(self):
        """
        Determines whether generator is slack.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cgen.GEN_is_slack(self._c_ptr)

    def is_regulator(self):
        """
        Determines whether generator provides voltage regulation.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cgen.GEN_is_regulator(self._c_ptr)

    def is_P_adjustable(self):
        """
        Determines whether generator has adjustable active power.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cgen.GEN_is_P_adjustable(self._c_ptr)

    def has_flags(self, flag_type, q):
        """
        Determines whether the generator has the flags associated with
        certain quantities set.

        Parameters
        ----------
        flag_type : string (|RefFlags|)
        q : string or list of strings (|RefGeneratorQuantities|)

        Returns
        -------
        flag : |TrueFalse|
        """

        q = q if isinstance(q,list) else [q]

        return cgen.GEN_has_flags(self._c_ptr,
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

        cdef char* info_string = cgen.GEN_get_var_info_string(self._c_ptr, index)
        if info_string:
            s = info_string.decode('UTF-8')
            free(info_string)
            return s
        else:
            raise GeneratorError('index does not correspond to any variable')

    def set_P(self, P, t=0):
        """
        Sets active power.

        Parameters
        ----------
        P : float
        t : int
        """

        cgen.GEN_set_P(self._c_ptr,P,t)

    def set_Q(self, Q, t=0):
        """
        Sets reactive power.

        Parameters
        ----------
        Q : float
        t : int
        """

        cgen.GEN_set_Q(self._c_ptr,Q,t)

    property name:
        """ Generator name (string). """
        def __get__(self):
            return cgen.GEN_get_name(self._c_ptr).decode('UTF-8')
        def __set__(self,name):
            name = name.encode('UTF-8')
            cgen.GEN_set_name(self._c_ptr,name)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cgen.GEN_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (string). """
        def __get__(self): return obj2str[cgen.GEN_get_obj_type(self._c_ptr)]

    property index:
        """ Generator index (int). """
        def __get__(self): return cgen.GEN_get_index(self._c_ptr)

    property index_P:
        """ Index of generator active power variable (int or |Array|). """
        def __get__(self):
            r = [cgen.GEN_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_Q:
        """ Index of generator reactive power variable (int or |Array|). """
        def __get__(self):
            r = [cgen.GEN_get_index_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property bus:
        """ |Bus| to which generator is connected. """
        def __get__(self):
            return new_Bus(cgen.GEN_get_bus(self._c_ptr))
        def __set__(self,bus):
            cdef Bus cbus
            if not isinstance(bus,Bus) and bus is not None:
                raise GeneratorError('Not a Bus type object')
            cbus = bus
            cgen.GEN_set_bus(self._c_ptr,cbus._c_ptr if bus is not None else NULL)

    property reg_bus:
        """ |Bus| whose voltage is regulated by this generator. """
        def __get__(self):
            return new_Bus(cgen.GEN_get_reg_bus(self._c_ptr))
        def __set__(self,reg_bus):
            cdef Bus creg_bus
            if not isinstance(reg_bus,Bus) and reg_bus is not None:
                raise GeneratorError('Not a Bus type object')
            creg_bus = reg_bus
            cgen.GEN_set_reg_bus(self._c_ptr,creg_bus._c_ptr if reg_bus is not None else NULL) 

    property P:
        """ Generator active power (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cgen.GEN_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_P)
        def __set__(self,P):
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cgen.GEN_set_P(self._c_ptr,Par[t],t)

    property P_prev:
        """ Generator active power during the previous time period (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_P_prev(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_P_prev(self._c_ptr,P)

    property dP_max:
        """ Generator active power ramping limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_dP_max(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_dP_max(self._c_ptr,P)

    property P_max:
        """ Generator active power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_P_max(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_P_max(self._c_ptr,P)

    property P_min:
        """ Generator active power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_P_min(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_P_min(self._c_ptr,P)

    property Q:
        """ Generator reactive power (p.u. system base MVA) (float or |Array|). """
        def __get__(self):
            r = [cgen.GEN_get_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_Q)
        def __set__(self,Q):
            cdef int t
            cdef np.ndarray Qar = np.array(Q).flatten()
            for t in range(np.minimum(Qar.size,self.num_periods)):
                cgen.GEN_set_Q(self._c_ptr,Qar[t],t)

    property Q_max:
        """ Generator reactive power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_Q_max(self._c_ptr)
        def __set__(self,Q): cgen.GEN_set_Q_max(self._c_ptr,Q)

    property Q_min:
        """ Generator reactive power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_Q_min(self._c_ptr)
        def __set__(self,Q): cgen.GEN_set_Q_min(self._c_ptr,Q)

    property P_cost:
        """ Active power generation cost ($/hr) (float or |Array|). """
        def __get__(self):
            r = [cgen.GEN_get_P_cost(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property cost_coeff_Q0:
        """ Coefficient for genertion cost function (constant term, units of $/hr) (float). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q0(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q0(self._c_ptr,c)

    property cost_coeff_Q1:
        """ Coefficient for genertion cost function (linear term, units of $/(hr p.u.)) (float). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q1(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q1(self._c_ptr,c)

    property cost_coeff_Q2:
        """ Coefficient for genertion cost function (quadratic term, units of $/(hr p.u.^2)) (float). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q2(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q2(self._c_ptr,c)

    property outage:
        """ Flag that indicates whehter generator is on outage (boolean). """
        def __get__(self): return cgen.GEN_is_on_outage(self._c_ptr)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = cgen.GEN_get_json_string(self._c_ptr, NULL)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power upper bound (float or |Array|). """
        def __get__(self): return DoubleArray(cgen.GEN_get_sens_P_u_bound_array(self._c_ptr),
                                              cgen.GEN_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_P_u_bound[:] = x

    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power lower bound (float or |Array|). """
        def __get__(self): return DoubleArray(cgen.GEN_get_sens_P_l_bound_array(self._c_ptr),
                                              cgen.GEN_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_P_l_bound[:] = x

    property sens_Q_u_bound:
        """ Objective function sensitivity with respect to reactive power upper bound (float or |Array|). """
        def __get__(self): return DoubleArray(cgen.GEN_get_sens_Q_u_bound_array(self._c_ptr),
                                              cgen.GEN_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_Q_u_bound[:] = x

    property sens_Q_l_bound:
        """ Objective function sensitivity with respect to reactive power lower bound (float or |Array|). """
        def __get__(self): return DoubleArray(cgen.GEN_get_sens_Q_l_bound_array(self._c_ptr),
                                              cgen.GEN_get_num_periods(self._c_ptr))
        def __set__(self,x):
            self.sens_Q_l_bound[:] = x

    property flags_vars:
        """ Flags associated with variable quantities (byte). """
        def __get__(self): return cgen.GEN_get_flags_vars(self._c_ptr)

    property flags_fixed:
        """ Flags associated with fixed quantities (byte). """
        def __get__(self): return cgen.GEN_get_flags_fixed(self._c_ptr)

    property flags_bounded:
        """ Flags associated with bounded quantities (byte). """
        def __get__(self): return cgen.GEN_get_flags_bounded(self._c_ptr)

    property flags_sparse:
        """ Flags associated with sparse quantities (byte). """
        def __get__(self): return cgen.GEN_get_flags_sparse(self._c_ptr)

cdef new_Generator(cgen.Gen* g):
    if g is not NULL:
        gen = Generator(alloc=False)
        gen._c_ptr = g
        return gen
    else:
        raise GeneratorError('no generator data')
