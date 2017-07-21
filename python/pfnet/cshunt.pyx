#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cshunt

# Infinite
SHUNT_INF_SUSC = cshunt.SHUNT_INF_SUSC

class ShuntError(Exception):
    """
    Shunt error exception.
    """

    pass

cdef class Shunt:
    """
    Shunt class.
    """

    cdef cshunt.Shunt* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Shunt class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cshunt.SHUNT_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_fixed(self):
        """
        Determines whether the shunt device is fixed (as opposed to switched).

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cshunt.SHUNT_is_fixed(self._c_ptr)

    def is_switched_v(self):
        """
        Determines whether the shunt is switchable and regulates
        bus voltage magnitude.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cshunt.SHUNT_is_switched_v(self._c_ptr)

    def has_flags(self,flag_type,q):
        """
        Determines whether the shunt devices has flags associated with
        certain quantities set.

        Parameters
        ----------
        flag_type : string (:ref:`ref_net_flag`)
        q : string or list of strings (:ref:`ref_bus_q`)

        Returns
        -------
        flag : {``True``, ``False``}
        """

        q = q if isinstance(q,list) else [q]

        return cshunt.SHUNT_has_flags(self._c_ptr,
                                      str2flag[flag_type],
                                      reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0))

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cshunt.SHUNT_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (string). """
        def __get__(self): return obj2str[cshunt.SHUNT_get_obj_type(self._c_ptr)]

    property index:
        """ Shunt index (int). """
        def __get__(self): return cshunt.SHUNT_get_index(self._c_ptr)

    property index_b:
        """ Index of shunt susceptance variable (int or array). """
        def __get__(self):
            r = [cshunt.SHUNT_get_index_b(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which the shunt devices is connected. """
        def __get__(self): return new_Bus(cshunt.SHUNT_get_bus(self._c_ptr))

    property reg_bus:
        """ :class:`Bus <pfnet.Bus>` whose voltage magnitude is regulated by this shunt device. """
        def __get__(self): return new_Bus(cshunt.SHUNT_get_reg_bus(self._c_ptr))

    property g:
        """ Shunt conductance (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_g(self._c_ptr)

    property b:
        """ Shunt susceptance (p.u.) (float or array). """
        def __get__(self):
            r = [cshunt.SHUNT_get_b(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property b_max:
        """ Shunt susceptance upper limit (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_b_max(self._c_ptr)
        def __set__(self,value): cshunt.SHUNT_set_b_max(self._c_ptr,value)

    property b_min:
        """ Shunt susceptance lower limit (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_b_min(self._c_ptr)
        def __set__(self,value): cshunt.SHUNT_set_b_min(self._c_ptr,value)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = cshunt.SHUNT_get_json_string(self._c_ptr, NULL)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

cdef new_Shunt(cshunt.Shunt* s):
    if s is not NULL:
        shunt = Shunt(alloc=False)
        shunt._c_ptr = s
        return shunt
    else:
        raise ShuntError('no shunt data')
