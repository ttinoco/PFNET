#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cbat

# Infinity
BAT_INF_P = cbat.BAT_INF_P
BAT_INF_E = cbat.BAT_INF_E

class BatteryError(Exception):
    """
    Battery error exception.
    """

    pass

cdef class Battery:
    """
    Battery class.
    """

    cdef cbat.Bat* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Battery class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cbat.BAT_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

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

        cdef char* info_string = cbat.BAT_get_var_info_string(self._c_ptr, index)
        if info_string:
            s = info_string.decode('UTF-8')
            free(info_string)
            return s
        else:
            raise BatteryError('index does not correspond to any variable')

    def has_flags(self,flag_type,q):
        """
        Determines whether the battery has the flags associated with
        certain quantities set.

        Parameters
        ----------
        flag_type : string (:ref:`ref_net_flag`)
        q : string or list of strings (:ref:`ref_bat_q`)

        Returns
        -------
        flag : {``True``, ``False``}
        """

        q = q if isinstance(q,list) else [q]

        return cbat.BAT_has_flags(self._c_ptr,
                                  str2flag[flag_type],
                                  reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0))

    def set_P(self,P,t=0):
        """
        Sets battery charging power.

        Parameters
        ----------
        P : float
        t : int
        """

        cbat.BAT_set_P(self._c_ptr,P,t)

    def set_E(self,E,t=0):
        """
        Sets battery energy level.

        Parameters
        ----------
        E : float
        t : int
        """

        cbat.BAT_set_E(self._c_ptr,E,t)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cbat.BAT_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (string). """
        def __get__(self): return obj2str[cbat.BAT_get_obj_type(self._c_ptr)]

    property index:
        """ Battery index (int). """
        def __get__(self): return cbat.BAT_get_index(self._c_ptr)

    property index_Pc:
        """ Index of battery charging power variable (int or array). """
        def __get__(self):
            r = [cbat.BAT_get_index_Pc(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_Pd:
        """ Index of battery discharging power variable (int or array). """
        def __get__(self):
            r = [cbat.BAT_get_index_Pd(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_E:
        """ Index of battery energy level variable (int or array). """
        def __get__(self):
            r = [cbat.BAT_get_index_E(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which battery is connected. """
        def __get__(self): return new_Bus(cbat.BAT_get_bus(self._c_ptr))

    property P:
        """ Battery charging power (p.u. system base MVA) (float or array). """
        def __get__(self):
            r = [cbat.BAT_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_P)
        def __set__(self,P):
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cbat.BAT_set_P(self._c_ptr,Par[t],t)

    property P_max:
        """ Battery charging power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cbat.BAT_get_P_max(self._c_ptr)
        def __set__(self,P): cbat.BAT_set_P_max(self._c_ptr,P)

    property P_min:
        """ Battery charging power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cbat.BAT_get_P_min(self._c_ptr)
        def __set__(self,P): cbat.BAT_set_P_min(self._c_ptr,P)

    property E:
        """ Battery energy level at the beginning of a period (p.u. system base MVA times time unit) (float or array). """
        def __get__(self):
            r = [cbat.BAT_get_E(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_E)
        def __set__(self,E):
            cdef int t
            cdef np.ndarray Ear = np.array(E).flatten()
            for t in range(np.minimum(Ear.size,self.num_periods)):
                cbat.BAT_set_E(self._c_ptr,Ear[t],t)

    property E_init:
        """ Initial battery energy level (p.u. system base MVA times time unit) (float). """
        def __get__(self): return cbat.BAT_get_E_init(self._c_ptr)
        def __set__(self,E): cbat.BAT_set_E_init(self._c_ptr,E)

    property E_final:
        """ Battery energy level at the end of the last period (p.u. system base MVA times time unit) (float). """
        def __get__(self): return cbat.BAT_get_E_final(self._c_ptr)
        def __set__(self,E): cbat.BAT_set_E_final(self._c_ptr,E)

    property E_max:
        """ Battery energy level upper limit (p.u. system base MVA times time unit) (float). """
        def __get__(self): return cbat.BAT_get_E_max(self._c_ptr)
        def __set__(self,E): cbat.BAT_set_E_max(self._c_ptr,E)

    property eta_c:
        """ Battery charging efficiency (unitless) (float). """
        def __get__(self): return cbat.BAT_get_eta_c(self._c_ptr)
        def __set__(self,eta_c): cbat.BAT_set_eta_c(self._c_ptr,eta_c)

    property eta_d:
        """ Battery discharging efficiency (unitless) (float). """
        def __get__(self): return cbat.BAT_get_eta_d(self._c_ptr)
        def __set__(self,eta_d): cbat.BAT_set_eta_d(self._c_ptr,eta_d)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = cbat.BAT_get_json_string(self._c_ptr, NULL)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

cdef new_Battery(cbat.Bat* b):
    if b is not NULL:
        bat = Battery(alloc=False)
        bat._c_ptr = b
        return bat
    else:
        raise BatteryError('no battery data')
