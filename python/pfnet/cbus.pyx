#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cbus

# Infinite
BUS_INF_V_MAG = cbus.BUS_INF_V_MAG
BUS_INF_V_ANG = cbus.BUS_INF_V_ANG

# Sensitivities
BUS_SENS_LARGEST = cbus.BUS_SENS_LARGEST
BUS_SENS_P_BALANCE = cbus.BUS_SENS_P_BALANCE
BUS_SENS_Q_BALANCE = cbus.BUS_SENS_Q_BALANCE
BUS_SENS_V_MAG_U_BOUND = cbus.BUS_SENS_V_MAG_U_BOUND
BUS_SENS_V_MAG_L_BOUND = cbus.BUS_SENS_V_MAG_L_BOUND
BUS_SENS_V_ANG_U_BOUND = cbus.BUS_SENS_V_ANG_U_BOUND
BUS_SENS_V_ANG_L_BOUND = cbus.BUS_SENS_V_ANG_L_BOUND
BUS_SENS_V_REG_BY_GEN = cbus.BUS_SENS_V_REG_BY_GEN
BUS_SENS_V_REG_BY_TRAN = cbus.BUS_SENS_V_REG_BY_TRAN
BUS_SENS_V_REG_BY_SHUNT = cbus.BUS_SENS_V_REG_BY_SHUNT

# Mismatches
BUS_MIS_LARGEST = cbus.BUS_MIS_LARGEST
BUS_MIS_ACTIVE = cbus.BUS_MIS_ACTIVE
BUS_MIS_REACTIVE = cbus.BUS_MIS_REACTIVE

class BusError(Exception):
    """
    Bus error exception.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Bus:
    """
    Bus class.
    """

    cdef cbus.Bus* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Bus class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cbus.BUS_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_equal(self,other):
        """
        Determines whether bus is equal to given bus.

        Parameters
        ----------
        other : :class:`Bus <pfnet.Bus>`
        """

        cdef Bus b_other

        if not isinstance(other,Bus):
            return False

        b_other = other

        return cbus.BUS_is_equal(self._c_ptr,b_other._c_ptr)

    def is_slack(self):
        """
        Determines whether the bus is a slack bus.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_slack(self._c_ptr)

    def is_regulated_by_gen(self):
        """
        Determines whether the bus is regulated by a generator.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_gen(self._c_ptr)

    def is_regulated_by_tran(self):
        """
        Determines whether the bus is regulated by a transformer.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_tran(self._c_ptr)

    def is_regulated_by_shunt(self):
        """
        Determines whether the bus is regulated by a shunt device.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_shunt(self._c_ptr)

    def has_flags(self,flag_type,q):
        """
        Determines whether the bus has the flags associated with
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

        return cbus.BUS_has_flags(self._c_ptr,
                                  str2flag[flag_type],
                                  reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0))

    def get_largest_sens(self,t=0):
        """
        Gets the bus sensitivity of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        sens : float
        """

        return cbus.BUS_get_largest_sens(self._c_ptr,t)

    def get_largest_sens_type(self,t=0):
        """
        Gets the type of bus sensitivity of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        type : int
        """

        return cbus.BUS_get_largest_sens_type(self._c_ptr,t)

    def get_largest_mis(self,t=0):
        """
        Gets the bus power mismatch of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        mis : float
        """

        return cbus.BUS_get_largest_mis(self._c_ptr,t)

    def get_largest_mis_type(self,t=0):
        """
        Gets the type of bus power mismatch of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        type : int
        """

        return cbus.BUS_get_largest_mis_type(self._c_ptr,t)

    def get_quantity(self,type,t=0):
        """
        Gets the bus quantity of the given type.

        Parameters
        ----------
        type : int (:ref:`ref_bus_sens`:, :ref:`ref_bus_mis`)
        t : int (time period)

        Returns
        -------
        value : float
        """

        return cbus.BUS_get_quantity(self._c_ptr,type,t)

    def get_total_gen_P(self,t=0):
        """
        Gets the total active power injected by generators
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        P : float
        """

        return cbus.BUS_get_total_gen_P(self._c_ptr,t)

    def get_total_gen_Q(self,t=0):
        """
        Gets the total reactive power injected by generators
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        Q : float
        """

        return cbus.BUS_get_total_gen_Q(self._c_ptr,t)

    def get_total_gen_Q_max(self):
        """
        Gets the largest total reactive power that can be
        injected by generators connected to this bus.

        Returns
        -------
        Q_max : float
        """

        return cbus.BUS_get_total_gen_Q_max(self._c_ptr)

    def get_total_gen_Q_min(self):
        """
        Gets the smallest total reactive power that can be
        injected by generators connected to this bus.

        Returns
        -------
        Q_min : float
        """

        return cbus.BUS_get_total_gen_Q_min(self._c_ptr)

    def get_total_load_P(self,t=0):
        """
        Gets the total active power consumed by loads
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        P : float
        """

        return cbus.BUS_get_total_load_P(self._c_ptr,t)

    def get_total_load_Q(self,t=0):
        """
        Gets the total reactive power consumed by loads
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        Q : float
        """

        return cbus.BUS_get_total_load_Q(self._c_ptr,t)

    def get_total_shunt_g(self):
        """
        Gets the combined conductance of shunt devices
        connected to this bus.

        Returns
        -------
        g : float
        """

        return cbus.BUS_get_total_shunt_g(self._c_ptr)

    def get_total_shunt_b(self,t=0):
        """
        Gets the combined susceptance of shunt devices
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        b : float
        """

        return cbus.BUS_get_total_shunt_b(self._c_ptr,t)

    def get_num_vars(self,q,t_start=0,t_end=None):
        """
        Gets number of variables associated with the
        given quantity.

        Parameters
        ----------
        q : string or list of strings (:ref:`ref_bus_q`)
        t_start : int
        t_end : int

        Returns
        -------
        num : int
        """

        q = q if isinstance(q,list) else [q]

        if t_end is None:
            t_end = self.num_periods-1
        return cbus.BUS_get_num_vars(self._c_ptr,
                                     reduce(lambda x,y: x|y,[str2q[self.obj_type][qq] for qq in q],0),
                                     t_start,
                                     t_end)

    def set_price(self,p,t=0):
        """
        Sets bus energy price.

        Parameters
        ----------
        p : float
        t : int
        """

        cbus.BUS_set_price(self._c_ptr,p,t)

    def set_v_mag(self,v,t=0):
        """
        Sets bus voltage magnitude.

        Parameters
        ----------
        v : float
        t : int
        """

        cbus.BUS_set_v_mag(self._c_ptr,v,t)

    def set_v_ang(self,v,t=0):
        """
        Sets bus voltage angle.

        Parameters
        ----------
        v : float
        t : int
        """

        cbus.BUS_set_v_ang(self._c_ptr,v,t)

    def show(self,t=0):
        """
        Shows bus properties.

        Parameters
        ----------
        t : int (time period)
        """
        cbus.BUS_show(self._c_ptr,t)

    def __richcmp__(self,other,op):
        """
        Compares two buses.

        Parameters
        ----------
        other : Bus
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

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cbus.BUS_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (string). """
        def __get__(self): return obj2str[cbus.BUS_get_obj_type(self._c_ptr)]

    property index:
        """ Bus index (int). """
        def __get__(self): return cbus.BUS_get_index(self._c_ptr)

    property index_v_mag:
        """ Index of voltage magnitude variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_v_mag(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_v_ang:
        """ Index of voltage angle variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_v_ang(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_y:
        """ Index of voltage magnitude positive deviation variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_y(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_z:
        """ Index of voltage magnitude negative deviation variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_z(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_vl:
        """ Index of voltage low limit violation variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_vl(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_vh:
        """ Index of voltage high limit violation variable (int or array). """
        def __get__(self):
            r = [cbus.BUS_get_index_vh(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property index_P:
        """ Index of bus active power mismatch (int). """
        def __get__(self): return cbus.BUS_get_index_P(self._c_ptr)

    property index_Q:
        """ Index for bus reactive power mismatch (int). """
        def __get__(self): return cbus.BUS_get_index_Q(self._c_ptr)

    property price:
        """ Bus energy price (float or array) ($ / (hr p.u.)). """
        def __get__(self):
            r = [cbus.BUS_get_price(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_price)
        def __set__(self,p):
            cdef int t
            cdef np.ndarray par = np.array(p).flatten()
            for t in range(np.minimum(par.size,self.num_periods)):
                cbus.BUS_set_price(self._c_ptr,par[t],t)

    property number:
        """ Bus number (int). """
        def __get__(self): return cbus.BUS_get_number(self._c_ptr)

    property name:
        """ Bus name (sting). """
        def __get__(self): return cbus.BUS_get_name(self._c_ptr).decode('UTF-8')
        def __set__(self,name):
            name = name.encode('UTF-8')
            cbus.BUS_set_name(self._c_ptr,name)

    property degree:
        """ Bus degree (number of incident branches) (float). """
        def __get__(self): return cbus.BUS_get_degree(self._c_ptr)

    property v_mag:
        """ Bus volatge magnitude (p.u. bus base kv) (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_v_mag(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_v_mag)
        def __set__(self,v):
            cdef int t
            cdef np.ndarray var = np.array(v).flatten()
            for t in range(np.minimum(var.size,self.num_periods)):
                cbus.BUS_set_v_mag(self._c_ptr,var[t],t)

    property v_ang:
        """ Bus voltage angle (radians) (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_v_ang(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return AttributeArray(r,self.set_v_ang)
        def __set__(self,v):
            cdef int t
            cdef np.ndarray var = np.array(v).flatten()
            for t in range(np.minimum(var.size,self.num_periods)):
                cbus.BUS_set_v_ang(self._c_ptr,var[t],t)

    property v_set:
        """ Bus voltage set point (p.u. bus base kv) (float or array). Equals one if bus is not regulated by a generator. """
        def __get__(self):
            r = [cbus.BUS_get_v_set(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property v_max:
        """ Bus volatge upper bound (p.u. bus base kv) (float). """
        def __get__(self): return cbus.BUS_get_v_max(self._c_ptr)

    property v_min:
        """ Bus voltage lower bound (p.u. bus base kv) (float). """
        def __get__(self): return cbus.BUS_get_v_min(self._c_ptr)

    property P_mis:
        """ Bus active power mismatch (p.u. system base MVA) (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_P_mis(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property Q_mis:
        """ Bus reactive power mismatch (p.u. system base MVA) (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_Q_mis(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_P_balance:
        """ Objective function sensitivity with respect to bus active power balance (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_P_balance(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_Q_balance:
        """ Objective function sensitivity with respect to bus reactive power balance (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_Q_balance(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_mag_u_bound:
        """ Objective function sensitivity with respect to voltage magnitude upper bound (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_mag_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_mag_l_bound:
        """ Objective function sensitivity with respect to voltage magnitude lower bound (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_mag_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_ang_u_bound:
        """ Objective function sensitivity with respect to voltage angle upper bound (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_ang_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_ang_l_bound:
        """ Objective function sensitivity with respect to voltage angle lower bound (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_ang_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_reg_by_gen:
        """ Objective function sensitivity with respect to bus voltage regulation by generators (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_reg_by_gen(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_reg_by_tran:
        """ Objective function sensitivity with respect to bus voltage regulation by transformers (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_reg_by_tran(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property sens_v_reg_by_shunt:
        """ Objective function sensitivity with respect to bus voltage regulation by shunts (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_v_reg_by_shunt(self._c_ptr,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property generators:
        """ List of :class:`generators <pfnet.Generator>` connected to this bus (list). """
        def __get__(self):
            gens = []
            cdef cgen.Gen* g = cbus.BUS_get_gen(self._c_ptr)
            while g is not NULL:
                gens.append(new_Generator(g))
                g = cgen.GEN_get_next(g)
            return gens

    property gens:
        """ Same as :attr:`generators <pfnet.Bus.generators>`. """
        def __get__(self): return self.generators

    property reg_generators:
        """ List of :class:`generators <pfnet.Generator>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_gens = []
            cdef cgen.Gen* g = cbus.BUS_get_reg_gen(self._c_ptr)
            while g is not NULL:
                reg_gens.append(new_Generator(g))
                g = cgen.GEN_get_reg_next(g)
            return reg_gens

    property reg_gens:
        """ Same as :attr:`reg_generators <pfnet.Bus.reg_generators>`. """
        def __get__(self): return self.reg_generators

    property reg_trans:
        """ List of :class:`tap-changing transformers <pfnet.Branch>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_trans = []
            cdef cbranch.Branch* br = cbus.BUS_get_reg_tran(self._c_ptr)
            while br is not NULL:
                reg_trans.append(new_Branch(br))
                br = cbranch.BRANCH_get_reg_next(br)
            return reg_trans

    property reg_shunts:
        """ List of :class:`switched shunt devices <pfnet.Shunt>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_shunts = []
            cdef cshunt.Shunt* s = cbus.BUS_get_reg_shunt(self._c_ptr)
            while s is not NULL:
                reg_shunts.append(new_Shunt(s))
                s = cshunt.SHUNT_get_reg_next(s)
            return reg_shunts

    property branches_from:
        """ .. deprecated:: 1.2.5  Same as :attr:`branches_k <pfnet.Bus.branches_k>`. """
        def __get__(self): return self.branches_k

    property branches_k:
        """ List of :class:`branches <pfnet.Branch>` that have this bus on the "k" (aka "from" or "i") side (list). """
        def __get__(self):
            branches = []
            cdef cbranch.Branch* br = cbus.BUS_get_branch_k(self._c_ptr)
            while br is not NULL:
                branches.append(new_Branch(br))
                br = cbranch.BRANCH_get_next_k(br)
            return branches

    property branches_to:
        """ .. deprecated:: 1.2.5  Same as :attr:`branches_m <pfnet.Bus.branches_m>`. """
        def __get__(self): return self.branches_m

    property branches_m:
        """ List of :class:`branches <pfnet.Branch>` that have this bus on the "m" (aka "to" or "j") side (list). """
        def __get__(self):
            branches = []
            cdef cbranch.Branch* br = cbus.BUS_get_branch_m(self._c_ptr)
            while br is not NULL:
                branches.append(new_Branch(br))
                br = cbranch.BRANCH_get_next_m(br)
            return branches

    property branches:
        """ List of :class:`branches <pfnet.Branch>` incident on this bus (list). """
        def __get__(self):
            # combine both "k"/"from" and "m"/"to" branches
            return self.branches_k+self.branches_m

    property loads:
        """ List of :class:`loads <pfnet.Load>` connected to this bus (list). """
        def __get__(self):
            loads = []
            cdef cload.Load* l = cbus.BUS_get_load(self._c_ptr)
            while l is not NULL:
                loads.append(new_Load(l))
                l = cload.LOAD_get_next(l)
            return loads

    property var_generators:
        """ List of :class:`variable generators <pfnet.VarGenerator>` connected to this bus (list). """
        def __get__(self):
            vargens = []
            cdef cvargen.Vargen* g = cbus.BUS_get_vargen(self._c_ptr)
            while g is not NULL:
                vargens.append(new_VarGenerator(g))
                g = cvargen.VARGEN_get_next(g)
            return vargens

    property var_gens:
        """ Same as :attr:`var_generators <pfnet.Bus.var_generators>`. """
        def __get__(self): return self.var_generators

    property batteries:
        """ List of :class:`batteries <pfnet.Battery>` connected to this bus (list). """
        def __get__(self):
            bats = []
            cdef cbat.Bat* b = cbus.BUS_get_bat(self._c_ptr)
            while b is not NULL:
                bats.append(new_Battery(b))
                b = cbat.BAT_get_next(b)
            return bats

    property bats:
        """ Same as :attr:`batteries <pfnet.Bus.batteries>`. """
        def __get__(self): return self.batteries


cdef new_Bus(cbus.Bus* b):
    if b is not NULL:
        bus = Bus(alloc=False)
        bus._c_ptr = b
        return bus
    else:
        raise BusError('no bus data')
