#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cnet

class NetworkError(Exception):
    """
    Network error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Network:
    """
    Network class.
    """

    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self,num_periods=1,alloc=True):
        """
        Network class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_net = cnet.NET_new(num_periods)
        else:
            self._c_net = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees network C data structure.
        """

        if self.alloc:
            cnet.NET_del(self._c_net)
            self._c_net = NULL

    def add_vargens(self,buses,penetration,uncertainty,corr_radius,corr_value):
        """
        Adds variable generators to the network.

        Parameters
        ----------

        buses : list of :class:`Buses <pfnet.Bus>`
        penetration : float
                      percentage
        uncertainty : float
                      percentage
        corr_radius : int
                      number of branches
        corr_value : float
                     correlation coefficient
        """

        cdef Bus head = buses[0] if buses else None
        cdef Bus prev = head
        cdef Bus curr
        for b in buses[1:]:
            curr = b
            cbus.BUS_set_next(prev._c_ptr,curr._c_ptr)
            prev = curr
        if prev is not None:
            cbus.BUS_set_next(prev._c_ptr,NULL)

        if head:
            cnet.NET_add_vargens(self._c_net,head._c_ptr,penetration,uncertainty,corr_radius,corr_value)
        else:
            cnet.NET_add_vargens(self._c_net,NULL,penetration,uncertainty,corr_radius,corr_value)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def adjust_generators(self):
        """
        Adjusts powers of slack and regulator generators connected to or regulating the
        same bus to correct generator participations without modifying the total power injected.
        """

        cnet.NET_adjust_generators(self._c_net);

    def clear_error(self):
        """
        Clear error flag and message string.
        """

        cnet.NET_clear_error(self._c_net);

    def clear_flags(self):
        """
        Clears all the flags of all the network components.
        """

        cnet.NET_clear_flags(self._c_net)

    def clear_properties(self):
        """
        Clears all the network properties.
        """

        cnet.NET_clear_properties(self._c_net)

    def clear_sensitivities(self):
        """
        Clears all sensitivity information.
        """

        cnet.NET_clear_sensitivities(self._c_net)

    def create_sorted_bus_list(self,sort_by,t=0):
        """
        Creates list of buses sorted in descending order according to a specific quantity.

        Parameters
        ----------
        sort_by : int (:ref:`ref_bus_sens`, :ref:`ref_bus_mis`).
        t : int

        Returns
        -------
        buses : list of :class:`Buses <pfnet.Bus>`
        """

        buses = []
        cdef cbus.Bus* b = cnet.NET_create_sorted_bus_list(self._c_net,sort_by,t)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def create_vargen_P_sigma(self,spread,corr):
        """
        Creates covariance matrix (lower triangular part) for
        variable vargen active powers.

        Parameters
        ----------
        spead : int
                Determines correlation neighborhood in terms of number of edges.
        corr : float
               Desired correlation coefficient for neighboring vargens.

        Returns
        -------
        sigma : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """

        sigma = Matrix(cnet.NET_create_vargen_P_sigma(self._c_net,spread,corr),
                       owndata=True)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))
        else:
            return sigma

    def get_bus_by_number(self,number):
        """
        Gets bus with the given number.

        Parameters
        ----------
        number : int

        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        ptr = cnet.NET_bus_hash_number_find(self._c_net,number)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_bus_by_name(self,name):
        """
        Gets bus with the given name.

        Parameters
        ----------
        name : string

        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_bus_hash_name_find(self._c_net,name)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_vargen_by_name(self,name):
        """
        Gets vargen with the given name.

        Parameters
        ----------
        name : string

        Returns
        -------
        vargen : :class:`VarGenerator <pfnet.VarGenerator>`
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_vargen_hash_name_find(self._c_net,name)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('vargen not found')

    def get_bus(self,index):
        """
        Gets bus with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        ptr = cnet.NET_get_bus(self._c_net,index)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('invalid bus index')

    def get_branch(self,index):
        """
        Gets branch with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        branch : :class:`Branch <pfnet.Branch>`
        """

        ptr = cnet.NET_get_branch(self._c_net,index)
        if ptr is not NULL:
            return new_Branch(ptr)
        else:
            raise NetworkError('invalid branch index')

    def get_gen(self,index):
        """
        Gets generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Generator <pfnet.Generator>`
        """

        ptr = cnet.NET_get_gen(self._c_net,index)
        if ptr is not NULL:
            return new_Generator(ptr)
        else:
            raise NetworkError('invalid gen index')

    def get_shunt(self,index):
        """
        Gets shunt with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Shunt <pfnet.Shunt>`
        """

        ptr = cnet.NET_get_shunt(self._c_net,index)
        if ptr is not NULL:
            return new_Shunt(ptr)
        else:
            raise NetworkError('invalid shunt index')

    def get_load(self,index):
        """
        Gets load with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Load <pfnet.Load>`
        """

        ptr = cnet.NET_get_load(self._c_net,index)
        if ptr is not NULL:
            return new_Load(ptr)
        else:
            raise NetworkError('invalid load index')

    def get_vargen(self,index):
        """
        Gets variable generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        vargen : :class:`VarGenerator <pfnet.VarGenerator>`
        """

        ptr = cnet.NET_get_vargen(self._c_net,index)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('invalid vargen index')

    def get_bat(self,index):
        """
        Gets battery with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bat : :class:`Battery <pfnet.Battery>`
        """

        ptr = cnet.NET_get_bat(self._c_net,index)
        if ptr is not NULL:
            return new_Battery(ptr)
        else:
            raise NetworkError('invalid battery index')

    def get_gen_buses(self):
        """
        Gets list of buses where generators are connected.

        Returns
        -------
        buses : list
        """

        buses = []
        cdef cbus.Bus* b = cnet.NET_get_gen_buses(self._c_net)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def get_load_buses(self):
        """
        Gets list of buses where loads are connected.

        Returns
        -------
        buses : list
        """

        buses = []
        cdef cbus.Bus* b = cnet.NET_get_load_buses(self._c_net)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def get_var_values(self,option='current'):
        """
        Gets network variable values.

        Parameters
        ----------
        option : string (See var values)

        Returns
        -------
        values : :class:`ndarray <numpy.ndarray>`
        """
        return Vector(cnet.NET_get_var_values(self._c_net,str2const[option]),owndata=True)

    def get_var_projection(self,obj_type,q,t_start=0,t_end=None):
        """
        Gets projection matrix for specific object variables.

        Parameters
        ----------
        obj_type : string (:ref:`ref_net_obj`)
        q : string or list of strings (:ref:`ref_bus_q`, :ref:`ref_branch_q`, :ref:`ref_gen_q`, :ref:`ref_shunt_q`, :ref:`ref_load_q`, :ref:`ref_vargen_q`, :ref:`ref_bat_q`)
        t_start : int
        t_end : int (inclusive)
        """

        q = q if isinstance(q,list) else [q]

        if t_end is None:
            t_end = self.num_periods-1
        m = Matrix(cnet.NET_get_var_projection(self._c_net,
                                               str2obj[obj_type],
                                               reduce(lambda x,y: x|y,[str2q[obj_type][qq] for qq in q],0),
                                               t_start,
                                               t_end),
                   owndata=True)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))
        else:
            return m

    def get_num_buses(self):
        """
        Gets number of buses in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_buses(self._c_net)

    def get_num_slack_buses(self):
        """
        Gets number of slack buses in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_slack_buses(self._c_net)

    def get_num_buses_reg_by_gen(self):
        """
        Gets number of buses whose voltage magnitudes are regulated by generators.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_buses_reg_by_gen(self._c_net)

    def get_num_buses_reg_by_tran(self,only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by tap-changing transformers.

        Returns
        -------
        num : int
        """

        if not only:
            return cnet.NET_get_num_buses_reg_by_tran(self._c_net)
        else:
            return cnet.NET_get_num_buses_reg_by_tran_only(self._c_net)

    def get_num_buses_reg_by_shunt(self,only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by switched shunt devices.

        Returns
        -------
        num : int
        """

        if not only:
            return cnet.NET_get_num_buses_reg_by_shunt(self._c_net)
        else:
            return cnet.NET_get_num_buses_reg_by_shunt_only(self._c_net)

    def get_num_branches(self):
        """
        Gets number of branches in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_branches(self._c_net)

    def get_num_branches_not_on_outage(self):
        """
        Gets number of branches in the network that are not on outage.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_branches_not_on_outage(self._c_net)

    def get_num_fixed_trans(self):
        """
        Gets number of fixed transformers in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_fixed_trans(self._c_net)

    def get_num_lines(self):
        """
        Gets number of transmission lines in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_lines(self._c_net)

    def get_num_phase_shifters(self):
        """
        Gets number of phase-shifting transformers in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_phase_shifters(self._c_net)

    def get_num_tap_changers(self):
        """
        Gets number of tap-changing transformers in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_tap_changers(self._c_net)

    def get_num_tap_changers_v(self):
        """
        Gets number of tap-changing transformers in the network that regulate voltage magnitudes.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_tap_changers_v(self._c_net)

    def get_num_tap_changers_Q(self):
        """
        Gets number of tap-changing transformers in the network that regulate reactive flows.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_tap_changers_Q(self._c_net)

    def get_num_generators(self):
        """
        Gets number of generators in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_gens(self._c_net)

    def get_num_gens(self):
        """ Same as :attr:`get_num_generators <pfnet.Network.get_num_generators>`. """
        return self.get_num_generators()

    def get_num_gens_not_on_outage(self):
        """
        Gets number of generators in the network that are not on outage.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_gens_not_on_outage(self._c_net)

    def get_num_reg_gens(self):
        """
        Gets number generators in the network that provide voltage regulation.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_reg_gens(self._c_net)

    def get_num_slack_gens(self):
        """
        Gets number of slack generators in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_slack_gens(self._c_net)

    def get_num_P_adjust_gens(self):
        """
        Gets number of generators in the network that have adjustable active powers.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_P_adjust_gens(self._c_net)

    def get_num_loads(self):
        """
        Gets number of loads in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_loads(self._c_net)

    def get_num_P_adjust_loads(self):
        """
        Gets number of loads in the network that have adjustable active powers.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_P_adjust_loads(self._c_net)

    def get_num_shunts(self):
        """
        Gets number of shunts in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_shunts(self._c_net)

    def get_num_fixed_shunts(self):
        """
        Gets number of fixed shunts in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_fixed_shunts(self._c_net)

    def get_num_switched_shunts(self):
        """
        Gets number of switched shunts in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_switched_shunts(self._c_net)

    def get_num_var_generators(self):
        """
        Gets number of variable generators in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_vargens(self._c_net)

    def get_num_var_gens(self):
        """ Same as :attr:`get_num_var_generators <pfnet.Network.get_num_var_generators>`. """
        return self.get_num_var_generators()

    def get_num_batteries(self):
        """
        Gets number of batteries in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_bats(self._c_net)

    def get_num_bats(self):
        """ Same as :attr:`get_num_batteries <pfnet.Network.get_num_batteries>`. """
        return self.get_num_bats()

    def get_properties(self):
        """
        Gets network properties.

        Returns
        -------
        properties : dict
        """

        return {'bus_v_max': self.bus_v_max,
                'bus_v_min': self.bus_v_min,
                'bus_v_vio': self.bus_v_vio,
                'bus_P_mis': self.bus_P_mis,
                'bus_Q_mis': self.bus_Q_mis,
                'gen_P_cost': self.gen_P_cost,
                'gen_v_dev': self.gen_v_dev,
                'gen_Q_vio': self.gen_Q_vio,
                'gen_P_vio': self.gen_P_vio,
                'tran_v_vio': self.tran_v_vio,
                'tran_r_vio': self.tran_r_vio,
                'tran_p_vio': self.tran_p_vio,
                'shunt_v_vio': self.shunt_v_vio,
                'shunt_b_vio': self.shunt_b_vio,
                'load_P_util': self.load_P_util,
                'load_P_vio': self.load_P_vio,
                'num_actions': self.num_actions}

    def has_error(self):
        """
        Indicates whether the network has the error flag set due to an
        invalid operation.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cnet.NET_has_error(self._c_net)

    def load(self,filename,output_level=0):
        """
        Loads a network data contained in a specific file.

        Parameters
        ----------
        filename : string
        output_level : int
        """

        filename = filename.encode('UTF-8')
        cnet.NET_load(self._c_net,filename,output_level)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def set_flags(self,obj_type,flags,props,q):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj_type : string (:ref:`ref_net_obj`)
        flags : string or list of strings (:ref:`ref_net_flag`)
        props : string or list of strings (:ref:`ref_bus_prop`, :ref:`ref_branch_prop`, :ref:`ref_gen_prop`, :ref:`ref_shunt_prop`, :ref:`ref_load_prop`, :ref:`ref_vargen_prop`, :ref:`ref_bat_prop`)
        q : string or list of strings (:ref:`ref_bus_q`, :ref:`ref_branch_q`, :ref:`ref_gen_q`, :ref:`ref_shunt_q`, :ref:`ref_load_q`, :ref:`ref_vargen_q`, :ref:`ref_bat_q`)
        """

        flags = flags if isinstance(flags,list) else [flags]
        props = props if isinstance(props,list) else [props]
        q = q if isinstance(q,list) else [q]
        cnet.NET_set_flags(self._c_net,
                           str2obj[obj_type],
                           reduce(lambda x,y: x|y,[str2flag[f] for f in flags],0),
                           reduce(lambda x,y: x|y,[str2prop[obj_type][pp] for pp in props],0),
                           reduce(lambda x,y: x|y,[str2q[obj_type][qq] for qq in q],0))
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def set_flags_of_component(self,obj,flags,q):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj : :class:`Bus <pfnet.Bus>`, :class:`Branch <pfnet.Branch>`, :class:`Generator <pfnet.Generator>`, :class:`Load <pfnet.Load>`, :class:`Shunt <pfnet.Shunt>`, :class:`VarGenerator <pfnet.VarGenerator>`, :class:`Battery <pfnet.Battery>`
        flags : string or list of strings (:ref:`ref_net_flag`)
        q : string or list of strings (:ref:`ref_bus_q`, :ref:`ref_branch_q`, :ref:`ref_gen_q`, :ref:`ref_shunt_q`, :ref:`ref_load_q`, :ref:`ref_vargen_q`, :ref:`ref_bat_q`)
        """

        cdef CPtr ptr = obj._get_c_ptr()
        flags = flags if isinstance(flags,list) else [flags]
        q = q if isinstance(q,list) else [q]
        cnet.NET_set_flags_of_component(self._c_net,
                                        ptr._c_ptr,
                                        str2obj[obj.obj_type],
                                        reduce(lambda x,y: x|y,[str2flag[f] for f in flags],0),
                                        reduce(lambda x,y: x|y,[str2q[obj.obj_type][qq] for qq in q],0))
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def set_var_values(self,values):
        """
        Sets network variable values.

        Parameters
        ----------
        values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if values.size else NULL
        cnet.NET_set_var_values(self._c_net,v)

    def show_components(self):
        """
        Shows information about the number of network components of each type.
        """

        print(cnet.NET_get_show_components_str(self._c_net).decode('UTF-8'))

    def show_properties(self,t=0):
        """
        Shows information about the state of the network component quantities.

        Parameters
        ----------
        t : int (time period)
        """

        print(cnet.NET_get_show_properties_str(self._c_net,t).decode('UTF-8'))

    def show_buses(self,number,sort_by,t=0):
        """
        Shows information about the most relevant network buses sorted by a specific quantity.

        Parameters
        ----------
        number : int
        sort_by : int (:ref:`ref_bus_sens`, :ref:`ref_bus_mis`)
        t : int (time period)
        """

        cnet.NET_show_buses(self._c_net,number,sort_by,t)

    def update_properties(self,values=None):
        """
        Re-computes the network properties using the given values
        of the network variables. If no values are given, then the
        current values of the network variables are used.

        Parameters
        ----------
        values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if (values is not None and values.size) else NULL
        cnet.NET_update_properties(self._c_net,v)

    def update_set_points(self):
        """
        Updates voltage magnitude set points of gen-regulated buses
        to be equal to the bus voltage magnitudes.
        """

        cnet.NET_update_set_points(self._c_net)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cnet.NET_get_num_periods(self._c_net)

    property base_power:
        """ System base power (MVA) (float). """
        def __get__(self): return cnet.NET_get_base_power(self._c_net)

    property buses:
        """ List of network :class:`buses <pfnet.Bus>` (list). """
        def __get__(self):
            return [self.get_bus(i) for i in range(self.num_buses)]

    property branches:
        """ List of network :class:`branches <pfnet.Branch>` (list). """
        def __get__(self):
            return [self.get_branch(i) for i in range(self.num_branches)]

    property generators:
        """ List of network :class:`generators <pfnet.Generator>` (list). """
        def __get__(self):
            return [self.get_gen(i) for i in range(self.num_generators)]

    property gens:
        """ Same as :attr:`generators <pfnet.Network.generators>` """
        def __get__(self): return self.generators

    property shunts:
        """ List of network :class:`shunts <pfnet.Shunt>` (list). """
        def __get__(self):
            return [self.get_shunt(i) for i in range(self.num_shunts)]

    property loads:
        """ List of network :class:`loads <pfnet.Load>` (list). """
        def __get__(self):
            return [self.get_load(i) for i in range(self.num_loads)]

    property var_generators:
        """ List of network :class:`variable generators <pfnet.VarGenerator>` (list). """
        def __get__(self):
            return [self.get_vargen(i) for i in range(self.num_var_generators)]

    property var_gens:
        """ Same as :attr:`var_generators <pfnet.Network.var_generators>`. """
        def __get__(self): return self.var_generators

    property batteries:
        """ List of network :class:`batteries <pfnet.Battery>` (list). """
        def __get__(self):
            return [self.get_bat(i) for i in range(self.num_batteries)]

    property bats:
        """ Same as :attr:`batteries <pfnet.Network.batteries>`. """
        def __get__(self): return self.batteries

    property num_buses:
        """ Number of buses in the network (int). """
        def __get__(self): return cnet.NET_get_num_buses(self._c_net)

    property num_branches:
        """ Number of branches in the network (int). """
        def __get__(self): return cnet.NET_get_num_branches(self._c_net)

    property num_generators:
        """ Number of generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_gens(self._c_net)

    property num_gens:
        """ Same as :attr:`num_generators <pfnet.Network.num_generators>`. """
        def __get__(self): return self.num_generators

    property num_loads:
        """ Number of loads in the network (int). """
        def __get__(self): return cnet.NET_get_num_loads(self._c_net)

    property num_shunts:
        """ Number of shunt devices in the network (int). """
        def __get__(self): return cnet.NET_get_num_shunts(self._c_net)

    property num_var_generators:
        """ Number of variable generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_vargens(self._c_net)

    property num_vargens:
        """ Same as :attr:`num_var_generators <pfnet.Network.num_var_generators>`. """
        def __get__(self): return self.num_var_generators

    property num_batteries:
        """ Number of batteries in the network (int). """
        def __get__(self): return cnet.NET_get_num_bats(self._c_net)

    property num_bats:
        """ Same as :attr:`num_batteries <pfnet.Network.num_batteries>`. """
        def __get__(self): return self.num_batteries

    property num_vars:
        """ Number of network quantities that have been set to variable (int). """
        def __get__(self): return cnet.NET_get_num_vars(self._c_net)

    property num_fixed:
        """ Number of network quantities that have been set to fixed (int). """
        def __get__(self): return cnet.NET_get_num_fixed(self._c_net)

    property num_bounded:
        """ Number of network quantities that have been set to bounded (int). """
        def __get__(self): return cnet.NET_get_num_bounded(self._c_net)

    property num_sparse:
        """ Number of network control quantities that have been set to sparse (int). """
        def __get__(self): return cnet.NET_get_num_sparse(self._c_net)

    property bus_v_max:
        """ Maximum bus voltage magnitude (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_max(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_v_min:
        """ Minimum bus voltage magnitude (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_min(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_v_vio:
        """ Maximum bus voltage magnitude limit violation (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_P_mis:
        """ Largest bus active power mismatch in the network (MW) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_bus_P_mis(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_Q_mis:
        """ Largest bus reactive power mismatch in the network (MVAr) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_bus_Q_mis(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_P_cost:
        """ Total active power generation cost ($/hr) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_gen_P_cost(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_v_dev:
        """ Largest voltage magnitude deviation from set point of bus regulated by generator (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_gen_v_dev(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_Q_vio:
        """ Largest generator reactive power limit violation (MVAr) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_gen_Q_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_P_vio:
        """ Largest generator active power limit violation (MW) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_gen_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by transformer (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_tran_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_r_vio:
        """ Largest transformer tap ratio limit violation (float or array). """
        def __get__(self):
            r = [cnet.NET_get_tran_r_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_p_vio:
        """ Largest transformer phase shift limit violation (float or array). """
        def __get__(self):
            r = [cnet.NET_get_tran_p_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property shunt_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by switched shunt device (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_shunt_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property shunt_b_vio:
        """ Largest switched shunt susceptance limit violation (p.u.) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_shunt_b_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property load_P_util:
        """ Total active power consumption utility ($/hr) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_load_P_util(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property load_P_vio:
        """ Largest load active power limit violation (MW) (float or array). """
        def __get__(self):
            r = [cnet.NET_get_load_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property num_actions:
        """ Number of control adjustments (int or array). """
        def __get__(self):
            r = [cnet.NET_get_num_actions(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property vargen_corr_radius:
        """ Correlation radius of variable generators (number of edges). """
        def __get__(self): return cnet.NET_get_vargen_corr_radius(self._c_net)

    property vargen_corr_value:
        """ Correlation value (coefficient) of variable generators. """
        def __get__(self): return cnet.NET_get_vargen_corr_value(self._c_net)

cdef new_Network(cnet.Net* n):
    if n is not NULL:
        net = Network(alloc=False)
        net._c_net = n
        return net
    else:
        raise NetworkError('no network data')
