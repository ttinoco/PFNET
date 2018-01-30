#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import os
cimport cnet
import tempfile

class NetworkError(Exception):
    """
    Network error exception.
    """

    pass

cdef class Network:
    """
    Network class.
    """

    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self, num_periods=1, alloc=True):
        """
        Network class.

        Parameters
        ----------
        num_periods : int
        alloc : |TrueFalse|
        """

        pass

    def __cinit__(self, num_periods=1, alloc=True):

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

    def __getstate__(self):
        
        return  self.json_string

    def __setstate__(self, state):
        
        cdef Network new_net
        if self._c_net != NULL:
            cnet.NET_del(self._c_net)
            self._c_net = NULL
            
        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            f.write(state.encode('UTF-8'))
            f.seek(0)
            f.close()
            new_net = ParserJSON().parse(f.name)
            self._c_net = new_net._c_net
            new_net.alloc = False
            os.remove(f.name)

    def extract_subnetwork(self, buses):
        """
        Extracts subnetwork containig the given buses.

        Parameters
        ----------
        buses : list of |Bus| objects

        Returns
        -------
        net : |Network|
        """

        cdef Network net
        cdef cnet.Bus** array
        cdef Bus bus

        array = <cnet.Bus**>malloc(len(buses)*sizeof(cnet.Bus*))
        for i in range(len(buses)):
            bus = buses[i]
            array[i] = bus._c_ptr

        net = new_Network(cnet.NET_extract_subnet(self._c_net, array, len(buses)))
        net.alloc = True
        free(array)

        return net       

    def add_buses(self, buses):
        """
        Adds buses to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        buses : list of |Bus| objects
        """

        cdef cnet.Bus** array
        cdef Bus bus

        old_num_buses = self.num_buses
        array = <cnet.Bus**>malloc(len(buses)*sizeof(cnet.Bus*))
        for i in range(len(buses)):
            bus = buses[i]
            array[i] = bus._c_ptr

        cnet.NET_add_buses(self._c_net, array, len(buses))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(buses)):
            bus = buses[i]
            if array[i] != NULL:
                if bus.alloc:
                    cbus.BUS_array_del(bus._c_ptr,1)
                bus._c_ptr = cnet.NET_get_bus(self._c_net, old_num_buses+index)
                bus.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_buses(self, buses):
        """
        Removes buses from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        buses : list of |Bus| objects
        """

        cdef cnet.Bus** array
        cdef Bus bus

        array = <cnet.Bus**>malloc(len(buses)*sizeof(cnet.Bus*))
        for i in range(len(buses)):
            bus = buses[i]
            array[i] = bus._c_ptr

        cnet.NET_del_buses(self._c_net, array, len(buses))

        # Update pointers and alloc flags
        for i in range(len(buses)):
            bus = buses[i]
            if array[i] != NULL:
                bus._c_ptr = NULL
                bus.alloc = False

        # Clean up
        free(array)
            
    def add_branches(self, branches):
        """
        Adds branches to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        branches : list of |Branch| objects
        """

        cdef cnet.Branch** array
        cdef Branch br

        old_num_branches = self.num_branches
        array = <cnet.Branch**>malloc(len(branches)*sizeof(cnet.Branch*))
        for i in range(len(branches)):
            br = branches[i]
            array[i] = br._c_ptr

        cnet.NET_add_branches(self._c_net, array, len(branches))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(branches)):
            br = branches[i]
            if array[i] != NULL:
                if br.alloc:
                    cbranch.BRANCH_array_del(br._c_ptr,1)
                br._c_ptr = cnet.NET_get_branch(self._c_net, old_num_branches+index)
                br.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_branches(self, branches):
        """
        Removes branches from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        branches : list of |Branch| objects
        """

        cdef cnet.Branch** array
        cdef Branch br

        array = <cnet.Branch**>malloc(len(branches)*sizeof(cnet.Branch*))
        for i in range(len(branches)):
            br = branches[i]
            array[i] = br._c_ptr

        cnet.NET_del_branches(self._c_net, array, len(branches))

        # Update pointers and alloc flags
        for i in range(len(branches)):
            br = branches[i]
            if array[i] != NULL:
                br._c_ptr = NULL
                br.alloc = False

        # Clean up
        free(array)

    def add_generators(self, generators):
        """
        Adds generators to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        generators : list of |Generator| objects
        """

        cdef cnet.Gen** array
        cdef Generator gen

        old_num_gens = self.num_generators
        array = <cnet.Gen**>malloc(len(generators)*sizeof(cnet.Gen*))
        for i in range(len(generators)):
            gen = generators[i]
            array[i] = gen._c_ptr

        cnet.NET_add_gens(self._c_net, array, len(generators))
        
        # Update pointers and alloc flags
        index = 0
        for i in range(len(generators)):
            gen = generators[i]
            if array[i] != NULL:
                if gen.alloc:
                    cgen.GEN_array_del(gen._c_ptr,1)
                gen._c_ptr = cnet.NET_get_gen(self._c_net, old_num_gens+index)
                gen.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_generators(self, generators):
        """
        Removes generators from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        generators : list of |Generator| objects
        """

        cdef cnet.Gen** array
        cdef Generator gen

        array = <cnet.Gen**>malloc(len(generators)*sizeof(cnet.Gen*))
        for i in range(len(generators)):
            gen = generators[i]
            array[i] = gen._c_ptr

        cnet.NET_del_gens(self._c_net, array, len(generators))

        # Update pointers and alloc flags
        for i in range(len(generators)):
            gen = generators[i]
            if array[i] != NULL:
                gen._c_ptr = NULL
                gen.alloc = False

        # Clean up
        free(array)

    def add_loads(self, loads):
        """
        Adds loads to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        loads : list of |Load| objects
        """

        cdef cnet.Load** array
        cdef Load load

        old_num_loads = self.num_loads
        array = <cnet.Load**>malloc(len(loads)*sizeof(cnet.Load*))
        for i in range(len(loads)):
            load = loads[i]
            array[i] = load._c_ptr

        cnet.NET_add_loads(self._c_net, array, len(loads))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(loads)):
            load = loads[i]
            if array[i] != NULL:
                if load.alloc:
                    cload.LOAD_array_del(load._c_ptr,1)
                load._c_ptr = cnet.NET_get_load(self._c_net, old_num_loads+index)
                load.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_loads(self, loads):
        """
        Removes loads from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        loads : list of |Load| objects
        """

        cdef cnet.Load** array
        cdef Load load

        array = <cnet.Load**>malloc(len(loads)*sizeof(cnet.Load*))
        for i in range(len(loads)):
            load = loads[i]
            array[i] = load._c_ptr

        cnet.NET_del_loads(self._c_net, array, len(loads))

        # Update pointers and alloc flags
        for i in range(len(loads)):
            load = loads[i]
            if array[i] != NULL:
                load._c_ptr = NULL
                load.alloc = False

        # Clean up
        free(array)

    def add_shunts(self, shunts):
        """
        Adds shunts to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        shunts : list of |Shunt| objects
        """

        cdef cnet.Shunt** array
        cdef Shunt shunt

        old_num_shunts = self.num_shunts
        array = <cnet.Shunt**>malloc(len(shunts)*sizeof(cnet.Shunt*))
        for i in range(len(shunts)):
            shunt = shunts[i]
            array[i] = shunt._c_ptr

        cnet.NET_add_shunts(self._c_net, array, len(shunts))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(shunts)):
            shunt = shunts[i]
            if array[i] != NULL:
                if shunt.alloc:
                    cshunt.SHUNT_array_del(shunt._c_ptr,1)
                shunt._c_ptr = cnet.NET_get_shunt(self._c_net, old_num_shunts+index)
                shunt.alloc = False
                index += 1

        # Clean up
        free(array)
        
    def remove_shunts(self, shunts):
        """
        Removes shunts from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        shunts : list of |Shunt| objects
        """

        cdef cnet.Shunt** array
        cdef Shunt shunt

        array = <cnet.Shunt**>malloc(len(shunts)*sizeof(cnet.Shunt*))
        for i in range(len(shunts)):
            shunt = shunts[i]
            array[i] = shunt._c_ptr

        cnet.NET_del_shunts(self._c_net, array, len(shunts))

        # Update pointers and alloc flags
        for i in range(len(shunts)):
            shunt = shunts[i]
            if array[i] != NULL:
                shunt._c_ptr = NULL
                shunt.alloc = False

        # Clean up
        free(array)

    def add_batteries(self, batteries):
        """
        Adds batteries to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        batteries : list of |Battery| objects
        """

        cdef cnet.Bat** array
        cdef Battery bat

        old_num_bats = self.num_batteries
        array = <cnet.Bat**>malloc(len(batteries)*sizeof(cnet.Bat*))
        for i in range(len(batteries)):
            bat = batteries[i]
            array[i] = bat._c_ptr

        cnet.NET_add_bats(self._c_net, array, len(batteries))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(batteries)):
            bat = batteries[i]
            if array[i] != NULL:
                if bat.alloc:
                    cbat.BAT_array_del(bat._c_ptr,1)
                bat._c_ptr = cnet.NET_get_bat(self._c_net, old_num_bats+index)
                bat.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_batteries(self, batteries):
        """
        Removes batteries from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        batteries : list of |Battery| objects
        """

        cdef cnet.Bat** array
        cdef Battery bat

        array = <cnet.Bat**>malloc(len(batteries)*sizeof(cnet.Bat*))
        for i in range(len(batteries)):
            bat = batteries[i]
            array[i] = bat._c_ptr

        cnet.NET_del_bats(self._c_net, array, len(batteries))

        # Update pointers and alloc flags
        for i in range(len(batteries)):
            bat = batteries[i]
            if array[i] != NULL:
                bat._c_ptr = NULL
                bat.alloc = False

        # Clean up
        free(array)

    def add_var_generators(self, var_generators):
        """
        Adds var generators to the network.
        Flags are not preserved (for now).

        Parameters
        ----------
        var_generators : list of |VarGenerator| objects
        """

        cdef cnet.Vargen** array
        cdef VarGenerator gen

        old_num_gens = self.num_var_generators
        array = <cnet.Vargen**>malloc(len(var_generators)*sizeof(cnet.Vargen*))
        for i in range(len(var_generators)):
            gen = var_generators[i]
            array[i] = gen._c_ptr

        cnet.NET_add_vargens(self._c_net, array, len(var_generators))

        # Update pointers and alloc flags
        index = 0
        for i in range(len(var_generators)):
            gen = var_generators[i]
            if array[i] != NULL:
                if gen.alloc:
                    cvargen.VARGEN_array_del(gen._c_ptr,1)
                gen._c_ptr = cnet.NET_get_vargen(self._c_net, old_num_gens+index)
                gen.alloc = False
                index += 1

        # Clean up
        free(array)

    def remove_var_generators(self, var_generators):
        """
        Removes var generators from the network.
        All network flags are cleared (for now). 

        Parameters
        ----------
        var_generators : list of |VarGenerator| objects
        """

        cdef cnet.Vargen** array
        cdef VarGenerator gen

        array = <cnet.Vargen**>malloc(len(var_generators)*sizeof(cnet.Vargen*))
        for i in range(len(var_generators)):
            gen = var_generators[i]
            array[i] = gen._c_ptr

        cnet.NET_del_vargens(self._c_net, array, len(var_generators))

        # Update pointers and alloc flags
        for i in range(len(var_generators)):
            gen = var_generators[i]
            if array[i] != NULL:
                gen._c_ptr = NULL
                gen.alloc = False

        # Clean up
        free(array)
            
    def add_var_generators_from_parameters(self, buses, power_capacity, power_base, power_std=0., corr_radius=0, corr_value=0.):
        """
        Adds variable generators to the network using the given parameters.
        The capacities of the generators are divided evenly.

        Parameters
        ----------
        buses : list of |Bus| objects
        power_capacity : percentage of max total load power (float)
        power_base : percentage of power capacity (float)
        power_std : percentage of power capacity (float)
        corr_radius : number of branches for correlation radius (int)
        corr_value : correlation coefficient for correlated generators (float)
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
            cnet.NET_add_vargens_from_params(self._c_net,head._c_ptr,power_capacity,power_base,power_std,corr_radius,corr_value)
        else:
            cnet.NET_add_vargens_from_params(self._c_net,NULL,power_capacity,power_base,power_std,corr_radius,corr_value)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def add_batteries_from_parameters(self, buses, power_capacity, energy_capacity, eta_c=1., etc_d=1.):
        """
        Adds batteries to the network using the given parameters. 
        The power and energy capacities of the batteries are divided evenly.

        Parameters
        ----------
        buses : list of |Bus| objects
        power_capacity : percentage of max total load power (float)
        energy_capacity : percentage of max total load energy during one time period (float)
        eta_c : charging efficiency in (0,1] (float)
        eta_d : discharging efficiency in (0,1] (float)
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
            cnet.NET_add_batteries_from_params(self._c_net,head._c_ptr,power_capacity,energy_capacity,eta_c,etc_d)
        else:
            cnet.NET_add_batteries_from_params(self._c_net,NULL,power_capacity,energy_capacity,eta_c,etc_d)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

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

    def copy_from_network(self, net):
        """
        Copies data from another network.

        Parameters
        ----------
        net : |Network|
        """
        
        cdef Network n = net
        if net is not None:
            cnet.NET_copy_from_net(self._c_net, n._c_net)

    def create_var_generators_P_sigma(self, spread, corr):
        """
        Creates covariance matrix (lower triangular part) for
        active powers of variable generators.

        Parameters
        ----------
        spead : correlation neighborhood in terms of number of edges (int)
        corr : correlation coefficient for neighboring generators (float)

        Returns
        -------
        sigma : |CooMatrix|
        """

        sigma = Matrix(cnet.NET_create_vargen_P_sigma(self._c_net,spread,corr),
                       owndata=True)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))
        else:
            return sigma

    def get_copy(self):
        """ 
        Gets deep copy of network.

        Returns
        -------
        net : |Network|
        """

        cdef Network net = new_Network(cnet.NET_get_copy(self._c_net))
        net.alloc = True
        return net

    def get_var_info_string(self, index):
        """
        Gets info string of variable associated with index.
        The info string has format ``obj_type:obj_index:quantity:time``.

        Parameters
        ----------
        index : int

        Returns
        -------
        info : string
        """

        cdef char* info_string = cnet.NET_get_var_info_string(self._c_net, index)
        if info_string:
            s = info_string.decode('UTF-8')
            free(info_string)
            return s
        else:
            raise NetworkError('index does not correspond to any variable')

    def get_bus_from_number(self, number):
        """
        Gets bus with the given number.

        Parameters
        ----------
        number : int

        Returns
        -------
        bus : |Bus|
        """

        ptr = cnet.NET_bus_hash_number_find(self._c_net,number)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_bus_from_name(self, name):
        """
        Gets bus with the given name.

        Parameters
        ----------
        name : string

        Returns
        -------
        bus : |Bus|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_bus_hash_name_find(self._c_net,name)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_bus(self, index):
        """
        Gets bus with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bus : |Bus|
        """

        ptr = cnet.NET_get_bus(self._c_net,index)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('invalid bus index')

    def get_branch(self, index):
        """
        Gets branch with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        branch : |Branch|
        """

        ptr = cnet.NET_get_branch(self._c_net,index)
        if ptr is not NULL:
            return new_Branch(ptr)
        else:
            raise NetworkError('invalid branch index')

    def get_generator(self,index):
        """
        Gets generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : |Generator|
        """

        ptr = cnet.NET_get_gen(self._c_net,index)
        if ptr is not NULL:
            return new_Generator(ptr)
        else:
            raise NetworkError('invalid gen index')

    def get_shunt(self, index):
        """
        Gets shunt with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        shunt : |Shunt|
        """

        ptr = cnet.NET_get_shunt(self._c_net,index)
        if ptr is not NULL:
            return new_Shunt(ptr)
        else:
            raise NetworkError('invalid shunt index')

    def get_load(self, index):
        """
        Gets load with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        load : |Load|
        """

        ptr = cnet.NET_get_load(self._c_net,index)
        if ptr is not NULL:
            return new_Load(ptr)
        else:
            raise NetworkError('invalid load index')

    def get_var_generator(self, index):
        """
        Gets variable generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        vargen : |VarGenerator|
        """

        ptr = cnet.NET_get_vargen(self._c_net,index)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('invalid vargen index')

    def get_battery(self, index):
        """
        Gets battery with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bat : |Battery|
        """

        ptr = cnet.NET_get_bat(self._c_net,index)
        if ptr is not NULL:
            return new_Battery(ptr)
        else:
            raise NetworkError('invalid battery index')

    def get_generator_from_name_and_bus_number(self, name, number):
        """
        Gets generator of given name connected to the bus of the 
        given number.

        Parameters
        ----------
        name : string
        number : integer

        Returns
        -------
        gen : |Generator|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_gen_from_name_and_bus_number(self._c_net, name, number)
        if ptr is not NULL:
            return new_Generator(ptr)
        else:
            raise NetworkError('generator not found')

    def get_branch_from_name_and_bus_numbers(self, name, number1, number2):
        """
        Gets branch of given name connected across buss of the 
        given numbers.

        Parameters
        ----------
        name : string
        number1 : integer
        number2 : integer

        Returns
        -------
        branch : |Branch|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_branch_from_name_and_bus_numbers(self._c_net, name, number1, number2)
        if ptr is not NULL:
            return new_Branch(ptr)
        else:
            raise NetworkError('branch not found')

    def get_load_from_name_and_bus_number(self, name, number):
        """
        Gets load of given name connected to the bus of the 
        given number.

        Parameters
        ----------
        name : string
        number : integer

        Returns
        -------
        load : |Load|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_load_from_name_and_bus_number(self._c_net, name, number)
        if ptr is not NULL:
            return new_Load(ptr)
        else:
            raise NetworkError('load not found')

    def get_shunt_from_name_and_bus_number(self, name, number):
        """
        Gets shunt of given name connected to the bus of the 
        given number.

        Parameters
        ----------
        name : string
        number : integer

        Returns
        -------
        shunt : |Shunt|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_shunt_from_name_and_bus_number(self._c_net, name, number)
        if ptr is not NULL:
            return new_Shunt(ptr)
        else:
            raise NetworkError('shunt not found')

    def get_var_generator_from_name_and_bus_number(self, name, number):
        """
        Gets variable generator of given name connected to the bus of the 
        given number.

        Parameters
        ----------
        name : string
        number : integer

        Returns
        -------
        vargen : |VarGenerator|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_vargen_from_name_and_bus_number(self._c_net, name, number)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('variable generator not found')

    def get_battery_from_name_and_bus_number(self, name, number):
        """
        Gets battery of given name connected to the bus of the 
        given number.

        Parameters
        ----------
        name : string
        number : integer

        Returns
        -------
        bat : |Battery|
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_get_bat_from_name_and_bus_number(self._c_net, name, number)
        if ptr is not NULL:
            return new_Battery(ptr)
        else:
            raise NetworkError('battery not found')
        
    def get_generator_buses(self):
        """
        Gets list of buses where generators are connected.

        Returns
        -------
        buses : list of |Bus| objects
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
        buses : list of |Bus| objects
        """

        buses = []
        cdef cbus.Bus* b = cnet.NET_get_load_buses(self._c_net)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def get_var_values(self, option='current'):
        """
        Gets network variable values.

        Parameters
        ----------
        option : string (|RefVarValueOptions|)

        Returns
        -------
        values : |Array|
        """
        return Vector(cnet.NET_get_var_values(self._c_net, str2const[option]), owndata=True)

    def get_var_projection(self, obj_type, props, q, t_start=0, t_end=None):
        """
        Gets projection matrix for specific object variables.

        Parameters
        ----------
        obj_type : string (|RefObjects|)
        props : string or list of strings (|RefProperties|)
        q : string or list of strings (|RefQuantities|)
        t_start : int
        t_end : int (inclusive)

        Returns
        -------
        P : |CooMatrix|
        """

        props = props if isinstance(props,list) else [props]
        q = q if isinstance(q,list) else [q]

        if t_end is None:
            t_end = self.num_periods-1
        m = Matrix(cnet.NET_get_var_projection(self._c_net,
                                               str2obj[obj_type],
					                           reduce(lambda x,y: x|y,[str2prop[obj_type][pp] for pp in props],0),
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

    def get_num_buses_reg_by_tran(self, only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by tap-changing transformers.

        Parameters
        ----------
        only : |TrueFalse|

        Returns
        -------
        num : int
        """

        if not only:
            return cnet.NET_get_num_buses_reg_by_tran(self._c_net)
        else:
            return cnet.NET_get_num_buses_reg_by_tran_only(self._c_net)

    def get_num_buses_reg_by_shunt(self, only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by switched shunt devices.

        Parameters
        ----------
        only : |TrueFalse|

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

    def get_num_batteries(self):
        """
        Gets number of batteries in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_bats(self._c_net)

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

    def has_same_ptr(self, Network other):
        """
        Checks whether network shares memory with another network.

        Parameters
        ----------
        other : |Network|

        Returns
        -------
        flag : |TrueFalse|
        """

        return self._c_net == other._c_net

    def has_error(self):
        """
        Indicates whether the network has the error flag set due to an
        invalid operation.

        Returns
        -------
        flag : |TrueFalse|
        """

        return cnet.NET_has_error(self._c_net)

    def set_flags(self, obj_type, flags, props, q):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj_type : string (|RefObjects|)
        flags : string or list of strings (|RefFlags|)
        props : string or list of strings (|RefProperties|)
        q : string or list of strings (|RefQuantities|)
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

    def set_flags_of_component(self, obj, flags, q):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj : |Bus|, |Branch|, |Generator|, |Load|, |Shunt|, |VarGenerator|, |Battery|
        flags : string or list of strings (|RefFlags|)
        q : string or list of strings (|RefQuantities|)
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

    def set_branch_array(self, size):
        """
        Allocates and sets branch array.

        Parameters
        ----------
        size : int
        """

        cdef cbranch.Branch* array = cbranch.BRANCH_array_new(size,self.num_periods)
        cnet.NET_set_branch_array(self._c_net,array,size)

    def set_bus_array(self, size):
        """
        Allocates and sets bus array.

        Parameters
        ----------
        size : int
        """

        cdef cbus.Bus* array = cbus.BUS_array_new(size,self.num_periods)
        cnet.NET_set_bus_array(self._c_net,array,size)

    def set_gen_array(self, size):
        """
        Allocates and sets generator array.

        Parameters
        ----------
        size : int
        """

        cdef cgen.Gen* array = cgen.GEN_array_new(size,self.num_periods)
        cnet.NET_set_gen_array(self._c_net,array,size)

    def set_load_array(self, size):
        """
        Allocates and sets load array.

        Parameters
        ----------
        size : int
        """

        cdef cload.Load* array = cload.LOAD_array_new(size,self.num_periods)
        cnet.NET_set_load_array(self._c_net,array,size)

    def set_shunt_array(self, size):
        """
        Allocates and sets shunt array.

        Parameters
        ----------
        size : int
        """

        cdef cshunt.Shunt* array = cshunt.SHUNT_array_new(size,self.num_periods)
        cnet.NET_set_shunt_array(self._c_net,array,size)

    def set_vargen_array(self, size):
        """
        Allocates and sets variable generator array.

        Parameters
        ----------
        size : int
        """

        cdef cvargen.Vargen* array = cvargen.VARGEN_array_new(size,self.num_periods)
        cnet.NET_set_vargen_array(self._c_net,array,size)

    def set_battery_array(self, size):
        """
        Allocates and sets battery array.

        Parameters
        ----------
        size : int
        """

        cdef cbat.Bat* array = cbat.BAT_array_new(size,self.num_periods)
        cnet.NET_set_bat_array(self._c_net,array,size)

    def set_var_values(self, values):
        """
        Sets network variable values.

        Parameters
        ----------
        values : |Array|
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cnet.REAL*>(x.data),x.size)
        cnet.NET_set_var_values(self._c_net,v)
        free(v)

    def show_components(self):
        """
        Shows information about the number of network components of each type.
        """

        print(cnet.NET_get_show_components_str(self._c_net).decode('UTF-8'))

    def show_properties(self, t=0):
        """
        Shows information about the state of the network component quantities.

        Parameters
        ----------
        t : time period (int)
        """

        print(cnet.NET_get_show_properties_str(self._c_net,t).decode('UTF-8'))

    def update_properties(self, values=None):
        """
        Re-computes the network properties using the given values
        of the network variables. If no values are given, then the
        current values of the network variables are used.

        Parameters
        ----------
        values : |Array|
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(<cnet.REAL*>(x.data),x.size) if values is not None else NULL
        cnet.NET_update_properties(self._c_net,v)
        if v != NULL:
            free(v)

    def propogate_data_in_time(self, start, end):
        """ 
        Propogates data from the first period through time.

        Parameters
        ----------
        start : int
        end : int
        """
        cnet.NET_propagate_data_in_time(self._c_net, start, end)

    def update_set_points(self):
        """
        Updates voltage magnitude set points of gen-regulated buses
        to be equal to the bus voltage magnitudes.
        """

        cnet.NET_update_set_points(self._c_net)

    def update_hashes(self):
        """
        Updates the bus name and number hash lists.
        """
        
        cdef Bus cb

        for bus in self.buses:
            cb = bus
            cnet.NET_bus_hash_number_add(self._c_net,cb._c_ptr)
            cnet.NET_bus_hash_name_add(self._c_net,cb._c_ptr)

    property error_string:
        """ Error string (string). """
        def __get__(self):
            return cnet.NET_get_error_string(self._c_net).decode('UTF-8')

    property json_string:
        """ JSON string (string). """
        def __get__(self):
            cdef char* json_string = cnet.NET_get_json_string(self._c_net)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cnet.NET_get_num_periods(self._c_net)

    property base_power:
        """ System base power (MVA) (float). """
        def __get__(self): return cnet.NET_get_base_power(self._c_net)
        def __set__(self,v): cnet.NET_set_base_power(self._c_net, v)

    property buses:
        """ List of |Bus| objects. """
        def __get__(self):
            return [self.get_bus(i) for i in range(self.num_buses)]

    property branches:
        """ List of |Branch| objects. """
        def __get__(self):
            return [self.get_branch(i) for i in range(self.num_branches)]

    property generators:
        """ List of |Generator| objects. """
        def __get__(self):
            return [self.get_generator(i) for i in range(self.num_generators)]

    property shunts:
        """ List of |Shunt| objects. """
        def __get__(self):
            return [self.get_shunt(i) for i in range(self.num_shunts)]

    property loads:
        """ List of |Load| objects. """
        def __get__(self):
            return [self.get_load(i) for i in range(self.num_loads)]

    property var_generators:
        """ List of |VarGenerator| objects. """
        def __get__(self):
            return [self.get_var_generator(i) for i in range(self.num_var_generators)]

    property batteries:
        """ List of |Battery| objects. """
        def __get__(self):
            return [self.get_battery(i) for i in range(self.num_batteries)]

    property num_buses:
        """ Number of buses in the network (int). """
        def __get__(self): return cnet.NET_get_num_buses(self._c_net)

    property num_branches:
        """ Number of branches in the network (int). """
        def __get__(self): return cnet.NET_get_num_branches(self._c_net)

    property num_generators:
        """ Number of generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_gens(self._c_net)

    property num_loads:
        """ Number of loads in the network (int). """
        def __get__(self): return cnet.NET_get_num_loads(self._c_net)

    property num_shunts:
        """ Number of shunt devices in the network (int). """
        def __get__(self): return cnet.NET_get_num_shunts(self._c_net)

    property num_var_generators:
        """ Number of variable generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_vargens(self._c_net)

    property num_batteries:
        """ Number of batteries in the network (int). """
        def __get__(self): return cnet.NET_get_num_bats(self._c_net)

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

    property total_load_P:
        """ Total load active power (MW) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_total_load_P(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_v_max:
        """ Maximum bus voltage magnitude (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_max(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_v_min:
        """ Minimum bus voltage magnitude (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_min(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_v_vio:
        """ Maximum bus voltage magnitude limit violation (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_bus_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_P_mis:
        """ Largest bus active power mismatch in the network (MW) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_bus_P_mis(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property bus_Q_mis:
        """ Largest bus reactive power mismatch in the network (MVAr) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_bus_Q_mis(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_P_cost:
        """ Total active power generation cost ($/hr) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_gen_P_cost(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_v_dev:
        """ Largest voltage magnitude deviation from set point of bus regulated by generator (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_gen_v_dev(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_Q_vio:
        """ Largest generator reactive power limit violation (MVAr) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_gen_Q_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property gen_P_vio:
        """ Largest generator active power limit violation (MW) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_gen_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by transformer (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_tran_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_r_vio:
        """ Largest transformer tap ratio limit violation (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_tran_r_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property tran_p_vio:
        """ Largest transformer phase shift limit violation (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_tran_p_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property shunt_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by switched shunt device (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_shunt_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property shunt_b_vio:
        """ Largest switched shunt susceptance limit violation (p.u.) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_shunt_b_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property load_P_util:
        """ Total active power consumption utility ($/hr) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_load_P_util(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property load_P_vio:
        """ Largest load active power limit violation (MW) (float or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_load_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeFloat(r[0])
            else:
                return np.array(r)

    property num_actions:
        """ Number of control adjustments (int or |Array|). """
        def __get__(self):
            r = [cnet.NET_get_num_actions(self._c_net,t) for t in range(self.num_periods)]
            if self.num_periods == 1:
                return AttributeInt(r[0])
            else:
                return np.array(r)

    property var_generators_corr_radius:
        """ Correlation radius of variable generators (number of edges). """
        def __get__(self): return cnet.NET_get_vargen_corr_radius(self._c_net)

    property var_generators_corr_value:
        """ Correlation value (coefficient) of variable generators. """
        def __get__(self): return cnet.NET_get_vargen_corr_value(self._c_net)

    property show_components_str:
        """ String with information about network components. """
        def __get__(self): return cnet.NET_get_show_components_str(self._c_net).decode('UTF-8')

cdef public new_Network(cnet.Net* n):
    if n is not NULL:
        net = Network(alloc=False)
        net._c_net = n
        return net
    else:
        raise NetworkError('no network data')
