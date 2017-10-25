.. include:: defs.hrst

.. _net:

**************
Power Networks
**************

This section describes how to use and analyze power networks in PFNET.

.. _net_overview:

Overview
========

Power networks in PFNET are represented by objects of type |Network|. These objects are created from power network data files using a :ref:`Parser <parsers>`, as described in the previous section. Once the network is created, it can be analyzed, visualized, and used to construct network optimization problems.

An important attribute of the |Network| class is :data:`base_power <pfnet.Network.base_power>`. This quantity, which has units of MVA, is useful for converting power quantities in per unit system base power to MW or MVAr.

.. _net_components:

Components
==========

Power networks have several components. These are :ref:`buses <net_bus>`, :ref:`branches <net_branch>`, :ref:`generators <net_gen>`, :ref:`shunt devices <net_shunt>`, :ref:`loads <net_load>`, :ref:`variable generators <net_vargen>`, *e.g.*, renewable energy sources, and :ref:`batteries <net_bat>`. For obtaining an overview of the components that form a network, the class method :func:`show_components() <pfnet.Network.show_components>` can be used, as illustrated in the following example::

  >>> import pfnet

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')
  >>> net.show_components()

  Network Components
  ------------------
  buses            : 14
    slack          : 1
    reg by gen     : 5
    reg by tran    : 0
    reg by shunt   : 0
  shunts           : 1
    fixed          : 1
    switched v     : 0
  branches         : 20
    lines          : 17
    fixed trans    : 3
    phase shifters : 0
    tap changers v : 0
    tap changers Q : 0
  generators       : 5
    slack          : 1
    reg            : 5
    P adjust       : 5
  loads            : 11
    P adjust       : 0
  vargens          : 0
  batteries        : 0

Again, this and subsequent examples assume that the Python interpreter was started from a directory that contains the sample case |ieee14|.

.. _net_bus:

Buses
-----

Buses in a power network are objects of type |Bus|. Each bus has an :data:`index <pfnet.Bus.index>`, a :data:`number <pfnet.Bus.number>`, and a :data:`name <pfnet.Bus.name>` attribute that can be used to identify this bus in a network. The :data:`index <pfnet.Bus.index>` is associated with the location of the bus in the underlying C array of bus structures, while the :data:`number <pfnet.Bus.number>` and :data:`name <pfnet.Bus.name>` attributes are specified in the input data. An :data:`index <pfnet.Bus.index>`, a :data:`number <pfnet.Bus.number>`, or a :data:`name <pfnet.Bus.name>` can be used to extract a specific bus from a network using the |Network| class methods :func:`get_bus() <pfnet.Network.get_bus>`, :func:`get_bus_from_number() <pfnet.Network.get_bus_from_number>`, and :func:`get_bus_from_name() <pfnet.Network.get_bus_from_name>`, respectively::

  >>> bus = net.get_bus(10)

  >>> print bus.index == 10
  True

  >>> other_bus = net.get_bus_from_number(bus.number)

  >>> print bus == other_bus
  True

For convenience, a list of all the buses in the network is contained in the :data:`buses <pfnet.Network.buses>` attribute of the |Network| class.

Buses in a network can have different properties. For example, some buses can be slack buses and others can have their voltage magnitudes regulated by generators, tap-changing transformers, or switched shunt devices. The |Bus| class provides methods for checking whether a bus has specific properties. The following example shows how to get a list of all the buses whose voltage magnitudes are regulated by generators::

  >>> reg_buses = [bus for bus in net.buses if bus.is_regulated_by_gen()]

  >>> print len(reg_buses), net.get_num_buses_reg_by_gen()
  5 5

A bus also has information about the devices that are connected to it or that are regulating its voltage magnitude. For example, the attributes :data:`generators <pfnet.Bus.generators>` and :data:`reg_trans <pfnet.Bus.reg_trans>` contain a list of generators connected to the bus and a list of tap-changing transformers regulating its voltage magnitude, respectively.

More information about network buses can be found in the :ref:`API reference <ref_bus>`.

.. _net_branch:

Branches
--------

Branches in a power network are objects of type |Branch| and are represented mathematically by the model described in Section 2.1.2 of [TT2015]_. Each branch has an :data:`index <pfnet.Branch.index>` and a :data:`name <pfnet.Branch.name>` attribute that can be used to identify this branch in a network. The |Network| class methods :func:`get_branch() <pfnet.Network.get_branch>` and :func:`get_branch_from_name_and_bus_numbers() <pfnet.Network.get_branch_from_name_and_bus_numbers>` can be used to get a specific branch of the network::

  >>> branch = net.get_branch(5)

  >>> print branch.index == 5
  True

For convenience, a list of all the branches in the network is contained in the :data:`branches <pfnet.Network.branches>` attribute of the |Network| class.

Branches in a power network can have different properties. For example, some branches can be transmission lines, fixed transformers, tap-changing transformers, or phase-shifting transformers. Tap-changing transformers in turn can control the reactive power flowing through the branch or the voltage magnitude of a bus. The |Branch| class provides methods for checking whether a branch has specific properties. The following example shows how to get a list of all the branches that are transmission lines::

  >>> lines = [br for br in net.branches if br.is_line()]

  >>> print len(lines), net.get_num_lines()
  17 17

For branches that are transformers, the |Branch| class attributes :data:`ratio <pfnet.Branch.ratio>` and :data:`phase <pfnet.Branch.phase>` correspond to the transformer's tap ratio and phase shift, respectively. These attributes correspond to the quantities :math:`a_{km}` and :math:`\phi_{km}` of the branch model described in Section 2.1.2 of [TT2015]_. The quantity :math:`a_{mk}` in this model is always one.

More information about network branches can be found in the :ref:`API reference <ref_branch>`. 

.. _net_gen:

Generators
----------

Generators in a power network are objects of type |Generator|. Each generator has an :data:`index <pfnet.Generator.index>` and a :data:`name <pfnet.Generator.name>` attribute that can be used to identify this generator in a network. The |Network| class methods :func:`get_generator() <pfnet.Network.get_generator>` and :func:`get_generator_from_name_and_bus_number() <pfnet.Network.get_generator_from_name_and_bus_number>` can be used to get a specific generator of the network::

  >>> gen = net.get_generator(2)

  >>> print gen.index == 2
  True

For convenience, a list of all the generators in the network is contained in the :data:`generators <pfnet.Network.generators>` attribute of the |Network| class.

Generators in a power network can also have different properties. For example, some generators can be slack generators and others can provide bus voltage magnitude regulation. The |Generator| class provides methods for checking whether a generator has specific properties. The following example shows how to get a list of all the slack generators::

  >>> slack_gens = [g for g in net.generators if g.is_slack()]

  >>> print len(slack_gens), net.get_num_slack_gens()
  1 1

The active and reactive powers that a generator injects into the bus to which it is connected are obtained from the :data:`P <pfnet.Generator.P>` and :data:`Q <pfnet.Generator.Q>` attributes of the  |Generator| class. These quantities are given in units of per unit :data:`system base power <pfnet.Network.base_power>`. The following example computes the total active power injected into the network by generators in units of MW::

  >>> print sum([g.P for g in net.generators])*net.base_power
  272.4

More information about network generators can be found in the :ref:`API reference <ref_gen>`. 
  
.. _net_shunt:

Shunt Devices
-------------

Shunt devices in a power network are objects of type |Shunt|. Each shunt has an :data:`index <pfnet.Shunt.index>` and a :data:`name <pfnet.Shunt.name>` attribute that can be used to identify this shunt in a network. The |Network| class methods :func:`get_shunt() <pfnet.Network.get_shunt>` and :func:`get_shunt_from_name_and_bus_number() <pfnet.Network.get_shunt_from_name_and_bus_number>` can be used to get a specific shunt of the network::

  >>> shunt = net.get_shunt(0)

  >>> print shunt.index == 0
  True

For convenience, a list of all the shunt devices in the network is contained in the :data:`shunts <pfnet.Network.shunts>` attribute of the |Network| class.

As with other network components, shunt devices can have different properties. Some shunt devices can be fixed while others can be switchable and configured to regulate a bus voltage magnitude.

More information about network shunts can be found in the :ref:`API reference <ref_shunt>`.

.. _net_load:

Loads
-----

Loads in a power network are objects of type |Load|. As with other components, the :data:`index <pfnet.Load.index>` and :data:`name <pfnet.Load.name>` attributes can be used to identify a load in the network. A list of all the loads in the network is contained in the :data:`loads <pfnet.Network.loads>` attribute of the |Network| class.

As with generators, the active and reactive powers that a load consumes from the bus to which it is connected are obtained from the :data:`P <pfnet.Load.P>` and :data:`Q <pfnet.Load.Q>` attributes of the |Load| class. They are also given in units of per unit :data:`system base power <pfnet.Network.base_power>`.

More information about network loads can be found in the :ref:`API reference <ref_load>`.

.. _net_vargen:

Variable Generators
-------------------

Variable generators in a power network are objects of type |VarGenerator|. They represent non-dispatchable energy sources such as wind generators or farms and photovoltaic power plants. As with other components, the :data:`index <pfnet.VarGenerator.index>` and :data:`name <pfnet.VarGenerator.name>` attributes can be used to identify a variable generator in the network. Also, a list of all the variable generators in the network is contained in the :data:`var_generators <pfnet.Network.var_generators>` attribute of the |Network| class.

As with generators and loads, the active and reactive powers produced by a variable generator are obtained from the :data:`P <pfnet.VarGenerator.P>` and :data:`Q <pfnet.VarGenerator.Q>` attributes of the |VarGenerator| class in units of per unit :data:`system base power <pfnet.Network.base_power>`. Output limits of a variable generator are given by the attributes :data:`P_min <pfnet.VarGenerator.P_min>`, :data:`P_max <pfnet.VarGenerator.P_max>`, :data:`Q_min <pfnet.VarGenerator.Q_min>`, and :data:`Q_max <pfnet.VarGenerator.Q_max>`.

The output of variable generators in a network is subject to random variations that can be correlated, especially for devices that are "nearby". The method :func:`create_var_generators_P_sigma() <pfnet.Network.create_var_generators_P_sigma>` of the |Network| class allows constructing a covariance matrix for these variations based on a "correlation distance" ``N`` and a given correlation coefficient. The cross-covariance between the variation of any two devices that are connected to buses that are at most ``N`` branches away from each other is set such that it is consistent with the given correlation coefficient. For other devices, the cross-covariance is set to zero.  

Lastly, since many power network input files do not have variable generator information, these devices can be manually added to a network using the :func:`add_var_generators() <pfnet.Network.add_var_generators>` method of the |Network| class.

More information about network variable generators can be found in the :ref:`API reference <ref_vargen>`.

.. _net_bat:

Batteries
---------

Batteries are objects of type |Battery| and have an :data:`index <pfnet.Battery.index>` and :data:`name <pfnet.Battery.name>` attribute like all the other network components. Other important attributes of these objects are energy level :data:`E <pfnet.Battery.E>` and charging power :data:`P <pfnet.Battery.P>`.  Since power network input files do not have variable generator information, these devices can be manually added to a network using the :func:`add_batteries() <pfnet.Network.add_batteries>` method of the :class:`Network <pfnet.Network>` class.

More information about network batteries can be found in the :ref:`API reference <ref_bat>`.

.. _net_properties:

Properties
==========

A |Network| object has several quantities or ``properties`` that provide important information about the state of the network. The following table provides a description of each of these properties.

================== ================================================================= ========
Names              Description                                                       Units
================== ================================================================= ========
``bus_v_max``      Maximum bus voltage magnitude                                     per unit
``bus_v_min``      Minimum bus voltage magnitude                                     per unit
``bus_v_vio``      Maximum bus voltage magnitude limit violation                     per unit
``bus_P_mis``      Maximum absolute bus active power mismatch                        MW
``bus_Q_mis``      Maximum absolute bus reactive power mismatch                      MVAr
``gen_P_cost``     Total active power generation cost                                $/hour
``gen_v_dev``      Maximum set point deviation of generator-regulated voltage        per unit
``gen_Q_vio``      Maximum generator reactive power limit violation                  MVAr
``gen_P_vio``      Maximum generator active power limit violation                    MW
``tran_v_vio``     Maximum band violation of transformer-regulated voltage           per unit
``tran_r_vio``     Maximum tap ratio limit violation of tap-changing transformer     unitless
``tran_p_vio``     Maximum phase shift limit violation of phase-shifting transformer radians
``shunt_v_vio``    Maximum band violation of shunt-regulated voltage                 per unit
``shunt_b_vio``    Maximum susceptance limit violation of switched shunt device      per unit
``load_P_util``    Total active power consumption utility                            $/hour
``load_P_vio``     Maximum load active power limit violation                         MW
``num_actions``    Number of control adjustments (greater than 2% of control range)  unitless
================== ================================================================= ========

All of these properties are attributes of the |Network| class. If there is a change in the network, *e.g.*, the voltage magnitude :data:`v_mag <pfnet.Bus.v_mag>` of a bus is changed, the class method :func:`update_properties() <pfnet.Network.update_properties>` needs to be called in order for the network properties to reflect the change. The following example shows how to update and extract properties::

  >>> print net.bus_v_max
  1.09

  >>> for bus in net.buses:
  ...     bus.v_mag = bus.v_mag + 0.1
  ...

  >>> print net.bus_v_max
  1.09

  >>> net.update_properties()

  >>> print net.bus_v_max
  1.19

For convenience, all the network properties can be extracted at once in a dictionary using the :func:`get_properties() <pfnet.Network.get_properties>` class method::

  >>> properties = net.get_properties()

  >>> print properties['bus_v_max']
  1.19

.. _net_variables:

Variables
=========

Network quantities can be specified to be ``variables``. This is useful to represent network quantities with vectors and move to the linear algebra domain to do some computations.

To set network quantities as variables, the |Network| class method :func:`set_flags() <pfnet.Network.set_flags>` is used. This method takes as arguments a :ref:`component name <ref_net_obj>`, one or more :ref:`flag names <ref_net_flag>`, one or more ``component properties``, and one or more ``component quantities``.

**Component properties** are component-specific. They can be combined into a list to make properties more complex and target a specific subset of components of a given type. More information can be found in the following sections:

* :ref:`ref_bus_prop`
* :ref:`ref_branch_prop`
* :ref:`ref_gen_prop`
* :ref:`ref_load_prop`
* :ref:`ref_shunt_prop`
* :ref:`ref_vargen_prop`

**Component quantities** are also component-specific. They can be combined into a list to specify all quantities that should be affected by the method :func:`set_flags() <pfnet.Network.set_flags>`. More information can be found in the following sections:

* :ref:`ref_bus_q`
* :ref:`ref_branch_q`
* :ref:`ref_gen_q`
* :ref:`ref_load_q`
* :ref:`ref_shunt_q`
* :ref:`ref_vargen_q`

The following example shows how to set as variables all the voltage magnitudes and angles of buses regulated by generators::

  >>> import pfnet

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> print net.num_vars
  0

  >>> net.set_flags('bus',
  ...               'variable',
  ...               'regulated by generator',
  ...               ['voltage magnitude', 'voltage angle'])

  >>> print net.num_vars, 2*net.get_num_buses_reg_by_gen()
  10 10

Network components have a :func:`has_flags() <pfnet.Bus.has_flags>` method that allows checking whether flags of a certain type associated with specific quantities are set.

Once variables have been set, the vector containing all the current variable values can be extracted using :func:`get_var_values() <pfnet.Network.get_var_values>`::

  >>> values = net.get_var_values()

  >>> print type(values)
  <type 'numpy.ndarray'>

  >>> print values.shape
  (10,)

The network components that have quantities set as variables have indices that can be used to locate these quantities in the vector of all variable values::

  >>> bus = [bus for bus in net.buses if bus.is_reg_by_gen()][0]

  >>> print bus.has_flags('variable','voltage magnitude')
  True

  >>> bus.has_flags('variable','voltage angle')
  True

  >>> print bus.v_mag, net.get_var_values()[bus.index_v_mag]
  1.09 1.09

  >>> print bus.v_ang, net.get_var_values()[bus.index_v_ang]
  -0.23 -0.23

A vector of variable values can be used to update the corresponding network quantities. This is done with the |Network| class method :func:`set_var_values() <pfnet.Network.set_var_values>`::

  >>> bus.has_flags('variable','voltage angle')
  True

  >>> values = net.get_var_values()

  >>> print bus.v_mag
  1.09

  >>> values[bus.index_v_mag] = 1.20
  >>> net.set_var_values(values)

  >>> print bus.v_mag
  1.20

As will be seen later, variables are also useful for constructing network optimization problems.

The class method :func:`get_var_values() <pfnet.Network.get_var_values>` can also be used to get upper or lower limits of the variables. To do this, a valid :ref:`variable value option <ref_var_values>` must be passed to this method.

In addition to the class method :func:`set_flags() <pfnet.Network.set_flags>`, which allows specifying variables of components having certain properties, one can also use the |Network| class method :func:`set_flags_of_component() <pfnet.Network.set_flags_of_component>` to specify variables of individual components. This is useful when the desired components cannot be targeted using the available ``component properties``. For example, the following code illustrates how to set as variables the voltage magnitudes of buses whose indices are multiples of three::

  >>> net.clear_flags()

  >>> for bus in net.buses:
  ...     if bus.index % 3 == 0:
  ...         net.set_flags_of_component(bus,'variable','voltage magnitude')

  >>> print net.num_vars, len([bus for bus in net.buses if bus.index % 3 == 0]), net.num_buses
  5 5 14

Lastly, a very useful method of the |Network| class is the method :func:`get_var_info_string() <pfnet.Network.get_var_info_string>`. This method takes as argument an index of the vector of all variable values, and returns a string with information about the network quantity associated with it. 
  
.. _net_var_projections:

Projections
===========

As explained above, once the network variables have been set, a vector with the current values of the selected variables is obtained with the class method :func:`get_var_values() <pfnet.Network.get_var_values>`. To extract subvectors that contain values of specific variables, projection matrices can be used. These matrices can be obtained using the class method :func:`get_var_projection() <pfnet.Network.get_var_projection>`, which takes as arguments a :ref:`component name <ref_net_obj>`, one or more ``component properties``, *e.g.*, :ref:`ref_bus_prop`, and one or more ``component quantities``, *e.g.*, :ref:`ref_bus_q`. The next example sets the variables of the network to be the bus voltage magnitudes and angles of all the buses, extracts the vector of values of all variables, and then extracts two subvectors having only voltage magnitudes and only voltage angles, respectively::

  >>> import pfnet
  >>> import numpy as np

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> net.set_flags('bus',
  ...               'variable',
  ...               'any',
  ...               ['voltage magnitude','voltage angle'])

  >>> print net.num_vars, 2*net.num_buses
  28 28

  >>> P1 = net.get_var_projection('bus', 'any', 'voltage magnitude')
  >>> P2 = net.get_var_projection('bus', 'any', 'voltage angle')

  >>> print type(P1)
  <class 'scipy.sparse.coo.coo_matrix'>

  >>> x = net.get_var_values()
  >>> v_mags = P1*x
  >>> v_angs = P2*x

  >>> print v_mags
  [ 1.036  1.05   1.055  1.057  1.051  1.056  1.09   1.062  1.07   1.02
    1.019  1.01   1.045  1.06 ]

  >>> print v_angs
  [-0.27995081 -0.26459191 -0.26302112 -0.2581342  -0.26354472 -0.26075219
   -0.23317599 -0.23335052 -0.24818582 -0.15323991 -0.18029251 -0.22200588
   -0.0869174   0. ]

  >>> print np.linalg.norm(x - (P1.T*v_mags+P2.T*v_angs))
  0.0

.. _net_cont:

Contingencies
=============

PFNET provides a way to specify and analyze network contingencies. A contingency is represented by an object of type :class:`Contingency <pfnet.Contingency>`, and is characterized by one or more :class:`generator <pfnet.Generator>` or :class:`branch <pfnet.Branch>` outages. The lists of generator and branch outages of a contingency can be specified at construction, or by using the class methods :func:`add_generator_outage() <pfnet.Contingency.add_generator_outage>` and :func:`add_branch_outage() <pfnet.Contingency.add_branch_outage>`, respectively. The following example shows how to construct a contingency::

  >>> import pfnet

  >>> pfnet.ParserMAT().parse('ieee14.mat')

  >>> gen = net.get_generator(3)
  >>> branch = net.get_branch(2)

  >>> c1 = pf.Contingency(generators=[gen],branches=[branch])

  >>> print c1.num_generator_outages, c1.num_branch_outages
  1 1
 
  >>> print c1.outages
  [('branch', 2), ('generator', 3)]

Once a contingency has been constructed, it can be applied and later cleared. This is done using the class methods :func:`apply() <pfnet.Contingency.apply>` and :func:`clear() <pfnet.Contingency.clear>`. The :func:`apply() <pfnet.Contingency.apply>` method sets the specified generator and branches on outage and **disconnects** them from the network. Voltage regulation and other controls provided by generators or transformers on outage are lost. The :func:`clear() <pfnet.Contingency.clear>` method undoes the changes made by the :func:`apply() <pfnet.Contingency.apply>` method. The following example shows how to apply and clear contingencies, and illustrates some of the side effects::

  >>> print c1.has_generator_outage(gen), c1.has_branch_outage(branch)
  True True

  >>> gen_bus = gen.bus
  >>> branch_bus = branch.bus_k

  >>> # generator and branch are connected to buses
  >>> print gen in gen_bus.generators, branch in branch_bus.branches
  True True

  >>> c1.apply(net)

  >>> print gen.is_on_outage(), branch.is_on_outage()
  True True

  >>> # generator and branch are disconnected from buses
  >>> print gen in gen_bus.generators, branch in branch_bus.branches
  False False

  >>> c1.clear(net)

  >>> print gen.is_on_outage(), branch.is_on_outage()
  False False

  >>> # generator and branch are connected to buses again
  >>> print gen in gen_bus.generators, branch in branch_bus.branches
  True True

More information about network contingencies can be found in the :ref:`API reference <ref_cont>`.
  
.. _net_multi_period:

Multiple Time Periods
=====================

PFNET can also be used to represent and analyze power networks over multiple time periods. By default, the networks created using most :ref:`parsers <parsers>`, as in all the examples above, have data corresponding to a single time period. To consider multiple time periods, an argument needs to be passed to the :func:`parse <pfnet.Parser>` method of a :class:`Parser <pfnet.ParserBase>`::

  >>> net = pfnet.ParserMAT().parse('ieee14.mat', num_periods=5)

  >>> print net.num_periods
  5

In "multi-period" networks, certain quantities can vary over time and hence are represented by vectors. Examples of such quantities are the :ref:`network properties <net_properties>`, generators powers, load powers, battery energy levels, bus voltage magnitudes, etc. The example below shows how to set the load profile over the time periods and extract the maximum active power mismatches in the network for each time::

  >>> import numpy as np

  >>> for load in net.loads:
  ...     load.P = np.random.rand(5)

  >>> print net.loads[0].P
  [ 0.84  0.47  0.62  0.65  0.36]

  >>> net.update_properties()

  >>> print [net.bus_P_mis[t] for t in range(5)]
  [81.92, 87.35, 86.71, 93.61, 89.90]

Lastly, for component quantities that can potentially vary over time, setting these quantities to be variables results in one variable for each time. For example, selecting the bus voltage magnitude of a bus to be variable leads to having one variable for each time period::

  >>> bus = net.buses[3]

  >>> net.set_flags_of_component(bus,'variable','voltage magnitude')

  >>> print(net.num_vars)
  5

  >>> print bus.index_v_mag
  [0 1 2 3 4]
