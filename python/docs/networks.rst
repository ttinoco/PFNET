.. _net:

**************
Power Networks
**************

This section describes how to load and analyze power networks using PFNET.

.. _net_overview:

Overview
========

Power networks in PFNET are represented by objects of type :class:`Network <pfnet.Network>`. These objects are initially empty and need to be loaded with data contained in specific types of files. Once the data is loaded, the network and its components can be analyzed, visualized, and used to construct network optimization problems. After a network optimization problem is solved, the network object can be updated with the solution to perform further analysis.

An important attribute of the :class:`Network <pfnet.Network>` class is :data:`base_power <pfnet.Network.base_power>`. This quantity, which has units of MVA, is useful for converting power quantities in per unit system base power to MW or MVAr.

.. _net_loading:

Loading Data
============

Power networks can be loaded with data using the :func:`load() <pfnet.Network.load>` class method. This method takes as input the filename of a supported power flow file. Information about the data parsers available in PFNET and the supported file formats can be found in Section :ref:`parsers`. The following simple example shows how to load data from a power flow ``mat`` file::

  >>> from pfnet import Network

  >>> net = Network()
  >>> print net.num_buses
  0

  >>> net.load('ieee14.mat')
  >>> print net.num_buses
  14

.. _net_components:

Components
==========

Power networks have several components. These are :ref:`buses <net_bus>`, :ref:`branches <net_branch>`, :ref:`generators <net_gen>`, :ref:`shunt devices <net_shunt>`, :ref:`loads <net_load>`, and :ref:`variable generators <net_vargen>` (*i.e.*, non-dispatchable). For obtaining an overview of the components that form a network, the class method :func:`show_components() <pfnet.Network.show_components>` can be used::

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
  loads            : 11
  vargens          : 0

.. _net_bus:

Buses
-----

Buses in a power network are objects of type :class:`Bus <pfnet.Bus>`. Each bus has an :data:`index <pfnet.Bus.index>`, a :data:`number <pfnet.Bus.number>`, and a :data:`name <pfnet.Bus.name>` attribute that can be used to identify this bus in a network. The :data:`index <pfnet.Bus.index>` is associated with the location of the bus in the underlying C array of bus structures, while the :data:`number <pfnet.Bus.number>` and :data:`name <pfnet.Bus.name>` attributes are specified in the input data. An :data:`index <pfnet.Bus.index>`, a :data:`number <pfnet.Bus.number>`, or a :data:`name <pfnet.Bus.name>` can be used to extract a specific bus from a network using the :class:`Network <pfnet.Network>` class methods :func:`get_bus() <pfnet.Network.get_bus>`, :func:`get_bus_by_number() <pfnet.Network.get_bus_by_number>`, and :func:`get_bus_by_name() <pfnet.Network.get_bus_by_name>`, respectively::
  
  >>> bus = net.get_bus(10)

  >>> print bus.index == 10
  True

  >>> other_bus = net.get_bus_by_number(bus.number)

  >>> print bus == other_bus
  True

For convenience, a list of all the buses in the network is contained in the :data:`buses <pfnet.Network.buses>` attribute of the :class:`Network <pfnet.Network>` class.

Buses in a network can have different properties. For example, some buses can be slack buses and others can have their voltage magnitudes regulated by generators, tap-changing transformers, or switched shunt devices. The :class:`Bus <pfnet.Bus>` class provides methods for checking whether a bus has specific properties. The following example shows how to get a list of all the buses whose voltage magnitudes are regulated by generators::

  >>> reg_buses = [b for b in net.buses if b.is_regulated_by_gen()]

  >>> print len(reg_buses), net.get_num_buses_reg_by_gen()
  5 5

A bus also has information about the devices that are connected to it or that are regulating its voltage magnitude. For example, the attributes :data:`gens <pfnet.Bus.gens>` and :data:`reg_trans <pfnet.Bus.reg_trans>` contain a list of generators connected to the bus and a list of tap-changing transformers regulating its voltage magnitude, respectively.

.. _net_branch:

Branches
--------

Branches in a power network are objects of type :class:`Branch <pfnet.Branch>` and are represented mathematically by the model described in Section 2.1.2 of [TTR2015]_. Each branch has an :data:`index <pfnet.Branch.index>` attribute that can be used to identify this branch in a network. The :class:`Network <pfnet.Network>` class method :func:`get_branch() <pfnet.Network.get_branch>` can be used to extract a branch of a given :data:`index <pfnet.Branch.index>`::
  
  >>> branch = net.get_branch(5)

  >>> print branch.index == 5
  True

For convenience, a list of all the branches in the network is contained in the :data:`branches <pfnet.Network.branches>` attribute of the :class:`Network <pfnet.Network>` class.

Branches in a power network can have different properties. Fore example, some branches can be transmission lines, fixed transformers, tap-changing transformers, or phase-shifting transformers. Tap-changing transformers in turn can control the reactive power flowing through the branch or the voltage magnitude of a bus. The :class:`Branch <pfnet.Branch>` class provides methods for checking whether a branch has specific properties. The following example shows how to get a list of all the branches that are transmission lines::

  >>> lines = [br for br in net.branches if br.is_line()]

  >>> print len(lines), net.get_num_lines()
  17 17

For branches that are transformers, the :class:`Branch <pfnet.Branch>` class attributes :data:`ratio <pfnet.Branch.ratio>` and :data:`phase <pfnet.Branch.phase>` correspond to the transformer's tap ratio and phase shift, respectively. These attributes correspond to the quantities :math:`a_{km}` and :math:`\phi_{km}` of the branch model described in Section 2.1.2 of [TTR2015]_. The quantity :math:`a_{mk}` in this model is always one.

.. _net_gen:

Generators
----------

Generators in a power network are objects of type :class:`Generator <pfnet.Generator>`. Each generator has an :data:`index <pfnet.Generator.index>` attribute that can be used to identify this generator in a network. The :class:`Network <pfnet.Network>` class method :func:`get_gen() <pfnet.Network.get_gen>` can be used to extract a generator of a given :data:`index <pfnet.Generator.index>`::
  
  >>> gen = net.get_gen(2)

  >>> print gen.index == 2
  True

For convenience, a list of all the generators in the network is contained in the :data:`generators <pfnet.Network.generators>` attribute of the :class:`Network <pfnet.Network>` class.

Generators in a power network can have different properties. Fore example, some generators can be slack generators and others can provide bus voltage magnitude regulation. The :class:`Generator <pfnet.Generator>` class provides methods for checking whether a generator has specific properties. The following example shows how to get a list of all the slack generators::

  >>> slack_gens = [g for g in net.generators if g.is_slack()]

  >>> print len(slack_gens), net.get_num_slack_gens()
  1 1

The active and reactive powers that a generator injects into the bus to which it is connected are obtained from the :data:`P <pfnet.Generator.P>` and :data:`Q <pfnet.Generator.Q>` attributes of the :class:`Generator <pfnet.Generator>` class. These quantities are given in units of per unit :data:`system base power <pfnet.Network.base_power>`. The following example computes the total active power injected into the network by generators in units of MW::

  >>> print sum([g.P for g in net.generators])*net.base_power
  272.4
  
.. _net_shunt:

Shunt Devices
-------------

Shunt devices in a power network are objects of type :class:`Shunt <pfnet.Shunt>`. Each shunt has an :data:`index <pfnet.Shunt.index>` attribute that can be used to identify this shunt in a network. The :class:`Network <pfnet.Network>` class method :func:`get_shunt() <pfnet.Network.get_shunt>` can be used to extract a shunt of a given :data:`index <pfnet.Shunt.index>`::
  
  >>> shunt = net.get_shunt(0)

  >>> print shunt.index == 0
  True

For convenience, a list of all the shunt devices in the network is contained in the :data:`shunts <pfnet.Network.shunts>` attribute of the :class:`Network <pfnet.Network>` class.

As other network components, shunt devices can have different properties. Some shunt devices can be fixed while others can be switchable and configured to regulate a bus voltage magnitude.

.. _net_load:

Loads
-----

Loads in a power network are objects of type :class:`Load <pfnet.Load>`. As other components, the :data:`index <pfnet.Load.index>` attribute is used to identify a load in the network. A list of all the loads in the network is contained in the :data:`loads <pfnet.Network.loads>` attribute of the :class:`Network <pfnet.Network>` class. 

Similar to generators, the active and reactive powers that a load consumes from the bus to which it is connected are obtained from the :data:`P <pfnet.Load.P>` and :data:`Q <pfnet.Load.Q>` attributes of the :class:`Load <pfnet.Load>` class. They are also given in units of per unit :data:`system base power <pfnet.Network.base_power>`.

.. _net_vargen:

Variable Generators
-------------------

Variable generators in a power network are objects of type :class:`VarGenerator <pfnet.VarGenerator>`. They represent non-dispatchable energy sources such as wind generators or farms and photovoltaic power plants. As with other components, the :data:`index <pfnet.VarGenerator.index>` attribute is used to identify a variable generator in the network. In addition to the :data:`index <pfnet.VarGenerator.index>` attribute, a :data:`name <pfnet.VarGenerator.name>` attribute is also available, which can be used to extract a specific variable generator from the network using the :class:`Network <pfnet.Network>` class method :func:`get_vargen_by_name() <pfnet.Network.get_vargen_by_name>`. A list of all the variable generators in the network is also contained in the :data:`var_generators <pfnet.Network.var_generators>` attribute of the :class:`Network <pfnet.Network>` class. 

Similar to generators, the active and reactive powers produced by a variable generator are obtained from the :data:`P <pfnet.VarGenerator.P>` and :data:`Q <pfnet.VarGenerator.Q>` attributes of the :class:`VarGenerator <pfnet.VarGenerator>` class in units of per unit :data:`system base power <pfnet.Network.base_power>`. This is the output of the device in the absence of uncertainty. When there is uncertainty, the output of the device is subject to variations about :data:`P <pfnet.VarGenerator.P>` that have a standard deviation given by the attribute :data:`P_std <pfnet.VarGenerator.P_std>`. Output limits of a variable generator are given by the :data:`P_min <pfnet.VarGenerator.P_min>`, :data:`P_max <pfnet.VarGenerator.P_max>`, :data:`Q_min <pfnet.VarGenerator.Q_min>`, and :data:`Q_max <pfnet.VarGenerator.Q_max>` attributes. 

The output of variable generators in a network are subject to random variations that can be correlated, especially for devices that are "nearby". The method :func:`create_vargen_P_sigma() <pfnet.Network.create_vargen_P_sigma>` of the :class:`Network <pfnet.Network>` class allows constructing a covariance matrix for these variations based on a "correlation distance" ``N`` and a given correlation coefficient. The cross-covariance between the variation of two devices that are connected to buses that are less than ``N`` branches away from each other are set such that they have the given correlation coefficient.

Lastly, since many power network input files do not have variable generator information, these devices can be added to the network by using the :func:`add_vargens() <pfnet.Network.add_vargens>` method of the :class:`Network <pfnet.Network>` class.

.. _net_properties:

Properties
==========

A :class:`Network <pfnet.Network>` object has several quantities or ``properties`` that provide important information about the state of the network. The following table provides a description of each of these properties.

=============== ================================================================= ========
Names           Description                                                       Units
=============== ================================================================= ========
``bus_v_max``   Maximum bus voltage magnitude                                     per unit
``bus_v_min``   Minimum bus voltage magnitude                                     per unit
``bus_v_vio``   Maximum bus voltage magnitude limit violation                     per unit
``bus_P_mis``   Maximum absolute bus active power mismatch                        MW
``bus_Q_mis``   Maximum absolute bus reactive power mismatch                      MVAr
``gen_v_dev``   Maximum set point deviation of generator-regulated voltage        per unit
``gen_Q_vio``   Maximum generator reactive power limit violation                  MVAr
``gen_P_vio``   Maximum generator active power limit violation                    MW
``tran_v_vio``  Maximum band violation of transformer-regulated voltage           per unit
``tran_r_vio``  Maximum tap ratio limit violation of tap-changing transformer     unitless
``tran_p_vio``  Maximum phase shift limit violation of phase-shifting transformer radians
``shunt_v_vio`` Maximum band violation of shunt-regulated voltage                 per unit
``shunt_b_vio`` Maximum susceptance limit violation of switched shunt device      per unit
``num_actions`` Number of control adjustments (greater than 2% of control range)  unitless
=============== ================================================================= ========

All of these properties are attributes of the :class:`Network <pfnet.Network>` class. If there is a change in the network, the class method :func:`update_properties() <pfnet.Network.update_properties>` needs to be called in order for the network properties to reflect the change. The following example shows how to update and extract properties::

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

Network quantities can be specified to be ``variables``. This is useful to represent network quantities with vectors and turn the network properties described above as functions of these vectors. 

To set network quantities as variables, the :class:`Network <pfnet.Network>` class method :func:`set_flags() <pfnet.Network.set_flags>` is used. This method takes as arguments a :ref:`component type <ref_net_obj>`, a :ref:`flag mask <ref_net_flag>` for specifying which flags types to set, a ``property mask`` for targeting objects with specific properties, and a ``variable mask`` for specifying which component quantities should be affected.

**Property masks** are component-specific. They can be combined using ``logical OR`` to make properties more complex. More information can be found in the following sections:

* :ref:`ref_bus_prop`
* :ref:`ref_branch_prop`
* :ref:`ref_gen_prop`
* :ref:`ref_shunt_prop`
* :ref:`ref_vargen_prop`

**Variable masks** are also component-specific. They can be combined using ``logical OR`` to target more than one component quantity. More information can be found in the following sections:

* :ref:`ref_bus_var`
* :ref:`ref_branch_var`
* :ref:`ref_gen_var`
* :ref:`ref_shunt_var`
* :ref:`ref_vargen_var`

The following example shows how to set as variables all the voltage magnitudes and angles of buses regulated by generators::

  >>> import pfnet as pf

  >>> net = pf.Network()
  >>> net.load('ieee14.mat')

  >>> print net.num_vars
  0

  >>> net.set_flags(pf.OBJ_BUS,
  ...               pf.FLAG_VARS,
  ...               pf.BUS_PROP_REG_BY_GEN,
  ...               pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

  >>> print net.num_vars, 2*net.get_num_buses_reg_by_gen()
  10 10

Network components have a :func:`has_flags() <pfnet.Bus.has_flags>` method that allows checking whether flags of a certain type associated with specific quantities are set.

Once variables have been set, the :ref:`vector <ref_vec>` containing all the current variable values can be extracted using :func:`get_var_values() <pfnet.Network.get_var_values>`::

  >>> values = net.get_var_values()
  
  >>> print type(values)
  <type 'numpy.ndarray'>

  >>> print values.shape
  (10,)

The components that have quantities set as variables have indices that can be used to locate these quantities in the vector of all variable values::

  >>> bus = [b for b in net.buses if b.is_reg_by_gen()][0]

  >>> print bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG)
  True

  >>> bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG)
  True

  >>> print bus.v_mag, net.get_var_values()[bus.index_v_mag]
  1.09 1.09

  >>> print bus.v_ang, net.get_var_values()[bus.index_v_ang]
  -0.23 -0.23

A vector of variable values can be used to update the corresponding network quantities. This is done with the :class:`Network <pfnet.Network>` class method :func:`set_var_values() <pfnet.Network.set_var_values>`::

  >>> bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG)
  True

  >>> values = net.get_var_values()

  >>> print bus.v_mag
  1.09

  >>> values[bus.index_v_mag] = 1.20
  >>> net.set_var_values(values)

  >>> print bus.v_mag
  1.20

As we will see in later, variables are also useful for constructing network optimization problems.

Lastly, the class method :func:`get_var_values() <pfnet.Network.get_var_values>` can also be used to get upper or lower limits of the variables. To do this, a valid :ref:`variable value code <ref_var_values>` must be passed to this method.

.. _net_var_projections:

Projections
===========

As explained above, once the network variables have been set, a vector with the current values of the selected variables is obtained with the class method :func:`get_var_values() <pfnet.Network.get_var_values>`. To extract subvectors that contain values of specific variables, projection matrices can be used. These :ref:`matrices <ref_mat>` can be obtained using the class method :func:`get_var_projection() <pfnet.Network.get_var_projection>`, which take as arguments a :ref:`component type <ref_net_obj>` and a ``variable mask``, *e.g.*, :ref:`bus variable masks <ref_bus_var>`. The next example sets the variables of the network to be the bus voltage magnitudes and angles of all the buses, extracts the vector of values of all variables, and then extracts two subvectors having only voltage magnitudes and only voltage angles, respectively::

  >>> import numpy as np
  >>> import pfnet as pf

  >>> net = pf.Network()
  >>> net.load('ieee14.mat')

  >>> net.set_flags(pf.OBJ_BUS,
  ...               pf.FLAG_VARS,
  ...               pf.BUS_PROP_ANY,
  ...               pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

  >>> print net.num_vars, 2*net.num_buses
  28 28

  >>> P1 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VMAG)
  >>> P2 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)

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


