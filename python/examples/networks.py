#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Power Networks - Components

import sys
sys.path.append('.')
import pfnet

net = pfnet.Parser(sys.argv[1]).parse(sys.argv[1])
net.show_components()

# Power Networks - Buses 

bus = net.get_bus(10)

print(bus.index == 10)

other_bus = net.get_bus_from_number(bus.number)

print(bus == other_bus)

reg_buses = [bus for bus in net.buses if bus.is_regulated_by_gen()]

print(len(reg_buses), net.get_num_buses_reg_by_gen())

# Power Networks - Branches

branch = net.get_branch(5)

print(branch.index == 5)

lines = [br for br in net.branches if br.is_line()]

print(len(lines), net.get_num_lines())

# Power Networks - Generators

gen = net.get_generator(2)

print(gen.index == 2)

slack_gens = [g for g in net.generators if g.is_slack()]

print(len(slack_gens), net.get_num_slack_gens())

print(sum([g.P for g in net.generators]) * net.base_power)

# Power Networks - Shunt Devices

shunt = net.get_shunt(0)

print(shunt.index == 0)

# Power Networks - Variable Generators

buses = net.get_generator_buses()
net.add_var_generators_from_parameters(buses, 80., 50., 50., 5, 0.05)

print(net.num_var_generators == len([bus for bus in net.buses if bus.generators]))

vargen = net.get_var_generator(3)

print(vargen.index == 3)

print(vargen.P, vargen.P_std, vargen.P_max)

# Power Networks - Properties

print(net.bus_v_max)

for bus in net.buses:
    bus.v_mag = bus.v_mag + 0.1

print(net.bus_v_max)

net.update_properties()

print(net.bus_v_max)

properties = net.get_properties()

print(properties['bus_v_max'])
