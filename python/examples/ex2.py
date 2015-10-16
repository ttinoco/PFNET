#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
from pfnet import Network

net = Network()
print net.num_buses

net.load(sys.argv[1])
print net.num_buses

net.show_components()

bus = net.get_bus(10)

print bus.index == 10

other_bus = net.get_bus_by_number(bus.number)

print bus == other_bus

reg_buses = [b for b in net.buses if b.is_regulated_by_gen()]

print len(reg_buses), net.get_num_buses_reg_by_gen()

branch = net.get_branch(5)

print branch.index == 5

lines = [br for br in net.branches if br.is_line()]

print len(lines), net.get_num_lines()

gen = net.get_gen(2)

print gen.index == 2

slack_gens = [g for g in net.generators if g.is_slack()]

print len(slack_gens), net.get_num_slack_gens()

print sum([g.P for g in net.generators])*net.base_power

shunt = net.get_shunt(0)

print shunt.index == 0

net.add_vargens(net.get_gen_buses(),50.,50.,5,0.05)

print net.num_vargens == len([b for b in net.buses if b.gens])

vargen = net.get_vargen(3)

print vargen.index == 3

print vargen.P,vargen.P_std,vargen.P_max

print net.bus_v_max

for bus in net.buses:

    bus.v_mag = bus.v_mag + 0.1

print net.bus_v_max

net.update_properties()

print net.bus_v_max

properties = net.get_properties()

print properties['bus_v_max']
