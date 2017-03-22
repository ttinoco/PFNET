#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Power Networks - Variables

import pfnet

net = pfnet.ParserMAT().parse('../../data/ieee14.mat')

print(net.num_vars)

net.set_flags('bus',
              'variable',
              'regulated by generator',
              ['voltage magnitude', 'voltage angle'])

print(net.num_vars, 2*net.get_num_buses_reg_by_gen())

values = net.get_var_values()

print(type(values))

print(values.shape)

bus = [bus for bus in net.buses if bus.is_regulated_by_gen()][0]

print(bus.has_flags('variable','voltage magnitude'))

print(bus.has_flags('variable','voltage angle'))

print(bus.v_mag, net.get_var_values()[bus.index_v_mag])

print(bus.v_ang, net.get_var_values()[bus.index_v_ang])

print(bus.has_flags('variable','voltage angle'))

values = net.get_var_values()

print(bus.v_mag)

values[bus.index_v_mag] = 1.20
net.set_var_values(values)

print(bus.v_mag)

net.clear_flags()

for bus in net.buses:
    if bus.index % 3 == 0:
        net.set_flags_of_component(bus,'variable','voltage magnitude')

print(net.num_vars, len([bus for bus in net.buses if bus.index % 3 == 0]), net.num_buses)
