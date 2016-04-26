#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import pfnet as pf

net = pf.Network()
net.load(sys.argv[1])

print(net.num_vars)

net.set_flags(pf.OBJ_BUS,
              pf.FLAG_VARS,
              pf.BUS_PROP_REG_BY_GEN,
              pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

print(net.num_vars, 2*net.get_num_buses_reg_by_gen())

values = net.get_var_values()

print(type(values))

print(values.shape)

bus = [b for b in net.buses if b.is_regulated_by_gen()][0]

print(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))

print(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG))

print(bus.v_mag, net.get_var_values()[bus.index_v_mag])

print(bus.v_ang, net.get_var_values()[bus.index_v_ang])

print(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG))

values = net.get_var_values()

print(bus.v_mag)

values[bus.index_v_mag] = 1.20
net.set_var_values(values)

print(bus.v_mag)

net.clear_flags()

for bus in net.buses:
    if bus.index % 3 == 0:
        net.set_flags_of_component(bus,pf.FLAG_VARS,pf.BUS_VAR_VMAG)

print(net.num_vars, net.num_buses, len([b for b in net.buses if b.index % 3 == 0]))
