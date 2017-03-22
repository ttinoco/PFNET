# **************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
# **************************************************#

# Power Networks - Multiperiod

import pfnet
import numpy as np

net = pfnet.ParserMAT().parse('../../data/ieee14.mat',5)

print(net.num_periods)

net.show_components()

for load in net.loads:
    load.P = np.random.rand(5)

print(net.loads[0].P)

net.update_properties()

print([net.bus_P_mis[t] for t in range(5)])

bus = net.buses[3]

net.set_flags_of_component(bus,'variable','voltage magnitude')

print(net.num_vars)

print bus.index_v_mag


