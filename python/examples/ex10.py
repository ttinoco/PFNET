#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import numpy as np
import pfnet as pf

net = pf.Network()
net.load(sys.argv[1])

gen = net.get_gen(3)
gen_bus = gen.bus
branch = net.get_branch(2)
branch_bus = branch.bus_from

c1 = pf.Contingency(gens=[gen],branches=[branch])

print c1.num_gen_outages, c1.num_branch_outages

print gen in gen_bus.gens, branch in branch_bus.branches

c1.apply()

print gen.is_on_outage(), branch.is_on_outage()

print gen in gen_bus.gens, branch in branch_bus.branches

c1.clear()

print gen.is_on_outage(), branch.is_on_outage()

print gen in gen_bus.gens, branch in branch_bus.branches
