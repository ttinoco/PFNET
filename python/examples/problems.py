#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Optimization Problems - Problems

import pfnet
from power_flow import NRsolve

net = pfnet.ParserMAT().parse('../../data/case3012wp.mat')

print('%.2e %.2e' %(net.bus_P_mis, net.bus_Q_mis))

NRsolve(net)

print('%.2e %.2e' %(net.bus_P_mis, net.bus_Q_mis))

