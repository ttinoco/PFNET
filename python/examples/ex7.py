#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
from pfnet import Network
from ex6 import NRsolve

net = Network()
net.load(sys.argv[1])

print '%.2e %.2e' %(net.bus_P_mis, net.bus_Q_mis)

NRsolve(net)

print '%.2e %.2e' %(net.bus_P_mis, net.bus_Q_mis)
