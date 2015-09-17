#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import pfnet as pf

net = pf.Network()

net.load(sys.argv[1])

g = pf.Graph(net)

g.set_layout()

#g.color_nodes_by_mismatch(pf.BUS_MIS_REACTIVE)

#g.set_edges_property("color","black")
#g.set_nodes_property("color","black")

g.write('pdf','graph.pdf')

g.view()
