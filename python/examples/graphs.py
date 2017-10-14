#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Visualization - Overview

import sys
sys.path.append('.')
import pfnet

net = pfnet.Parser(sys.argv[1]).parse(sys.argv[1])

g = pfnet.Graph(net)

for bus in net.buses:
    g.set_node_property(bus,"label",str(bus.number))

g.color_nodes_by_mismatch(pfnet.BUS_MIS_REACTIVE)

g.set_layout()

g.write('png','graph.png')

g.view()
