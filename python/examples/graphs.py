#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Visualization - Overview

import pfnet

net = pfnet.ParserMAT().parse('../../data/ieee14.mat')

g = pfnet.Graph(net)

for bus in net.buses:
    g.set_node_property(bus,"label",str(bus.number))

g.color_nodes_by_mismatch(pfnet.BUS_MIS_REACTIVE)

g.set_layout()

g.write('png','graph.png')

g.view()
