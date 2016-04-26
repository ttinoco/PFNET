#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import numpy as np
from pfnet import Network

net = Network()
net.load(sys.argv[1])

print(np.average([b.degree for b in net.buses]))
