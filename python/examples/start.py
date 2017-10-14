#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Getting started - Installation

import pfnet
print(pfnet.info)

# Getting Started - Example

import sys
import pfnet
import numpy as np

net = pfnet.Parser(sys.argv[1]).parse(sys.argv[1])

print(np.average([bus.degree for bus in net.buses]))
