#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Getting started - Installation

import pfnet
print(pfnet.info)

# Getting Started - Example

import pfnet
import numpy as np

net = pfnet.ParserMAT().parse('../../data/ieee14.mat')

print(np.average([bus.degree for bus in net.buses]))
