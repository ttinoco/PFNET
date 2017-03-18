import pfnet
import numpy as np
from scipy.sparse import coo_matrix

net = pfnet.Network()
net.load('../data/aesoSL2014.raw')

f = pfnet.functions.DummyGenCost(1.,net)

net.set_flags('generator',
              'variable',
              'any',
              'active power')

x = net.get_var_values()

f.analyze()
