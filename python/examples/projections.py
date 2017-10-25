#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Power Networks - Projections

import sys
sys.path.append('.')
import pfnet
import numpy as np

net = pfnet.Parser(sys.argv[1]).parse(sys.argv[1])

net.set_flags('bus',
              'variable',
              'any',
              ['voltage magnitude','voltage angle'])

print(net.num_vars, 2*net.num_buses)

P1 = net.get_var_projection('bus', 'any', 'voltage magnitude')
P2 = net.get_var_projection('bus', 'any', 'voltage angle')

print(type(P1))

x = net.get_var_values()
v_mags = P1*x
v_angs = P2*x

print(v_mags)
print(v_angs)

print(np.linalg.norm(x - (P1.T*v_mags+P2.T*v_angs)))
