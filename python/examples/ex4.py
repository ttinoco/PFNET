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

net.set_flags('bus',
              'variable',
              'any',
              'voltage magnitude')

print(net.num_vars == net.num_buses)

func = pf.Function('voltage magnitude regularization',0.3,net)

print(func.type == 'voltage magnitude regularization')

print(func.weight)

x = net.get_var_values()

func.analyze()
func.eval(x)

print(x.shape)

print(func.phi)

print(type(func.gphi), func.gphi.shape)

print(type(func.Hphi), func.Hphi.shape)

print(func.Hphi)


