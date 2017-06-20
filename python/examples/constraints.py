#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Optimization Problems - Constraints

import pfnet

net = pfnet.ParserMAT().parse('../../data/ieee14.mat')

net.set_flags('bus',
              'variable',
              'any',
              ['voltage magnitude','voltage angle'])

print(net.num_vars == 2*net.num_buses)

constr = pfnet.Constraint('AC power balance',net)

print(constr.name == 'AC power balance')

x = net.get_var_values()

constr.analyze()

print(constr.num_extra_vars)

constr.eval(x + 0.01)
constr.eval(x)

import numpy as np

f = constr.f

print(type(f), f.shape)

print(np.linalg.norm(f,np.inf))

bus = net.get_bus(5)

Hi = constr.get_H_single(bus.index_P)

print(type(Hi), Hi.shape, Hi.nnz)

coefficients = np.random.randn(f.size)

constr.combine_H(coefficients)

H = constr.H_combined

print(type(H), H.shape, H.nnz)

  


