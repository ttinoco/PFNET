#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2016, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import numpy as np
import pfnet as pf

net = pf.Network()
net.load(sys.argv[1])

net.set_flags(pf.OBJ_BUS,
              pf.FLAG_VARS,
              pf.BUS_PROP_ANY,
              pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

print net.num_vars, 2*net.num_buses

P1 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VMAG)
P2 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)

print type(P1)

x = net.get_var_values()
v_mags = P1*x
v_angs = P2*x

print v_mags
print v_angs

print np.linalg.norm(x - (P1.T*v_mags+P2.T*v_angs))
