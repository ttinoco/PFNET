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

net.set_flags(pf.OBJ_BUS,
              pf.FLAG_VARS,
              pf.BUS_PROP_ANY,
              pf.BUS_VAR_VMAG)

print(net.num_vars == net.num_buses)

func = pf.Function(pf.FUNC_TYPE_REG_VMAG,0.3,net)

print(func.type == pf.FUNC_TYPE_REG_VMAG)

print(func.weight)

x = net.get_var_values()

func.analyze()
func.eval(x)

print(x.shape)

print(func.phi)

print(type(func.gphi), func.gphi.shape)

print(type(func.Hphi), func.Hphi.shape)

print(func.Hphi)


