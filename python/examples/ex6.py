#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
from numpy import hstack
from numpy.linalg import norm
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve

def NRsolve(net):

    net.clear_flags()

    # bus voltage angles
    net.set_flags(pf.OBJ_BUS,
                  pf.FLAG_VARS,
                  pf.BUS_PROP_NOT_SLACK,
                  pf.BUS_VAR_VANG)
    
    # bus voltage magnitudes
    net.set_flags(pf.OBJ_BUS,
                  pf.FLAG_VARS,
                  pf.BUS_PROP_NOT_REG_BY_GEN,
                  pf.BUS_VAR_VMAG)
    
    # slack gens active powers
    net.set_flags(pf.OBJ_GEN,
                  pf.FLAG_VARS,
                  pf.GEN_PROP_SLACK,
                  pf.GEN_VAR_P)
    
    # regulator gens reactive powers
    net.set_flags(pf.OBJ_GEN,
                  pf.FLAG_VARS,
                  pf.GEN_PROP_REG,
                  pf.GEN_VAR_Q)

    p = pf.Problem()
    p.set_network(net)
    p.add_constraint(pf.CONSTR_TYPE_PF)       # power flow
    p.add_constraint(pf.CONSTR_TYPE_PAR_GEN)  # generator participation
    p.analyze()
    
    x = p.get_init_point()
    p.eval(x)

    residual = lambda x: hstack((p.A*x-p.b,p.f))

    while norm(residual(x)) > 1e-4:
        x = x + spsolve(bmat([[p.A],[p.J]],format='csr'),-residual(x))
        p.eval(x)

    net.set_var_values(x)
    net.update_properties()

