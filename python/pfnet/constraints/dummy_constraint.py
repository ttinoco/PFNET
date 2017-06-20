#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import numpy as np
from scipy.sparse import coo_matrix
from pfnet import CustomConstraint

class DummyDCPF(CustomConstraint):

    def init(self):

        self.name = "dummy DC power balance"

    def count_step(self,branch,t):

        buses = [branch.bus_k,branch.bus_m]
        
        if branch.is_on_outage():
            return

        for k in range(2):
            m = 1 if k == 0 else 0
            if buses[k].has_flags('variable','voltage angle'):
                self.A_nnz = self.A_nnz + 1
            if buses[m].has_flags('variable','voltage angle'):
                self.A_nnz = self.A_nnz + 1
            if branch.has_flags('variable','phase shift'):
                self.A_nnz = self.A_nnz + 1
                
        for bus in buses:
            index = bus.index+t*self.network.num_buses
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    if gen.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for load in bus.loads:
                    if load.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for vargen in bus.var_generators:
                    if vargen.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for bat in bus.batteries:
                    if bat.has_flags('variable','charging power'):
                        self.A_nnz = self.A_nnz + 2
            self.bus_counted[index] = True

    def allocate(self):
        
        nnz = self.A_nnz
        num_constr = self.network.num_buses*self.network.num_periods

        self.set_b(np.zeros(num_constr))
        self.set_A(coo_matrix((np.zeros(nnz),(nnz*[0],nnz*[0])),
                              shape=(num_constr,self.network.num_vars)))

        self.set_f(np.zeros(0))
        self.set_J(coo_matrix((0,self.network.num_vars)))

        self.set_l(np.zeros(0))
        self.set_u(np.zeros(0))
        self.set_G(coo_matrix((0,self.network.num_vars)))

    def clear(self):

        self.A_nnz = 0
        self.bus_counted[:] = False

    def analyze_step(self,branch,t):
        
        buses = [branch.bus_k,branch.bus_m]
        
        if branch.is_on_outage():
            return

        for k in range(2):
            m = 1 if k == 0 else 0
            sign_phi = 1. if k == 0 else -1.
            index = buses[k].index+t*self.network.num_buses
            if buses[k].has_flags('variable','voltage angle'):
                self.A.row[self.A_nnz] = index
                self.A.col[self.A_nnz] = buses[k].index_v_ang[t]
                self.A.data[self.A_nnz] = branch.b
                self.A_nnz = self.A_nnz + 1
            else:
                self.b[index] += -branch.b*buses[k].v_ang[t]
            if buses[m].has_flags('variable','voltage angle'):
                self.A.row[self.A_nnz] = index
                self.A.col[self.A_nnz] = buses[m].index_v_ang[t]
                self.A.data[self.A_nnz] = -branch.b
                self.A_nnz = self.A_nnz + 1
            else:
                self.b[index] += branch.b*buses[m].v_ang[t]
            if branch.has_flags('variable','phase shift'):
                self.A.row[self.A_nnz] = index
                self.A.col[self.A_nnz] = branch.index_phase[t]
                self.A.data[self.A_nnz] = -branch.b*sign_phi
                self.A_nnz = self.A_nnz + 1
            else:
                self.b[index] += branch.b*branch.phase[t]*sign_phi
            
        for bus in buses:
            index = bus.index+t*self.network.num_buses
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    if gen.has_flags('variable','active power'):
                        self.A.row[self.A_nnz] = index
                        self.A.col[self.A_nnz] = gen.index_P[t]
                        self.A.data[self.A_nnz] = 1.
                        self.A_nnz = self.A_nnz + 1
                    else:
                        self.b[index] += -gen.P[t]
                for load in bus.loads:
                    if load.has_flags('variable','active power'):
                        self.A.row[self.A_nnz] = index
                        self.A.col[self.A_nnz] = load.index_P[t]
                        self.A.data[self.A_nnz] = -1.
                        self.A_nnz = self.A_nnz + 1
                    else:
                        self.b[index] += load.P[t]
                for vargen in bus.var_generators:
                    if vargen.has_flags('variable','active power'):
                        self.A.row[self.A_nnz] = index
                        self.A.col[self.A_nnz] = vargen.index_P[t]
                        self.A.data[self.A_nnz] = 1.
                        self.A_nnz = self.A_nnz + 1
                    else:
                        self.b[index] += -vargen.P[t]
                for bat in bus.batteries:
                    if bat.has_flags('variable','charging power'):
                        self.A.row[self.A_nnz] = index
                        self.A.col[self.A_nnz] = bat.index_Pc[t]
                        self.A.data[self.A_nnz] = -1.
                        self.A_nnz = self.A_nnz + 1
                        self.A.row[self.A_nnz] = index
                        self.A.col[self.A_nnz] = bat.index_Pd[t]
                        self.A.data[self.A_nnz] = 1.
                        self.A_nnz = self.A_nnz + 1
                    else:
                        self.b[index] += bat.P[t]
            self.bus_counted[index] = True

    def eval_step(self,branch,t,x,y=None):
 
        pass
        
    def store_sens_step(self,branch,t,sA,sf,sGu,sGl):
        
        pass
        
