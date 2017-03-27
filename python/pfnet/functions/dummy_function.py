#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import numpy as np
from scipy.sparse import coo_matrix
from pfnet import CustomFunction

class DummyGenCost(CustomFunction):

    def init(self):
        
        self.name = "dummy generation cost"
    
    def count_step(self,branch,t):

        T = self.network.num_periods
        buses = [branch.bus_k,branch.bus_m]
        
        if branch.is_on_outage():
            return

        for bus in buses:
            index = bus.index*T+t
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    if gen.has_flags('variable','active power'):
                        self.Hphi_nnz = self.Hphi_nnz+1
            self.bus_counted[index] = True
        
    def allocate(self):

        nnz = self.Hphi_nnz
        num_vars = self.network.num_vars
        self.set_gphi(np.zeros(self.network.num_vars))
        self.set_Hphi(coo_matrix((np.zeros(nnz),(nnz*[0],nnz*[0])),
                                 shape=(num_vars,num_vars)))

    def clear(self):
        
        self.phi = 0
        self.gphi[:] = 0
        self.Hphi_nnz = 0
        self.bus_counted[:] = False
                
    def analyze_step(self,branch,t):
        
        T = self.network.num_periods
        buses = [branch.bus_k,branch.bus_m]
        
        if branch.is_on_outage():
            return

        for bus in buses:
            index = bus.index*T+t
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    if gen.has_flags('variable','active power'):
                        self.Hphi.row[self.Hphi_nnz] = gen.index_P[t]
                        self.Hphi.col[self.Hphi_nnz] = gen.index_P[t]
                        self.Hphi.data[self.Hphi_nnz] = 2.*gen.cost_coeff_Q2
                        self.Hphi_nnz = self.Hphi_nnz+1
            self.bus_counted[index] = True

    def eval_step(self,branch,t,x):

        T = self.network.num_periods
        buses = [branch.bus_k,branch.bus_m]
        
        if branch.is_on_outage():
            return

        for bus in buses:
            index = bus.index*T+t
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    Q0 = gen.cost_coeff_Q0
                    Q1 = gen.cost_coeff_Q1
                    Q2 = gen.cost_coeff_Q2
                    if gen.has_flags('variable','active power'):
                        P = x[gen.index_P[t]]
                        self.gphi[gen.index_P[t]] = Q1+2.*Q2*P
                    else:
                        P = gen.P[t]
                    self.phi = self.phi + Q0+Q1*P+Q2*(P**2.)
            self.bus_counted[index] = True
