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
            index = bus.index*t+self.network.num_buses
            if not self.bus_counted[index]:
                for gen in bus.generators:
                    if gen.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for load in bus.load:
                    if load.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for vargen in bus.var_generators:
                    if vargen.has_flags('variable','active power'):
                        self.A_nnz = self.A_nnz + 1
                for bat in bus.batteries:
                    if bat.has_flags('variable','charging power'):
                        self.A_nnz = self.A_nnz + 1
            self.bus_counted[index] = True

    def allocate(self):

        pass

    def clear(self):

        pass

    def analyze_step(self,branch,t):
        
        pass

    def eval_step(self,branch,t,x):
 
        pass

    def store_sens_step(self,branch,t,sA,sf,sGu,sGl):
        
        pass
