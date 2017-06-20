import numpy as np
from scipy.sparse import coo_matrix
from pfnet import CustomConstraint

class DummyConstraint(CustomConstraint):

    def init(self):
        self.name = "dummy constraint"

    def count_step(self,branch,t):
        pass

    def allocate(self):
        
        Annz = self.A_nnz
        Gnnz = self.G_nnz
        Arows = self.A_row
        Grows = self.G_row
        num_vars = self.network.num_vars

        self.set_b(np.zeros(Arows))
        self.set_A(coo_matrix((np.zeros(Annz),(Annz*[0],Annz*[0])),
                              shape=(Arows,num_vars)))

        self.set_l(np.zeros(Grows))
        self.set_u(np.zeros(Grows))
        self.set_G(coo_matrix((np.zeros(Gnnz),(Gnnz*[0],Gnnz*[0])),
                              shape=(Grows,num_vars)))

    def clear(self):
        self.A_nnz = 0
        self.G_nnz = 0
        self.A_row = 0
        self.G_row = 0
        self.bus_counted[:] = False

    def analyze_step(self,branch,t):
        pass
        
    def eval_step(self,branch,t,x):
        pass
        
    def store_sens_step(self,branch,t,sA,sf,sGu,sGl):
        pass
        
