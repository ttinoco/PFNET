import numpy as np
from scipy.sparse import coo_matrix
from pfnet import CustomFunction

class DummyFunction(CustomFunction):

    def init(self):        
        self.name = "dummy function"
    
    def count_step(self,branch,t):
        pass

    def clear(self):
        self.phi = 0
        self.gphi[:] = 0
        self.Hphi_nnz = 0
        self.bus_counted[:] = False
        
    def allocate(self):
        nnz = self.Hphi_nnz
        num_vars = self.network.num_vars
        self.set_gphi(np.zeros(self.network.num_vars))
        self.set_Hphi(coo_matrix((np.zeros(nnz),(nnz*[0],nnz*[0])),
                                 shape=(num_vars,num_vars)))
                
    def analyze_step(self,branch,t):
        pass

    def eval_step(self,branch,t,x):
        pass
