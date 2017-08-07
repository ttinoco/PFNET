import numpy as np
from scipy.sparse import coo_matrix
from pfnet import CustomConstraint

class DummyConstraint(CustomConstraint):

    def init(self):
        self.name = "dummy constraint"
        self.set_H_nnz(np.zeros(large_enough_number,dtype='int32'))

    def count_step(self,branch,t):
        pass

    def allocate(self):
        
        Annz = self.A_nnz
        Gnnz = self.G_nnz
        Jnnz = self.J_nnz
        Arows = self.A_row
        Grows = self.G_row
        Jrows = self.J_rows
        Hnnz = self.H_nnz
        num_vars = self.network.num_vars

        self.set_b(np.zeros(Arows))
        self.set_A(coo_matrix((np.zeros(Annz),(Annz*[0],Annz*[0])),
                              shape=(Arows,num_vars)))

        self.set_l(np.zeros(Grows))
        self.set_u(np.zeros(Grows))
        self.set_G(coo_matrix((np.zeros(Gnnz),(Gnnz*[0],Gnnz*[0])),
                              shape=(Grows,num_vars)))

        self.set_f(np.zeros(Jrows))
        self.set_J(coo_matrix((np.zeros(Jnnz),(Jnnz*[0],Jnnz*[0])),
                              shape=(Jrows,num_vars)))
        self.allocate_H_array(Jrows)
        for i in range(Jrows):
            self.set_H_single(i,
                              coo_matrix((np.zeros(Hnnz[i]),(Hnnz[i]*[0],Hnnz[i]*[0])),
                                         shape=(num_vars,num_vars)))

    def clear(self):
        self.A_nnz = 0
        self.G_nnz = 0
        self.J_nnz = 0
        self.A_row = 0
        self.G_row = 0
        self.J_row = 0
        self.H_nnz[:] = 0
        self.bus_counted[:] = False

    def analyze_step(self,branch,t):
        pass
        
    def eval_step(self,branch,t,x):
        pass
        
    def store_sens_step(self,branch,t,sA,sf,sGu,sGl):
        pass
        
