#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
import test_cases
import numpy as np
from scipy.sparse import coo_matrix,triu,tril

NUM_TRIALS = 25
EPS = 1e0 # %

class TestConstraints(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_constr_FIX(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             2*net.num_buses +
                             net.get_num_slack_gens() + 
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
            
            # Fixed
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_SLACK,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_REG_BY_GEN,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_FIXED,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_FIXED,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertGreater(net.num_fixed,0)
            self.assertEqual(net.num_fixed,
                             2*(net.get_num_slack_buses()) +
                             (net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            constr = pf.Constraint(pf.CONSTR_TYPE_FIX,net)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)

            Acounter = net.num_fixed+net.get_num_buses_reg_by_gen()
            constr.analyze()
            self.assertEqual(Acounter,constr.Acounter)
            constr.eval(x0)
            self.assertEqual(0,constr.Acounter)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(net.num_fixed,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(net.num_fixed,net.num_vars))
            self.assertEqual(A.nnz,net.num_fixed+net.get_num_buses_reg_by_gen())
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            
    def test_constr_PAR_GEN(self):
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            self.assertEqual(net.num_vars,0)
            
            # Vars
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          [pf.GEN_VAR_P,pf.GEN_VAR_Q])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_vars,net.get_num_slack_gens()+net.get_num_reg_gens())
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_PAR_GEN,net)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)

            # Manual count
            nnz = 0
            num_constr = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_slack():
                    num_constr += len(bus.gens)-1 # P participation
                    nnz += 2*(len(bus.gens)-1)
                if bus.is_regulated_by_gen():
                    num_constr += len(bus.reg_gens)-1 # Q participation
                    nnz += 2*(len(bus.reg_gens)-1)
            
            constr.analyze()
            self.assertEqual(nnz,constr.Acounter)
            constr.eval(x0)
            self.assertEqual(0,constr.Acounter)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(num_constr,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(num_constr,net.num_vars))
            self.assertEqual(A.nnz,nnz)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Detailed check
            Ai = A.row
            Aj = A.col
            Ad = A.data
            self.assertEqual(Ai.size,nnz)
            self.assertEqual(Aj.size,nnz)
            self.assertEqual(Ad.size,nnz)
            i = 0
            row = 0
            counted = {}
            for k in range(net.num_branches):
                br = net.get_branch(k)
                for bus in [br.bus_from,br.bus_to]:
                    if counted.has_key(bus.number):
                        continue
                    counted[bus.number] = True
                    if bus.is_slack():
                        gens = bus.gens
                        self.assertGreater(len(gens),0)
                        g1 = gens[0]
                        for g2 in gens[1:]:
                            self.assertEqual(b[row],0.)
                            self.assertEqual(Ai[i],row)
                            self.assertEqual(Aj[i],g1.index_P)
                            self.assertEqual(Ad[i],1.)
                            i += 1
                            self.assertEqual(Ai[i],row)
                            self.assertEqual(Aj[i],g2.index_P)
                            self.assertEqual(Ad[i],-1.)
                            i += 1
                            row += 1
                    if bus.is_regulated_by_gen():
                        reg_gens = bus.reg_gens
                        self.assertGreater(len(reg_gens),0)
                        g1 = reg_gens[0]
                        self.assertGreater(g1.Q_max,g1.Q_min)
                        for g2 in reg_gens[1:]:
                            self.assertTrue(np.abs(b[row]-(g1.Q_min/(g1.Q_max-g1.Q_min)-g2.Q_min/(g2.Q_max-g2.Q_min))) < 1e-10)
                            self.assertGreater(g2.Q_max,g2.Q_min)
                            self.assertEqual(Ai[i],row)
                            self.assertEqual(Aj[i],g1.index_Q)
                            self.assertTrue(np.abs(Ad[i]-1./(g1.Q_max-g1.Q_min)) < 1e-10)
                            i += 1
                            self.assertEqual(Ai[i],row)
                            self.assertEqual(Aj[i],g2.index_Q)
                            self.assertTrue(np.abs(Ad[i]+1./(g2.Q_max-g2.Q_min)) < 1e-10)
                            i += 1
                            row += 1
            self.assertEqual(i,nnz)

            # Last check
            x = np.zeros(net.num_vars)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_slack():
                    self.assertGreater(len(bus.gens),0)
                    for g in bus.gens:
                        self.assertTrue(g.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                        x[g.index_P] = 10.
                if bus.is_regulated_by_gen():
                    self.assertGreater(len(bus.reg_gens),0)
                    for g in bus.reg_gens:
                        self.assertTrue(g.has_flags(pf.FLAG_VARS,pf.GEN_VAR_Q))
                        x[g.index_Q] = (g.Q_max+g.Q_min)/2.
            self.assertGreater(np.linalg.norm(x),0)
            self.assertTrue(np.linalg.norm(A*x-b) < 1e-10)

    def test_constr_PF(self):
        
        # Constants
        h = 1e-10
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_vars,
                             2*net.get_num_buses() +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_PF,net)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)

            num_Jnnz = (net.get_num_buses()*4 +
                        net.get_num_branches()*8 +
                        net.get_num_tap_changers()*4 +
                        net.get_num_phase_shifters()*4 +
                        net.get_num_switched_shunts() +
                        net.get_num_slack_gens() +
                        net.get_num_reg_gens())
            
            constr.analyze()
            self.assertEqual(num_Jnnz,constr.Jcounter)
            constr.eval(x0)
            self.assertEqual(num_Jnnz,constr.Jcounter)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            constr.combine_H(np.ones(f.size),False)
            Hcomb = constr.H_combined            
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(2*net.num_buses,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(2*net.num_buses,net.num_vars))
            self.assertEqual(J.nnz,num_Jnnz)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            self.assertTupleEqual(Hcomb.shape,(net.num_vars,net.num_vars))
            self.assertEqual(Hcomb.nnz,2*(net.get_num_buses()*3 +
                                          net.get_num_branches()*12 +
                                          net.get_num_tap_changers()*9 +
                                          net.get_num_phase_shifters()*10 +
                                          net.get_num_switched_shunts()))
            
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))

            # Cross check mismatches
            net.update_properties(x0)
            dP_list = []
            dQ_list = []
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                dP = f[bus.index_P]
                dQ = f[bus.index_Q]
                dP_list.append(dP)
                dQ_list.append(dQ)
                self.assertLess(np.abs(dP-bus.P_mis),1e-10)
                self.assertLess(np.abs(dQ-bus.Q_mis),1e-10)
            self.assertLess(np.abs(net.bus_P_mis-np.max(np.abs(dP_list))*net.base_power),1e-10)
            self.assertLess(np.abs(net.bus_Q_mis-np.max(np.abs(dQ_list))*net.base_power),1e-10)

            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                f1 = constr.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.linalg.norm(Jd_exact)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):
                
                j = np.random.randint(0,f.shape[0])

                constr.eval(x0)
                
                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                
                H0 = (H0 + H0.T) - triu(H0)
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                
                g1 = constr.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
            sens = np.zeros(2*net.num_buses)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                sens[2*bus.index] = 3.5*(2*bus.index)+0.33
                sens[2*bus.index+1] = 3.5*(2*bus.index+1)+0.33
            constr.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,3.5*(2*bus.index)+0.33)
                self.assertEqual(bus.sens_Q_balance,3.5*(2*bus.index+1)+0.33)

    def test_constr_REG_GEN(self):
        
        # Constants
        h = 1e-8

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          [pf.BUS_PROP_NOT_SLACK,pf.BUS_PROP_REG_BY_GEN],
                          pf.BUS_VAR_VDEV)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-net.get_num_slack_buses()) +
                              2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) + 
                              net.get_num_slack_gens() + 
                              net.get_num_reg_gens()))
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_REG_GEN,net)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            
            #Jnnz = 4*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            Jnnz = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    Jnnz += 2 + 2*len(bus.reg_gens)
                    
            Annz = 3*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            self.assertGreater(Jnnz,0)
            self.assertGreater(Annz,0)
            
            rowsJ = 2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            rowsA = net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()
            self.assertGreater(rowsJ,0)
            self.assertGreater(rowsA,0)
                        
            constr.analyze()
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,Annz)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,rowsA)
            constr.eval(x0)
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,0)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ,net.num_vars))
            self.assertEqual(J.nnz,Jnnz)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA,net.num_vars))
            self.assertEqual(A.nnz,Annz)
            
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            # Ax=b check
            self.assertEqual(np.linalg.norm(A.data,1),rowsA*3)
            self.assertEqual(np.sum(A.data),net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            self.assertLess(np.linalg.norm(A*x0-b),1e-10)

            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                f1 = constr.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.linalg.norm(Jd_exact)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                j = np.random.randint(0,f.shape[0])

                constr.eval(x0)
                
                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                
                H0 = (H0 + H0.T) - triu(H0)
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                
                g1 = constr.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_v_reg_by_gen,0.)
            sens = np.zeros(constr.f.size)
            self.assertEqual(sens.size,rowsJ)
            Ji = constr.J.row
            Jj = constr.J.col
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    indices = Ji[np.where(np.logical_or(Jj == bus.index_y,
                                                        Jj == bus.index_z))[0]]
                    self.assertEqual(indices.size,2)
                    sens[indices[0]] = -bus.index-10
                    sens[indices[1]] = bus.index+11*(bus.index % 2)
            constr.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    if bus.index % 2 == 1:
                        self.assertEqual(bus.sens_v_reg_by_gen,bus.index+11)
                    else:
                        self.assertEqual(bus.sens_v_reg_by_gen,-bus.index-10)

    def test_constr_BOUND(self):
        
        # Constants
        h = 1e-8
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_vars,
                             2*net.num_buses +
                             2*net.num_gens +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())

            # Bound constraints
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_BOUNDED,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_BOUNDED,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_bounded,
                             2*net.num_buses +
                             2*net.num_gens +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_BOUND,net)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            
            Jcounter = 2*net.num_vars
            constr.analyze()
            self.assertEqual(Jcounter,constr.Jcounter)
            constr.eval(x0)
            self.assertEqual(Jcounter,constr.Jcounter)
            constr.store_sensitivities(np.zeros(0))
            self.assertEqual(Jcounter,constr.Jcounter)
            constr.eval(x0)
            self.assertEqual(Jcounter,constr.Jcounter)
                        
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(2*net.num_vars,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(2*net.num_vars,net.num_vars))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))

            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                f1 = constr.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.linalg.norm(Jd_exact)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):
                
                j = np.random.randint(0,f.shape[0])

                constr.eval(x0)
                
                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)
                
                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                self.assertEqual(H0.nnz,1)
                
                H0 = (H0 + H0.T) - triu(H0)
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),1e-4)
                self.assertLessEqual(error,EPS)
                
            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            self.assertEqual(H0.nnz,2*net.num_vars)
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                
                g1 = constr.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_v_mag_u_bound,0.)
                self.assertEqual(bus.sens_v_mag_l_bound,0.)
            sens = np.zeros(constr.f.size)
            Ji = constr.J.row
            Jj = constr.J.col
            for i in range(net.num_buses): # buses
                bus = net.get_bus(i)
                indices = Ji[np.where(Jj == bus.index_v_mag)[0]]
                self.assertEqual(indices.size,2)
                sens[indices[0]] = bus.index*10.
                sens[indices[1]] = -bus.index*10.
            constr.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_v_mag_u_bound,bus.index*10.)
                self.assertEqual(bus.sens_v_mag_l_bound,-bus.index*10.)
            
    def test_constr_REG_TRAN(self):
        
        # Constants
        h = 1e-8
        norm = 1e0
        eta = 1e-8

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VVIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO_DEV)
            self.assertEqual(net.num_vars,
                             3*net.get_num_buses_reg_by_tran() +
                             3*net.get_num_tap_changers_v())

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_REG_TRAN,net)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)

            Jnnz = 10*net.get_num_tap_changers_v()
            Annz = 3*net.get_num_tap_changers_v()
            self.assertGreaterEqual(Jnnz,0)
            self.assertGreaterEqual(Annz,0)
            
            rowsJ = 4*net.get_num_tap_changers_v()
            rowsA = net.get_num_tap_changers_v()
            self.assertGreaterEqual(rowsJ,0)
            self.assertGreaterEqual(rowsA,0)
                        
            constr.analyze()
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,Annz)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,rowsA)
            constr.eval(x0)
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ,net.num_vars))
            self.assertEqual(J.nnz,Jnnz)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA,net.num_vars))
            self.assertEqual(A.nnz,Annz)
            
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            # Ax=b check
            self.assertEqual(np.linalg.norm(A.data,1),rowsA*3)
            self.assertEqual(np.sum(A.data),net.get_num_tap_changers_v())
            self.assertLess(np.linalg.norm(A*x0-b),1e-10)
            
            # f check
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.is_tap_changer_v():
                    self.assertTrue(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
                    self.assertTrue(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO_DEV))
                    bus = br.reg_bus
                    fvmin = ((bus.v_mag-bus.v_min) - np.sqrt((bus.v_mag-bus.v_min)**2. + 2*eta))*norm
                    fvmax = ((bus.v_max-bus.v_mag) - np.sqrt((bus.v_max-bus.v_mag)**2. + 2*eta))*norm
                    ftmax = ((br.ratio_max-br.ratio) - np.sqrt((br.ratio_max-br.ratio)**2. + 2*eta))*norm
                    ftmin = ((br.ratio-br.ratio_min) - np.sqrt((br.ratio-br.ratio_min)**2. + 2*eta))*norm
                    self.assertLess(np.abs(fvmin-f[index]),1e-10)
                    self.assertLess(np.abs(fvmax-f[index+1]),1e-10)
                    self.assertLess(np.abs(ftmax-f[index+2]),1e-10)
                    self.assertLess(np.abs(ftmin-f[index+3]),1e-10)
                    index += 4
            
            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()

            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                f1 = constr.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),1e-4)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if not f.size:
                    continue

                j = np.random.randint(0,f.shape[0])
                
                constr.eval(x0)
                
                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                
                H0 = (H0 + H0.T) - triu(H0)
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                
                g1 = constr.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),1e-4)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_v_reg_by_tran,0.)
            sens = np.zeros(constr.f.size)
            Ji = constr.J.row
            Jj = constr.J.col
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_tran():
                    indices = Ji[np.where(np.logical_or(Jj == bus.index_vh,
                                                        Jj == bus.index_vl))[0]]
                    self.assertEqual(indices.size,4*len(bus.reg_trans))
                    j = np.random.randint(0,indices.size)                        
                    sens[indices[j]] = -bus.index - max([tran.index for tran in bus.reg_trans])
            constr.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_tran():
                    self.assertEqual(bus.sens_v_reg_by_tran,-bus.index-max([tran.index for tran in bus.reg_trans]))

    def test_constr_REG_SHUNT(self):
        
        # Constants
        h = 1e-8
        norm = 1e0
        eta = 1e-8

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_SHUNT,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_SHUNT,
                          pf.BUS_VAR_VVIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC_DEV)
            self.assertEqual(net.num_vars,
                             3*net.get_num_buses_reg_by_shunt() +
                             3*net.get_num_switched_shunts())

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_REG_SHUNT,net)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)

            Jnnz = 10*net.get_num_switched_shunts()
            Annz = 3*net.get_num_switched_shunts()
            self.assertGreaterEqual(Jnnz,0)
            self.assertGreaterEqual(Annz,0)
            
            rowsJ = 4*net.get_num_switched_shunts()
            rowsA = net.get_num_switched_shunts()
            self.assertGreaterEqual(rowsJ,0)
            self.assertGreaterEqual(rowsA,0)
                        
            constr.analyze()
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,Annz)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,rowsA)
            constr.eval(x0)
            self.assertEqual(constr.Jcounter,Jnnz)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,rowsJ)
            self.assertEqual(constr.Aconstr_index,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ,net.num_vars))
            self.assertEqual(J.nnz,Jnnz)
            self.assertTrue(np.all(J.row <= rowsJ-1))
            self.assertTrue(np.all(J.col <= net.num_vars-1))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA,net.num_vars))
            self.assertEqual(A.nnz,Annz)
            self.assertTrue(np.all(A.row <= rowsA-1))
            self.assertTrue(np.all(A.col <= net.num_vars-1))
            
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            # Ax=b check
            self.assertEqual(np.linalg.norm(A.data,1),rowsA*3)
            self.assertEqual(np.sum(A.data),net.get_num_switched_shunts())
            self.assertLess(np.linalg.norm(A*x0-b),1e-10)

            # f check
            index = 0
            counted = np.zeros(net.num_buses)
            for i in range(net.num_branches):
                br = net.get_branch(i)
                for bus in [br.bus_from,br.bus_to]:
                    if not counted[bus.index]:
                        for s in bus.reg_shunts:
                            self.assertEqual(bus.number,s.reg_bus.number)
                            self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))
                            self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VVIO))
                            self.assertTrue(s.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC))
                            self.assertTrue(s.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC_DEV))
                            self.assertEqual(x0[s.index_y],0.)
                            self.assertEqual(x0[s.index_z],0.)
                            self.assertEqual(x0[bus.index_vl],0.)
                            self.assertEqual(x0[bus.index_vh],0.)
                            fvmin = ((bus.v_mag-bus.v_min) - np.sqrt((bus.v_mag-bus.v_min)**2. + 2.*eta))*norm
                            fvmax = ((bus.v_max-bus.v_mag) - np.sqrt((bus.v_max-bus.v_mag)**2. + 2.*eta))*norm
                            fbmax = ((s.b_max-s.b) - np.sqrt((s.b_max-s.b)**2. + 2*eta))*norm
                            fbmin = ((s.b-s.b_min) - np.sqrt((s.b-s.b_min)**2. + 2*eta))*norm
                            self.assertLess(np.abs(fvmin-f[index]),1e-10)
                            self.assertLess(np.abs(fvmax-f[index+1]),1e-10)
                            self.assertLess(np.abs(fbmax-f[index+2]),1e-10)
                            self.assertLess(np.abs(fbmin-f[index+3]),1e-10)
                            index += 4
                    counted[bus.index] = 1
            
            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()

            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                f1 = constr.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),1e-4)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if not f.size:
                    continue

                j = np.random.randint(0,f.shape[0])
                
                constr.eval(x0)
                
                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                
                H0 = (H0 + H0.T) - triu(H0)
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                constr.eval(x)
                
                g1 = constr.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),1e-4)
                self.assertLessEqual(error,EPS)
            
            # Sensitivities
            net.clear_sensitivities()
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_v_reg_by_shunt,0.)
            sens = np.zeros(constr.f.size)
            Ji = constr.J.row
            Jj = constr.J.col
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_shunt():
                    indices = Ji[np.where(np.logical_or(Jj == bus.index_vh,
                                                        Jj == bus.index_vl))[0]]
                    self.assertEqual(indices.size,4*len(bus.reg_shunts))
                    j = np.random.randint(0,indices.size)           
                    sens[indices[j]] = -bus.index - max([s.index for s in bus.reg_shunts])
            constr.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_shunt():
                    self.assertEqual(bus.sens_v_reg_by_shunt,-bus.index-max([s.index for s in bus.reg_shunts]))

    def test_robustness(self):

        for case in test_cases.CASES:

            net = pf.Network()

            constraints = [pf.Constraint(pf.CONSTR_TYPE_BOUND,net),
                           pf.Constraint(pf.CONSTR_TYPE_FIX,net),
                           pf.Constraint(pf.CONSTR_TYPE_PAR_GEN,net),
                           pf.Constraint(pf.CONSTR_TYPE_PF,net),
                           pf.Constraint(pf.CONSTR_TYPE_REG_GEN,net),
                           pf.Constraint(pf.CONSTR_TYPE_REG_SHUNT,net),
                           pf.Constraint(pf.CONSTR_TYPE_REG_TRAN,net)]

            x0 = net.get_var_values()
        
            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))
                self.assertEqual(c.b.size,0)
                self.assertEqual(c.A.nnz,0)
                self.assertTupleEqual(c.A.shape,(0,0))
                self.assertEqual(c.f.size,0)
                self.assertEqual(c.J.nnz,0)
                self.assertTupleEqual(c.J.shape,(0,0))

            map(lambda c: c.eval(x0),constraints)
            map(lambda c: c.analyze(),constraints)
            map(lambda c: c.eval(x0),constraints)

            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))
                self.assertEqual(c.b.size,0)
                self.assertEqual(c.A.nnz,0)
                self.assertTupleEqual(c.A.shape,(0,0))
                self.assertEqual(c.f.size,0)
                self.assertEqual(c.J.nnz,0)
                self.assertTupleEqual(c.J.shape,(0,0))

            # Network changes
            net.load(case)

            # Before updating network
            map(lambda c: c.clear_error(),constraints)
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            map(lambda c: c.clear_error(),constraints)
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            map(lambda c: c.clear_error(),constraints)
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.analyze)
            map(lambda c: c.clear_error(),constraints)
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            map(lambda c: c.clear_error(),constraints)

            # Update network
            map(lambda c: c.update_network(),constraints)
            
            # After updating network
            map(lambda c: c.analyze(),constraints)
            map(lambda c: c.eval(x0),constraints)

            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))
                self.assertEqual(c.A.nnz,0)
                self.assertEqual(c.J.nnz,0)

            # Add variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_vars,
                             (2*net.num_buses + 
                              2*net.num_gens +
                              net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts()))

            x0 = net.get_var_values()

            # Before analyzing
            map(lambda c: c.clear_error(),constraints)
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            map(lambda c: c.clear_error(),constraints)

            # Do it right
            map(lambda c: c.analyze(),constraints)
            map(lambda c: c.eval(x0),constraints)
            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))
                self.assertEqual(c.A.shape[1],net.num_vars)
                self.assertEqual(c.J.shape[1],net.num_vars)
                if c.f.size:
                    self.assertTupleEqual(c.get_H_single(0).shape,(net.num_vars,net.num_vars))
                else:
                    self.assertTupleEqual(c.get_H_single(0).shape,(0,0))

    def test_constr_DCPF(self):
        
        # Constants
        h = 1e-10
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_vargens,0)
            
            # Add vargens
            load_buses = net.get_load_buses()
            vargen_array = pf.VarGeneratorArray(len(load_buses))
            net.set_vargen_array(vargen_array)
            net.set_vargen_buses(load_buses)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len([b for b in net.buses if b.loads]))
            for b in net.buses:
                if b.loads:
                    self.assertGreater(len(b.vargens),0)
                    for vargen in b.vargens:
                        self.assertEqual(vargen.bus,b)

            # Variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            self.assertEqual(net.num_vars,
                             (net.num_buses-net.get_num_slack_buses() +
                              net.num_gens + 
                              net.num_vargens +
                              net.get_num_phase_shifters()))
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_DCPF,net)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b

            # Before 
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)

            r = 0
            for b in net.buses:
                if b.is_slack():
                    r += len(b.branches)
            
            # Analyze
            constr.analyze()
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            self.assertEqual(constr.Acounter,
                             (net.num_gens +
                              net.num_vargens +
                              4*net.num_branches - 
                              2*r +
                              2*net.get_num_phase_shifters()))
            self.assertTupleEqual(b.shape,(net.num_buses,))
            self.assertTupleEqual(f.shape,(0,))
            self.assertTupleEqual(A.shape,(net.num_buses,net.num_vars))
            self.assertEqual(A.nnz,constr.Acounter)
            self.assertTupleEqual(J.shape,(0,net.num_vars))            
            
            constr.eval(x0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(A.nnz,
                             (net.num_gens +
                              net.num_vargens +
                              4*net.num_branches - 
                              2*r +
                              2*net.get_num_phase_shifters()))
            
            # Extract pieces
            P1 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)
            P2 = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P)
            P3 = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)
            P4 = net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_PHASE)
            
            G = A*P2.T
            R = A*P3.T
            Atheta = -A*P1.T
            Aphi = -A*P4.T
            x = np.random.randn(net.num_vars)
            p = P2*x
            r = P3*x
            theta = P1*x
            phi = P4*x
            self.assertLess(np.linalg.norm((G*p+R*r-Atheta*theta-Aphi*phi)-A*x),1e-10)            

    def tearDown(self):
        
        pass
