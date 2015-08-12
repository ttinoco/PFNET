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
from scipy.sparse import coo_matrix,triu

NUM_TRIALS = 25
EPS = 2e0 # %

class TestProblem(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()
        self.p = pf.Problem()

    def test_problem_LSNR(self):

        # Constants
        h = 1e-9
        
        p = self.p
        net = self.net

        for case in test_cases.CASES:
            
            p.clear()
            net.load(case)
            p.network = net
            
            # Variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
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
                          pf.BRANCH_PROP_TAP_CHANGER_V,
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
                             2*(net.num_buses-net.get_num_slack_buses()) +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers_v() + 
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
                             
            # Fixed
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_REG_BY_GEN,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_FIXED,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_fixed,
                             net.get_num_buses_reg_by_gen() +
                             net.get_num_tap_changers_v() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
                             
            # Constraints
            p.add_constraint(pf.CONSTR_TYPE_PF)
            p.add_constraint(pf.CONSTR_TYPE_PAR_GEN_P)
            p.add_constraint(pf.CONSTR_TYPE_PAR_GEN_Q)
            p.add_constraint(pf.CONSTR_TYPE_FIX)
            self.assertEqual(len(p.constraints),4)

            # Check adding redundant constraints
            p.add_constraint(pf.CONSTR_TYPE_PAR_GEN_P)
            self.assertEqual(len(p.constraints),4)
            
            # Functions
            self.assertEqual(len(p.functions),0)
                
            # Init point
            x0 = p.get_init_point()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Before
            phi = p.phi
            gphi = p.gphi
            Hphi = p.Hphi

            f = p.f
            b = p.b
            A = p.A
            Z = p.Z
            J = p.J
            
            self.assertTrue(type(phi) is float)
            self.assertEqual(phi,0.)
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(0,))
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
            self.assertTrue(type(Z) is coo_matrix)
            self.assertTupleEqual(Z.shape,(0,0))
            self.assertEqual(Z.nnz,0)
            self.assertTrue(type(Hphi) is coo_matrix)
            self.assertTupleEqual(Hphi.shape,(0,0))
            self.assertEqual(Hphi.nnz,0)
            self.assertTrue(np.all(Hphi.row >= Hphi.col))
            
            p.analyze()
            p.eval(x0)
            
            # After
            phi = p.phi
            gphi = p.gphi.copy()
            Hphi = p.Hphi.copy()

            f = p.f.copy()
            b = p.b.copy()
            A = p.A.copy()
            J = p.J.copy()
            Z = p.Z
                        
            # phi
            self.assertTrue(type(phi) is float)
            self.assertEqual(phi,0.)
            
            # gphi
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(net.num_vars,))
            self.assertLess(np.linalg.norm(gphi),1e-10)

            # Hphi
            self.assertTrue(type(Hphi) is coo_matrix)
            self.assertTupleEqual(Hphi.shape,(net.num_vars,net.num_vars))
            self.assertEqual(Hphi.nnz,0)
    
            # f
            self.assertTrue(type(f) is np.ndarray)
            f_size = sum(c.f.shape[0] for c in p.constraints)
            self.assertTupleEqual(f.shape,(f_size,))

            # b
            self.assertTrue(type(b) is np.ndarray)
            b_size = sum(c.b.shape[0] for c in p.constraints)
            self.assertTupleEqual(b.shape,(b_size,))

            # J
            self.assertTrue(type(J) is coo_matrix)
            J_size = sum(c.J.shape[0] for c in p.constraints)
            self.assertTupleEqual(J.shape,(J_size,net.num_vars))
            self.assertGreater(J.nnz,0)
            
            # A
            self.assertTrue(type(A) is coo_matrix)
            A_size = sum(c.A.shape[0] for c in p.constraints)
            self.assertTupleEqual(A.shape,(A_size,net.num_vars))
            self.assertGreater(A.nnz,0)

            # Z            
            Znnz = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if not bus.is_slack():
                    Znnz += 1 # w
                if not bus.is_slack() and not bus.is_regulated_by_gen():
                    Znnz += 1 # v
                if bus.is_slack():
                    for g in bus.gens:
                        Znnz += 1 # P
                if bus.is_regulated_by_gen():
                    for g in bus.reg_gens:
                        Znnz += 1 # Q
            self.assertTrue(type(Z) is coo_matrix)
            self.assertTupleEqual(Z.shape,(net.num_vars,net.num_vars-A.shape[0]))
            self.assertGreater(Z.nnz,0)
            self.assertEqual(Z.nnz,Znnz)
            self.assertLess(np.linalg.norm((A*Z).data),1e-10)
            
            # Check J
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                p.eval(x)
                f1 = p.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.linalg.norm(Jd_exact)
                self.assertLessEqual(error,EPS)
                
            # Check Hcombined
            coeff = np.random.randn(f.shape[0])
            p.eval(x0)
            p.combine_H(coeff,False)
            J0 = p.J.copy()
            g0 = J0.T*coeff
            H0 = p.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
                
                x = x0 + h*d
                
                p.eval(x)
                
                g1 = p.J.T*coeff
                
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
            sens = np.random.randn(p.f.size)
            offset = 0
            for c in p.constraints:
                if c.type == pf.CONSTR_TYPE_PF:
                    break
                else:
                    offset += c.f.size
            p.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,sens[2*bus.index+offset])
                self.assertEqual(bus.sens_Q_balance,sens[2*bus.index+1+offset])
            self.assertRaises(pf.ProblemError,p.store_sensitivities,np.zeros(p.f.size+5))

    def test_problem_vPF(self):

        # Constants
        h = 1e-9

        p = self.p
        net = self.net

        for case in test_cases.CASES:
            
            p.clear()
            net.load(case)
            p.network = net
            
            # Variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK|pf.BUS_PROP_REG_BY_GEN,
                          pf.BUS_VAR_VDEV)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VVIO)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_SHUNT,
                          pf.BUS_VAR_VVIO)
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
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV)                          

            reg_by_tran_or_shunt = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_tran() or bus.is_regulated_by_shunt():
                    reg_by_tran_or_shunt += 1
            
            self.assertEqual(net.num_vars,
                             2*(net.num_buses-net.get_num_slack_buses()) + 
                             2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                             2*(reg_by_tran_or_shunt) + 
                             net.get_num_slack_gens() + 
                             net.get_num_reg_gens() + 
                             3*net.get_num_tap_changers_v()+
                             3*net.get_num_switched_shunts())
                             
            # Constraints
            p.add_constraint(pf.CONSTR_TYPE_PF)
            p.add_constraint(pf.CONSTR_TYPE_PAR_GEN_P)
            p.add_constraint(pf.CONSTR_TYPE_PAR_GEN_Q)
            p.add_constraint(pf.CONSTR_TYPE_REG_GEN)
            p.add_constraint(pf.CONSTR_TYPE_REG_TRAN)
            p.add_constraint(pf.CONSTR_TYPE_REG_SHUNT)
            self.assertEqual(len(p.constraints),6)

            # Check adding redundant constraints
            p.add_constraint(pf.CONSTR_TYPE_PF)
            self.assertEqual(len(p.constraints),6)
            
            # Functions
            p.add_function(pf.FUNC_TYPE_REG_VMAG,1.)
            p.add_function(pf.FUNC_TYPE_REG_VANG,5.)
            p.add_function(pf.FUNC_TYPE_REG_PQ,8.)
            p.add_function(pf.FUNC_TYPE_REG_RATIO,3.)            
            p.add_function(pf.FUNC_TYPE_REG_SUSC,1.)
            self.assertEqual(len(p.functions),5)
                
            # Init point
            x0 = p.get_init_point()+np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Before
            phi = p.phi
            gphi = p.gphi
            Hphi = p.Hphi

            f = p.f
            b = p.b
            A = p.A
            Z = p.Z
            J = p.J
            
            self.assertTrue(type(phi) is float)
            self.assertEqual(phi,0.)
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(0,))
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
            self.assertTrue(type(Z) is coo_matrix)
            self.assertTupleEqual(Z.shape,(0,0))
            self.assertEqual(Z.nnz,0)
            self.assertTrue(type(Hphi) is coo_matrix)
            self.assertTupleEqual(Hphi.shape,(0,0))
            self.assertEqual(Hphi.nnz,0)
            self.assertTrue(np.all(Hphi.row >= Hphi.col))
            
            p.analyze()
            p.eval(x0)
            
            # After
            phi = p.phi
            gphi = p.gphi.copy()
            Hphi = p.Hphi.copy()

            f = p.f.copy()
            b = p.b.copy()
            A = p.A.copy()
            J = p.J.copy()
            Z = p.Z
                        
            # phi
            self.assertTrue(type(phi) is float)
            self.assertGreater(phi,0.)
            man_phi = sum(f.weight*f.phi for f in p.functions)
            self.assertLess(np.abs(man_phi-phi),1e-10)

            # gphi
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(net.num_vars,))
            man_gphi = sum(f.weight*f.gphi for f in p.functions)
            self.assertLess(np.linalg.norm(man_gphi-gphi),1e-10)

            # Hphi
            self.assertTrue(type(Hphi) is coo_matrix)
            self.assertTupleEqual(Hphi.shape,(net.num_vars,net.num_vars))
            self.assertGreater(Hphi.nnz,0)
    
            # f
            self.assertTrue(type(f) is np.ndarray)
            f_size = sum(c.f.shape[0] for c in p.constraints)
            self.assertTupleEqual(f.shape,(f_size,))

            # b
            self.assertTrue(type(b) is np.ndarray)
            b_size = sum(c.b.shape[0] for c in p.constraints)
            self.assertTupleEqual(b.shape,(b_size,))

            # J
            self.assertTrue(type(J) is coo_matrix)
            J_size = sum(c.J.shape[0] for c in p.constraints)
            self.assertTupleEqual(J.shape,(J_size,net.num_vars))
            self.assertGreater(J.nnz,0)
            
            # A
            self.assertTrue(type(A) is coo_matrix)
            A_size = sum(c.A.shape[0] for c in p.constraints)
            self.assertTupleEqual(A.shape,(A_size,net.num_vars))
            self.assertGreater(A.nnz,0)

            # Z            
            Znnz = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if not bus.is_slack():
                    Znnz += 1 # w
                if not bus.is_slack() and not bus.is_regulated_by_gen():
                    Znnz += 1 # v
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    Znnz += 2 # vy
                    Znnz += 2 # vz
                if bus.is_regulated_by_tran() or bus.is_regulated_by_shunt():
                    Znnz += 1 # vl
                    Znnz += 1 # vh
                if bus.is_slack():
                    for g in bus.gens:
                        Znnz += 1
                if bus.is_regulated_by_gen():
                    for g in bus.reg_gens:
                        Znnz += 1
            Znnz += 4*net.get_num_tap_changers_v()
            Znnz += 4*net.get_num_switched_shunts()
            self.assertTrue(type(Z) is coo_matrix)
            self.assertTupleEqual(Z.shape,(net.num_vars,net.num_vars-A.shape[0]))
            self.assertGreater(Z.nnz,0)
            self.assertEqual(Z.nnz,Znnz)
            self.assertLess(np.linalg.norm((A*Z).data),1e-10)
            
            # Check gphi
            phi0 = phi
            gphi0 = gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                p.eval(x)
                phi1 = p.phi
                
                gd_exact = np.dot(gphi0,d)
                gd_approx = (phi1-phi0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.linalg.norm(gd_exact)
                self.assertLessEqual(error,EPS)

            # Check J
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                p.eval(x)
                f1 = p.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.linalg.norm(Jd_exact)
                self.assertLessEqual(error,EPS)
                
            # Check Hphi
            gphi0 = gphi.copy()
            Hphi0 = Hphi.copy()
            Hphi0 = Hphi0 + Hphi0.T - triu(Hphi0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                p.eval(x)

                gphi1 = p.gphi.copy()
                
                Hd_exact = Hphi0*d
                Hd_approx = (gphi1-gphi0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.linalg.norm(Hd_exact)
                self.assertLessEqual(error,EPS)

            # Check Hcombined
            coeff = np.random.randn(f.shape[0])
            p.eval(x0)
            p.combine_H(coeff,False)
            J0 = p.J.copy()
            g0 = J0.T*coeff
            H0 = p.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars,net.num_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
                
                x = x0 + h*d
                
                p.eval(x)
                
                g1 = p.J.T*coeff
                
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
            sens = np.random.randn(p.f.size)
            offset = 0
            for c in p.constraints:
                if c.type == pf.CONSTR_TYPE_PF:
                    break
                else:
                    offset += c.f.size
            p.store_sensitivities(sens)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,sens[2*bus.index+offset])
                self.assertEqual(bus.sens_Q_balance,sens[2*bus.index+1+offset])
            self.assertRaises(pf.ProblemError,p.store_sensitivities,np.zeros(p.f.size+5))

    def test_problem_Z_shunts(self):

        p = self.p
        net = self.net

        for case in test_cases.CASES:
            
            p.clear()
            net.load(case)
            p.network = net

            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV)

            p.analyze()
            
            Z = p.Z

            self.assertTupleEqual(Z.shape,(3*net.get_num_switched_shunts(),
                                           3*net.get_num_switched_shunts()))
            self.assertTrue(np.all(Z.data == 1))
            self.assertTrue(np.all(Z.row == np.array(range(3*net.get_num_switched_shunts()))))
            self.assertTrue(np.all(Z.col == np.array(range(3*net.get_num_switched_shunts()))))

    def tearDown(self):
        
        pass
