#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
from . import test_cases
import numpy as np
from scipy.sparse import coo_matrix,triu,tril

NUM_TRIALS = 25
EPS = 2e0 # %
TOL = 1e-4

class TestFunctions(unittest.TestCase):

    def setUp(self):
        
        # Network
        self.net = pf.Network()

        # Random
        np.random.seed(1)
        
    def test_func_REG_VMAG(self):
        
        # Constants
        h = 1e-9
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            nb = net.num_buses
            nrg = net.get_num_buses_reg_by_gen()
            ns = net.get_num_slack_buses()
            nrt = net.get_num_buses_reg_by_tran()
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          [pf.BUS_PROP_REG_BY_GEN,pf.BUS_PROP_NOT_SLACK],
                          pf.BUS_VAR_VDEV)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VVIO)
            self.assertEqual(net.num_vars,2*nb+2*(nrg-ns)+2*nrt)
             
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_VMAG,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_VMAG)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,nb+2*(nrg-ns)+2*nrt)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,nb+2*(nrg-ns)+2*nrt)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                if np.linalg.norm(gd_exact) == 0.:
                    self.assertLessEqual(np.linalg.norm(gd_approx),2.)
                else:
                    error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                    self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_REG_PQ(self):
        
        # Constants
        h = 1e-9
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            ng = net.num_gens
            nrg = net.get_num_reg_gens()
            ns = net.get_num_slack_gens()
            
            # Vars
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            self.assertEqual(net.num_vars,ns+nrg)
             
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_PQ,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_PQ)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,nrg+ns)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,nrg+ns)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h

                if np.linalg.norm(gd_exact) == 0.:
                    self.assertLessEqual(np.linalg.norm(gd_approx),2.)
                else:
                    error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                    self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_REG_VANG(self):
        
        # Constants
        h = 1e-8
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            nb = net.num_buses
            ns = net.get_num_slack_buses()
            
            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            self.assertEqual(net.num_vars,2*(nb-ns))
             
            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_VANG,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_VANG)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)

            # Manual count
            man_Hcounter = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                bk = br.bus_from
                bm = br.bus_to
                if not bk.is_slack():
                    man_Hcounter +=1
                if not bm.is_slack():
                    man_Hcounter +=1
                if not bk.is_slack() and not bm.is_slack():
                    man_Hcounter += 1
            man_Hcounter += nb-ns
            
            func.analyze()
            self.assertEqual(func.Hcounter,man_Hcounter)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,man_Hcounter)
            self.assertTrue(np.all(H.row >= H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            H0 = H0 + H0.T - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_REG_RATIO(self):
        
        # Constants
        h = 1e-8
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV)
            self.assertEqual(net.num_vars,3*net.get_num_tap_changers_v())
             
            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_RATIO,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_RATIO)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,3*net.get_num_tap_changers_v())
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,3*net.get_num_tap_changers_v())
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.is_tap_changer_v():
                    self.assertTrue(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
                    t = x0[br.index_ratio]
                    tmax = br.ratio_max
                    tmin = br.ratio_min
                    t0 = br.ratio
                    dt = np.maximum(tmax-tmin,1e-4)
                    phi_manual += 0.5*((t-t0)/dt)**2.
                    phi_manual += 0.5*(x0[br.index_ratio_y]/dt)**2.
                    phi_manual += 0.5*(x0[br.index_ratio_z]/dt)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-9)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            H0 = H0 + H0.T - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_REG_SUSC(self):
        
        # Constants
        h = 1e-8
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV)
            self.assertEqual(net.num_vars,3*net.get_num_switched_shunts())
             
            x0 = net.get_var_values() + np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_SUSC,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_SUSC)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,3*net.get_num_switched_shunts())
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,3*net.get_num_switched_shunts())
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for i in range(net.num_shunts):
                sh = net.get_shunt(i)
                if sh.is_switched_v():
                    self.assertTrue(sh.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC))
                    b = x0[sh.index_b]
                    bmax = sh.b_max
                    bmin = sh.b_min
                    b0 = sh.b
                    db = np.maximum(bmax-bmin,1e-4)
                    phi_manual += 0.5*((b-b0)/db)**2.
                    phi_manual += 0.5*(x0[sh.index_y]/db)**2.
                    phi_manual += 0.5*(x0[sh.index_z]/db)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-9)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            H0 = H0 + H0.T - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_GEN_COST(self):
        
        # Constants
        h = 1e-9
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          [pf.GEN_VAR_P,pf.GEN_VAR_Q])
            self.assertEqual(net.num_vars,2*net.get_num_gens())
            self.assertGreater(net.num_vars,0)
             
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_GEN_COST,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_GEN_COST)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,net.get_num_gens())
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,net.get_num_gens())
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))
            
            # value check
            val = 0
            for gen in net.generators:
                self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                val += (gen.cost_coeff_Q0 + 
                        gen.cost_coeff_Q1*gen.P +
                        gen.cost_coeff_Q2*(gen.P**2.))
            self.assertLess(np.abs(val-f),1e-8)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_SP_CONTROLS(self):
        
        # Constants
        h = 1e-9
        
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
                          pf.GEN_VAR_P)
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
                              net.num_gens +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()))

            # Sparse
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_SPARSE,
                          pf.BUS_PROP_ANY,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_SPARSE,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_SPARSE,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_SPARSE,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_SPARSE,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_sparse,
                             (2*net.num_buses +
                              net.num_gens +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()))
             
            x0 = net.get_var_values() + np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_SP_CONTROLS,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_SP_CONTROLS)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)

            Hcounter_manual = (net.get_num_buses_reg_by_gen() +
                               net.get_num_gens() + 
                               net.get_num_tap_changers() + 
                               net.get_num_phase_shifters() + 
                               net.get_num_switched_shunts())
            
            func.analyze()
            self.assertEqual(func.Hcounter,Hcounter_manual)
            func.eval(x0)
            self.assertEqual(func.Hcounter,Hcounter_manual)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreater(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,Hcounter_manual)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual f value
            eps = 1e-6
            ceps = 1e-4
            f_manual = 0
            for branch in net.branches:
                if branch.is_tap_changer():
                    self.assertTrue(branch.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
                    val = x0[branch.index_ratio]
                    val0 = branch.ratio
                    dval = np.maximum(branch.ratio_max-branch.ratio_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_PHASE))
                    val = x0[branch.index_phase]
                    val0 = branch.phase
                    dval = np.maximum(branch.phase_max-branch.phase_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))
                    val = x0[bus.index_v_mag]
                    val0 = bus.v_set
                    dval = np.maximum(bus.v_max-bus.v_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for gen in net.generators:
                self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                val = x0[gen.index_P]
                val0 = gen.P
                dval = np.maximum(gen.P_max-gen.P_min,ceps)
                f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC))
                    val = x0[shunt.index_b]
                    val0 = shunt.b
                    dval = np.maximum(shunt.b_max-shunt.b_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            self.assertLess(np.abs(f-f_manual),1e-8)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            h = 1e-6
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_SLIM_VMAG(self):
        
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
            self.assertEqual(net.num_vars,2*net.num_buses)
             
            x0 = net.get_var_values() + 1.*np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_SLIM_VMAG,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_SLIM_VMAG)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,net.num_buses)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreater(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,net.num_buses)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # f manual
            f_manual = 0
            eps = 1e-4
            for bus in net.buses:
                self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))
                dv = np.maximum(bus.v_max-bus.v_min,eps)
                vmid = 0.5*(bus.v_max+bus.v_min)
                f_manual += 0.5*(((x0[bus.index_v_mag]-vmid)/dv)**2.)
            self.assertLess(np.abs(f-f_manual)/f_manual,1e-4)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                if np.linalg.norm(gd_exact) == 0.:
                    self.assertLessEqual(np.linalg.norm(gd_approx),2.)
                else:
                    error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_REG_PHASE(self):
        
        # Constants
        h = 1e-8
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            self.assertEqual(net.num_vars,net.get_num_phase_shifters())
             
            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_REG_PHASE,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_REG_PHASE)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,net.get_num_phase_shifters())
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertGreaterEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,net.get_num_phase_shifters())
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.is_phase_shifter():
                    self.assertTrue(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_PHASE))
                    p = x0[br.index_phase]
                    pmax = br.phase_max
                    pmin = br.phase_min
                    p0 = br.phase
                    dp = np.maximum(pmax-pmin,1e-4)
                    phi_manual += 0.5*((p-p0)/dp)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-9)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            H0 = H0 + H0.T - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_LOAD_UTIL(self):
        
        # Constants
        h = 1e-9
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            # Vars
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            self.assertEqual(net.num_vars,net.num_loads)
            self.assertGreater(net.num_vars,0)
             
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_LOAD_UTIL,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_LOAD_UTIL)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,net.num_loads)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertNotEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,net.num_loads)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))
            
            # value check
            val = 0
            for load in net.loads:
                self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                val += (load.util_coeff_Q0 + 
                        load.util_coeff_Q1*load.P +
                        load.util_coeff_Q2*(load.P**2.))
            net.update_properties()
            self.assertLess(np.abs(val-f),1e-7)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Hessian check
            func.eval(x0)
            g0 = func.gphi.copy()
            H0 = func.Hphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)

                g1 = func.gphi.copy()
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_func_NETCON_COST(self):
        
        # Constants
        h = 1e-9
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            # prices
            for bus in net.buses:
                bus.price = (bus.index%10)*0.5123

            # vargens
            net.add_vargens(net.get_load_buses(),50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = (vargen.index%10)*0.3233+0.1

            # batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P *= -1.
            
            # Vars
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P)
            self.assertEqual(net.num_vars,
                             (net.num_loads+
                              net.num_gens+
                              net.num_vargens+
                              2*net.num_bats))
            self.assertGreater(net.num_vars,0)
             
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Function
            func = pf.Function(pf.FUNC_TYPE_NETCON_COST,1.,net)

            self.assertEqual(func.type,pf.FUNC_TYPE_NETCON_COST)

            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # Before 
            self.assertTrue(type(f) is float)
            self.assertEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(0,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(0,0))
            self.assertEqual(H.nnz,0)
            
            self.assertEqual(func.Hcounter,0)
            
            func.analyze()
            self.assertEqual(func.Hcounter,0)
            func.eval(x0)
            self.assertEqual(func.Hcounter,0)
            
            f = func.phi
            g = func.gphi
            H = func.Hphi
            
            # After
            self.assertTrue(type(f) is float)
            self.assertNotEqual(f,0.)
            self.assertTrue(type(g) is np.ndarray)
            self.assertTupleEqual(g.shape,(net.num_vars,))
            self.assertTrue(type(H) is coo_matrix)
            self.assertTupleEqual(H.shape,(net.num_vars,net.num_vars))
            self.assertEqual(H.nnz,0)

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            
            # value check
            val = 0
            for bus in net.buses:
                for load in bus.loads:
                    self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                    val += bus.price*load.P
                for bat in bus.bats:
                    self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P))
                    val += bus.price*bat.P
                for gen in bus.gens:
                    self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                    val -= bus.price*gen.P
                for vargen in bus.vargens:
                    self.assertTrue(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P))
                    val -= bus.price*vargen.P
            self.assertLess(100*np.abs(val-f)/np.abs(f),1e-8)

            # Gradient check
            f0 = func.phi
            g0 = func.gphi.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars)
    
                x = x0 + h*d
                
                func.eval(x)
                f1 = func.phi
                
                gd_exact = np.dot(g0,d)
                gd_approx = (f1-f0)/h
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # No variables
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            
            func = pf.Function(pf.FUNC_TYPE_NETCON_COST,1.,net)
            self.assertEqual(func.type,pf.FUNC_TYPE_NETCON_COST)

            x0 = net.get_var_values()

            func.analyze()
            func.eval(x0)
            
            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))
            
            # value check
            val = 0
            for bus in net.buses:
                for load in bus.loads:
                    self.assertFalse(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                    val += bus.price*load.P
                for bat in bus.bats:
                    self.assertFalse(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P))
                    val += bus.price*bat.P
                for gen in bus.gens:
                    self.assertFalse(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                    val -= bus.price*gen.P
                for vargen in bus.vargens:
                    self.assertFalse(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P))
                    val -= bus.price*vargen.P
            self.assertLess(100*np.abs(val-func.phi)/np.abs(f),1e-8)
            
    def test_robustness(self):

        for case in test_cases.CASES:

            net = pf.Network()
            
            functions = [pf.Function(pf.FUNC_TYPE_GEN_COST,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_PHASE,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_PQ,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_RATIO,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_SUSC,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_VANG,1.,net),
                         pf.Function(pf.FUNC_TYPE_REG_VMAG,1.,net),
                         pf.Function(pf.FUNC_TYPE_SLIM_VMAG,1.,net),
                         pf.Function(pf.FUNC_TYPE_SP_CONTROLS,1.,net),
                         pf.Function(pf.FUNC_TYPE_LOAD_UTIL,1.,net)]
            
            x0 = net.get_var_values()
        
            for f in functions:
                self.assertEqual(f.phi,0.)
                self.assertTrue(isinstance(f.gphi,np.ndarray))
                self.assertTrue(isinstance(f.Hphi,coo_matrix))
                self.assertEqual(f.gphi.size,0)
                self.assertEqual(f.Hphi.nnz,0)

            list(map(lambda f: f.eval(x0),functions))
            list(map(lambda f: f.analyze(),functions))
            list(map(lambda f: f.eval(x0),functions))

            for f in functions:
                self.assertEqual(f.phi,0.)
                self.assertTrue(isinstance(f.gphi,np.ndarray))
                self.assertTrue(isinstance(f.Hphi,coo_matrix))
                self.assertEqual(f.gphi.size,0)
                self.assertEqual(f.Hphi.nnz,0)

            # Network changes
            net.load(case)

            # Before updating network
            list(map(lambda f: f.clear_error(),functions))
            for f in functions:
                self.assertRaises(pf.FunctionError,f.eval,x0)
            list(map(lambda f: f.clear_error(),functions))
            for f in functions:
                self.assertRaises(pf.FunctionError,f.eval,x0)
            list(map(lambda f: f.clear_error(),functions))
            for f in functions:
                self.assertRaises(pf.FunctionError,f.analyze)
            list(map(lambda f: f.clear_error(),functions))
            for f in functions:
                self.assertRaises(pf.FunctionError,f.eval,x0)
            list(map(lambda f: f.clear_error(),functions))

            # Update network
            list(map(lambda f: f.update_network(),functions))
            
            # After updating network
            list(map(lambda f: f.analyze(),functions))
            list(map(lambda f: f.eval(x0),functions))

            for f in functions:
                self.assertEqual(f.phi,0.)
                self.assertTrue(isinstance(f.gphi,np.ndarray))
                self.assertTrue(isinstance(f.Hphi,coo_matrix))
                self.assertEqual(f.gphi.size,0)
                self.assertEqual(f.Hphi.nnz,0)

            # Add variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
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
                              net.num_loads +
                              net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts()))

            x0 = net.get_var_values()

            # Before analyzing
            list(map(lambda f: f.clear_error(),functions))
            for f in functions:
                self.assertRaises(pf.FunctionError,f.eval,x0)
            list(map(lambda f: f.clear_error(),functions))

            # Do it right
            list(map(lambda f: f.analyze(),functions))
            list(map(lambda f: f.eval(x0),functions))
            for f in functions:
                self.assertTrue(isinstance(f.gphi,np.ndarray))
                self.assertTrue(isinstance(f.Hphi,coo_matrix))
                self.assertEqual(f.gphi.size,net.num_vars)
                self.assertTupleEqual(f.Hphi.shape,(net.num_vars,net.num_vars))
                self.assertGreaterEqual(f.Hphi.nnz,0)
            self.assertTrue(any([f.phi > 0 for f in functions]))
            self.assertTrue(any([f.Hphi.nnz > 0 for f in functions]))
            
    def tearDown(self):
        
        pass
