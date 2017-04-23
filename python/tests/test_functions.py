#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
from . import test_cases
import numpy as np
from scipy.sparse import coo_matrix,triu,tril

NUM_TRIALS = 25
EPS = 5. # %
TOL = 1e-3

class TestFunctions(unittest.TestCase):

    def setUp(self):

        # Network
        self.T = 4

        # Random
        np.random.seed(1)

    def test_func_REG_VMAG(self):

        # Constants
        h = 1e-9

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            nb = net.num_buses
            nrg = net.get_num_buses_reg_by_gen()
            ns = net.get_num_slack_buses()
            nrt = net.get_num_buses_reg_by_tran()

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('bus',
                          'variable',
                          'regulated by transformer',
                          'voltage magnitude violation')
            self.assertEqual(net.num_vars,(2*nb+2*nrt)*net.num_periods)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Perturbation
            net.set_var_values(x0 + np.random.randn(x0.size))
            x0 = net.get_var_values()

            # Function
            func = pf.Function('voltage magnitude regularization',1.,net)

            self.assertEqual(func.name,'voltage magnitude regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,(nb+2*nrt)*net.num_periods)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,(nb+2*nrt)*net.num_periods)
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

            # Value
            dv = 0.2
            func.eval(x0)
            phi = 0
            for t in range(net.num_periods):
                for bus in net.buses:
                    phi += 0.5*(((bus.v_mag[t]-bus.v_set[t])/dv)**2.)
            self.assertLess(np.abs(func.phi-phi),1e-10*(func.phi+1))
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            phi = 0
            for t in range(net.num_periods):
                for bus in net.buses:
                    phi += 0.5*(((bus.v_mag[t]-bus.v_set[t])/dv)**2.)
            self.assertLess(np.abs(func.phi-phi),1e-10*(func.phi+1))

    def test_func_REG_PQ(self):

        # Constants
        h = 1e-9

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            ng = net.num_generators
            nrg = net.get_num_reg_gens()
            ns = net.get_num_slack_gens()

            # Vars
            net.set_flags('generator',
                          'variable',
                          'slack',
                          ['active power','reactive power'])
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            self.assertEqual(net.num_vars,(ns+nrg)*self.T)

            x0 = net.get_var_values()
            net.set_var_values(x0+np.random.randn(x0.size))
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('generator powers regularization',1.,net)

            self.assertEqual(func.name,'generator powers regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,(nrg+ns)*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,(nrg+ns)*self.T)
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

            # Value
            func.eval(x0)
            phi = 0
            for t in range(self.T):
                for gen in net.generators:
                    dP = np.max([1e-4,gen.P_max-gen.P_min])
                    dQ = np.max([1e-4,gen.Q_max-gen.Q_min])
                    Pmid = (gen.P_max+gen.P_min)/2.
                    Qmid = (gen.Q_max+gen.Q_min)/2.
                    phi += 0.5*(((gen.P[t]-Pmid)/dP)**2.)
                    phi += 0.5*(((gen.Q[t]-Qmid)/dQ)**2.)
            self.assertLess(np.abs(func.phi-phi),1e-10*(func.phi+1.))
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            self.assertLess(np.abs(func.phi-phi),1e-10*(func.phi+1.))

    def test_func_REG_VANG(self):

        # Constants
        h = 1e-8

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            nb = net.num_buses
            ns = net.get_num_slack_buses()

            # Vars
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])
            self.assertEqual(net.num_vars,2*(nb-ns)*self.T)

            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            net.set_var_values(x0)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('voltage angle regularization',1.,net)

            self.assertEqual(func.name,'voltage angle regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            # Manual count
            man_Hphi_nnz = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                bk = br.bus_k
                bm = br.bus_m
                if not bk.is_slack():
                    man_Hphi_nnz +=1
                if not bm.is_slack():
                    man_Hphi_nnz +=1
                if not bk.is_slack() and not bm.is_slack():
                    man_Hphi_nnz += 1
            man_Hphi_nnz += nb-ns

            func.analyze()
            self.assertEqual(func.Hphi_nnz,man_Hphi_nnz*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,man_Hphi_nnz*self.T)
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

            # Value
            func.eval(net.get_var_values())
            phi = 0
            dw = 3.1416
            for t in range(self.T):
                for bus in net.buses:
                    phi += 0.5*((bus.v_ang[t]/dw)**2.)
                for branch in net.branches:
                    phi += 0.5*(((branch.bus_k.v_ang[t]-branch.bus_m.v_ang[t]-branch.phase[t])/dw)**2.)
            self.assertLess(np.abs(func.phi-phi),1e-10*(phi+1.))
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            self.assertLess(np.abs(func.phi-phi),1e-10*(phi+1.))

    def test_func_REG_RATIO(self):

        # Constants
        h = 1e-8

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          ['tap ratio','tap ratio deviation'])
            self.assertEqual(net.num_vars,3*net.get_num_tap_changers_v()*self.T)

            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            net.set_var_values(x0)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('tap ratio regularization',1.,net)

            self.assertEqual(func.name,'tap ratio regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,3*net.get_num_tap_changers_v()*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,3*net.get_num_tap_changers_v()*self.T)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for tau in range(self.T):
                for i in range(net.num_branches):
                    br = net.get_branch(i)
                    if br.is_tap_changer_v():
                        self.assertTrue(br.has_flags('variable','tap ratio'))
                        t = x0[br.index_ratio[tau]]
                        tmax = br.ratio_max
                        tmin = br.ratio_min
                        t0 = br.ratio[tau]
                        dt = np.maximum(tmax-tmin,1e-4)
                        phi_manual += 0.5*((t-t0)/dt)**2.
                        phi_manual += 0.5*(x0[br.index_ratio_y[tau]]/dt)**2.
                        phi_manual += 0.5*(x0[br.index_ratio_z[tau]]/dt)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-10*(phi_manual+1.))

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

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          ['susceptance','susceptance deviation'])
            self.assertEqual(net.num_vars,3*net.get_num_switched_shunts()*self.T)

            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            net.set_var_values(x0)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('susceptance regularization',1.,net)

            self.assertEqual(func.name,'susceptance regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,3*net.get_num_switched_shunts()*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,3*net.get_num_switched_shunts()*self.T)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for t in range(self.T):
                for i in range(net.num_shunts):
                    sh = net.get_shunt(i)
                    if sh.is_switched_v():
                        self.assertTrue(sh.has_flags('variable','susceptance'))
                        b = x0[sh.index_b[t]]
                        bmax = sh.b_max
                        bmin = sh.b_min
                        b0 = sh.b[t]
                        db = np.maximum(bmax-bmin,1e-4)
                        phi_manual += 0.5*((b-b0)/db)**2.
                        phi_manual += 0.5*(x0[sh.index_y[t]]/db)**2.
                        phi_manual += 0.5*(x0[sh.index_z[t]]/db)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-10*(phi_manual+1.))

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

        # Single period
        h = 1e-9
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # Vars
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            self.assertEqual(net.num_vars,2*net.get_num_generators())
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('generation cost',1.,net)
            
            self.assertEqual(func.name,'generation cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.get_num_generators())
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,net.get_num_generators())
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

            # value check
            val = 0
            for gen in net.generators:
                self.assertTrue(gen.has_flags('variable','active power'))
                val += (gen.cost_coeff_Q0 +
                        gen.cost_coeff_Q1*gen.P +
                        gen.cost_coeff_Q2*(gen.P**2.))
            self.assertLess(np.abs(val-f),1e-8)
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            val = 0
            for gen in net.generators:
                self.assertFalse(gen.has_flags('variable','active power'))
                val += (gen.cost_coeff_Q0 +
                        gen.cost_coeff_Q1*gen.P +
                        gen.cost_coeff_Q2*(gen.P**2.))
            self.assertLess(np.abs(val-func.phi),1e-8)

        # Multi period
        h = 1e-9
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Gen curves
            data = {}
            for gen in net.generators:
                P = np.random.rand(self.T)
                gen.P = P
                data[gen.index] = P
                self.assertEqual(gen.num_periods,self.T)

            # Vars
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            self.assertEqual(net.num_vars,2*net.get_num_generators()*self.T)
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('generation cost',1.,net)

            self.assertEqual(func.name,'generation cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.get_num_generators()*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,net.get_num_generators()*self.T)
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

            # value check
            val = 0
            for t in range(self.T):
                for gen in net.generators:
                    self.assertTrue(gen.has_flags('variable','active power'))
                    self.assertEqual(gen.P[t],data[gen.index][t])
                    self.assertEqual(x0[gen.index_P[t]],gen.P[t])
                    val += (gen.cost_coeff_Q0 +
                            gen.cost_coeff_Q1*gen.P[t] +
                            gen.cost_coeff_Q2*(gen.P[t]**2.))
            self.assertLess(np.abs(val-f),1e-10*np.abs(f))
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            val = 0
            for t in range(self.T):
                for gen in net.generators:
                    self.assertFalse(gen.has_flags('variable','active power'))
                    val += (gen.cost_coeff_Q0 +
                            gen.cost_coeff_Q1*gen.P[t] +
                            gen.cost_coeff_Q2*(gen.P[t]**2.))
            self.assertLess(np.abs(val-func.phi),1e-10*np.abs(func.phi))

    def test_func_SP_CONTROLS(self):

        # Constants
        h = 1e-9

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_vars,
                             (2*net.num_buses +
                              net.num_generators +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()))

            # Sparse
            net.set_flags('bus',
                          'sparse',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'sparse',
                          'any',
                          'active power')
            net.set_flags('branch',
                          'sparse',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'sparse',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'sparse',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_sparse,
                             (2*net.num_buses +
                              net.num_generators +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()))

            x0 = net.get_var_values() + np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('sparse controls penalty',1.,net)

            self.assertEqual(func.name,'sparse controls penalty')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            Hphi_nnz_manual = (net.get_num_buses_reg_by_gen() +
                               net.get_num_generators() +
                               net.get_num_tap_changers() +
                               net.get_num_phase_shifters() +
                               net.get_num_switched_shunts())

            func.analyze()
            self.assertEqual(func.Hphi_nnz,Hphi_nnz_manual)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,Hphi_nnz_manual)

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
            self.assertEqual(H.nnz,Hphi_nnz_manual)
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
                    self.assertTrue(branch.has_flags('variable','tap ratio'))
                    val = x0[branch.index_ratio]
                    val0 = branch.ratio
                    dval = np.maximum(branch.ratio_max-branch.ratio_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags('variable','phase shift'))
                    val = x0[branch.index_phase]
                    val0 = branch.phase
                    dval = np.maximum(branch.phase_max-branch.phase_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags('variable','voltage magnitude'))
                    val = x0[bus.index_v_mag]
                    val0 = bus.v_set
                    dval = np.maximum(bus.v_max-bus.v_min,ceps)
                    f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for gen in net.generators:
                self.assertTrue(gen.has_flags('variable','active power'))
                val = x0[gen.index_P]
                val0 = gen.P
                dval = np.maximum(gen.P_max-gen.P_min,ceps)
                f_manual += np.sqrt(((val-val0)/dval)**2. + eps)
            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags('variable','susceptance'))
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

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            self.assertEqual(net.num_vars,2*net.num_buses*self.T)

            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            net.set_var_values(x0)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('soft voltage magnitude limits',1.,net)

            self.assertEqual(func.name,'soft voltage magnitude limits')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.num_buses*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,net.num_buses*self.T)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # f manual
            f_manual = 0
            eps = 1e-4
            for t in range(self.T):
                for bus in net.buses:
                    self.assertTrue(bus.has_flags('variable','voltage magnitude'))
                    dv = np.maximum(bus.v_max-bus.v_min,eps)
                    vmid = 0.5*(bus.v_max+bus.v_min)
                    f_manual += 0.5*(((x0[bus.index_v_mag[t]]-vmid)/dv)**2.)
            self.assertLess(np.abs(f-f_manual),1e-10*(f_manual+1.))

            # f manual 1
            f_manual = 0
            eps = 1e-4
            for t in range(self.T):
                for bus in net.buses:
                    self.assertTrue(bus.has_flags('variable','voltage magnitude'))
                    dv = np.maximum(bus.v_max-bus.v_min,eps)
                    vmid = 0.5*(bus.v_max+bus.v_min)
                    f_manual += 0.5*(((bus.v_mag[t]-vmid)/dv)**2.)
            self.assertLess(np.abs(f-f_manual),1e-10*(f_manual+1.))

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

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            self.assertEqual(net.num_vars,net.get_num_phase_shifters()*self.T)

            # values
            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            net.set_var_values(x0)
            x0 = net.get_var_values()+np.random.randn(net.num_vars)
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('phase shift regularization',1.,net)

            self.assertEqual(func.name,'phase shift regularization')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.get_num_phase_shifters()*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,net.get_num_phase_shifters()*self.T)
            self.assertTrue(np.all(H.row == H.col))

            self.assertTrue(not np.any(np.isinf(g)))
            self.assertTrue(not np.any(np.isnan(g)))
            self.assertTrue(not np.any(np.isinf(H.data)))
            self.assertTrue(not np.any(np.isnan(H.data)))

            # manual phi check
            phi_manual = 0;
            for t in range(self.T):
                for i in range(net.num_branches):
                    br = net.get_branch(i)
                    if br.is_phase_shifter():
                        self.assertTrue(br.has_flags('variable','phase shift'))
                        p = x0[br.index_phase[t]]
                        pmax = br.phase_max
                        pmin = br.phase_min
                        p0 = br.phase[t]
                        dp = np.maximum(pmax-pmin,1e-4)
                        phi_manual += 0.5*((p-p0)/dp)**2.
            self.assertLess(np.abs(func.phi-phi_manual),1e-10*(phi_manual+1.))

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

        # Single period
        h = 1e-9
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # Vars
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            self.assertEqual(net.num_vars,net.num_loads)
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('consumption utility',1.,net)

            self.assertEqual(func.name,'consumption utility')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.num_loads)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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

            # value check
            val = 0
            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                val += (load.util_coeff_Q0 +
                        load.util_coeff_Q1*load.P +
                        load.util_coeff_Q2*(load.P**2.))
            self.assertLess(np.abs(val-f),1e-7)
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            val = 0
            for load in net.loads:
                self.assertFalse(load.has_flags('variable','active power'))
                val += (load.util_coeff_Q0 +
                        load.util_coeff_Q1*load.P +
                        load.util_coeff_Q2*(load.P**2.))
            self.assertLess(np.abs(val-func.phi),1e-7)

        # Multi period
        h = 1e-9
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            self.assertEqual(net.num_vars,net.num_loads*self.T)
            self.assertGreater(net.num_vars,0)

            # Load curves
            data = {}
            for load in net.loads:
                P = np.random.rand(self.T)
                load.P = P
                data[load.index] = P
                self.assertEqual(load.num_periods,self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('consumption utility',1.,net)

            self.assertEqual(func.name,'consumption utility')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,net.num_loads*self.T)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
            self.assertEqual(H.nnz,net.num_loads*self.T)
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

            # value check
            val = 0
            for t in range(self.T):
                for load in net.loads:
                    self.assertTrue(load.has_flags('variable','active power'))
                    self.assertEqual(load.P[t],data[load.index][t])
                    self.assertEqual(x0[load.index_P[t]],load.P[t])
                    val += (load.util_coeff_Q0 +
                            load.util_coeff_Q1*load.P[t] +
                            load.util_coeff_Q2*(load.P[t]**2.))
            self.assertLess(np.abs(val-f),1e-10*np.abs(f))
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            func.analyze()
            func.eval(np.zeros(0))
            val = 0
            for t in range(self.T):
                for load in net.loads:
                    self.assertFalse(load.has_flags('variable','active power'))
                    val += (load.util_coeff_Q0 +
                            load.util_coeff_Q1*load.P[t] +
                            load.util_coeff_Q2*(load.P[t]**2.))
            self.assertLess(np.abs(val-func.phi),1e-10*np.abs(func.phi))

    def test_func_NETCON_COST(self):

        # Single period
        h = 1e-7
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # prices
            for bus in net.buses:
                bus.price = (bus.index%10)*0.5123

            # vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = (vargen.index%10)*0.3233+0.1

            # batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P *= -1.

            # Vars
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('battery',
                          'variable',
                          'any',
                          'charging power')
            self.assertEqual(net.num_vars,
                             (net.num_loads+
                              net.num_generators+
                              net.num_var_generators+
                              2*net.num_batteries))
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('net consumption cost',1.,net)

            self.assertEqual(func.name,'net consumption cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,0)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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
                    self.assertTrue(load.has_flags('variable','active power'))
                    val += bus.price*load.P
                for bat in bus.batteries:
                    self.assertTrue(bat.has_flags('variable','charging power'))
                    val += bus.price*bat.P
                for gen in bus.generators:
                    self.assertTrue(gen.has_flags('variable','active power'))
                    val -= bus.price*gen.P
                for vargen in bus.var_generators:
                    self.assertTrue(vargen.has_flags('variable','active power'))
                    val -= bus.price*vargen.P
            self.assertLess(np.abs(val-f),1e-10*np.abs(f))

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

            func = pf.Function('net consumption cost',1.,net)
            self.assertEqual(func.name,'net consumption cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

            x0 = net.get_var_values()

            func.analyze()
            func.eval(x0)

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

            # value check
            val = 0
            for bus in net.buses:
                for load in bus.loads:
                    self.assertFalse(load.has_flags('variable','active power'))
                    val += bus.price*load.P
                for bat in bus.batteries:
                    self.assertFalse(bat.has_flags('variable','charging power'))
                    val += bus.price*bat.P
                for gen in bus.generators:
                    self.assertFalse(gen.has_flags('variable','active power'))
                    val -= bus.price*gen.P
                for vargen in bus.var_generators:
                    self.assertFalse(vargen.has_flags('variable','active power'))
                    val -= bus.price*vargen.P
            self.assertLess(np.abs(val-func.phi),1e-10*np.abs(f))

        # Multi period
        h = 1e-7
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # gens
            for gen in net.generators:
                self.assertEqual(gen.num_periods,self.T)
                gen.P = np.random.rand(self.T)*10.

            # loads
            for load in net.loads:
                self.assertEqual(load.num_periods,self.T)
                load.P = np.random.rand(self.T)*10.

            # prices
            for bus in net.buses:
                self.assertEqual(bus.num_periods,self.T)
                bus.price = np.random.rand(self.T)*10.

            # vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                self.assertEqual(vargen.num_periods,self.T)
                vargen.P = np.random.randn(self.T)*10.

            # batteries
            for bat in net.batteries:
                self.assertEqual(bat.num_periods,self.T)
                bat.P = np.random.randn(self.T)*10.

            # Vars
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('battery',
                          'variable',
                          'any',
                          'charging power')
            self.assertEqual(net.num_vars,
                             (net.num_loads+
                              net.num_generators+
                              net.num_var_generators+
                              2*net.num_batteries)*self.T)
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Function
            func = pf.Function('net consumption cost',1.,net)

            self.assertEqual(func.name,'net consumption cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

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

            self.assertEqual(func.Hphi_nnz,0)

            func.analyze()
            self.assertEqual(func.Hphi_nnz,0)
            func.eval(x0)
            self.assertEqual(func.Hphi_nnz,0)

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

            # value and grad check
            val = 0
            for t in range(self.T):
                for bus in net.buses:
                    for load in bus.loads:
                        self.assertTrue(load.has_flags('variable','active power'))
                        self.assertEqual(load.P[t],x0[load.index_P[t]])
                        self.assertEqual(g[load.index_P[t]],bus.price[t])
                        val += bus.price[t]*load.P[t]
                    for bat in bus.batteries:
                        self.assertTrue(bat.has_flags('variable','charging power'))
                        self.assertEqual(bat.P[t],x0[bat.index_Pc[t]]-x0[bat.index_Pd[t]])
                        self.assertEqual(g[bat.index_Pc[t]],bus.price[t])
                        self.assertEqual(g[bat.index_Pd[t]],-bus.price[t])
                        val += bus.price[t]*bat.P[t]
                    for gen in bus.generators:
                        self.assertTrue(gen.has_flags('variable','active power'))
                        self.assertEqual(gen.P[t],x0[gen.index_P[t]])
                        self.assertEqual(g[gen.index_P[t]],-bus.price[t])
                        val -= bus.price[t]*gen.P[t]
                    for vargen in bus.var_generators:
                        self.assertTrue(vargen.has_flags('variable','active power'))
                        self.assertEqual(vargen.P[t],x0[vargen.index_P[t]])
                        self.assertEqual(g[vargen.index_P[t]],-bus.price[t])
                        val -= bus.price[t]*vargen.P[t]
            self.assertLess(np.abs(val-f),1e-10*np.abs(f))

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

            func = pf.Function('net consumption cost',1.,net)
            self.assertEqual(func.name,'net consumption cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

            x0 = net.get_var_values()

            func.analyze()
            func.eval(x0)

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))

            # value check
            val = 0
            for t in range(self.T):
                for bus in net.buses:
                    for load in bus.loads:
                        self.assertFalse(load.has_flags('variable','active power'))
                        val += bus.price[t]*load.P[t]
                    for bat in bus.batteries:
                        self.assertFalse(bat.has_flags('variable','charging power'))
                        val += bus.price[t]*bat.P[t]
                    for gen in bus.generators:
                        self.assertFalse(gen.has_flags('variable','active power'))
                        val -= bus.price[t]*gen.P[t]
                    for vargen in bus.var_generators:
                        self.assertFalse(vargen.has_flags('variable','active power'))
                        val -= bus.price[t]*vargen.P[t]
            self.assertLess(np.abs(val-func.phi),1e-10*np.abs(f))

    def test_robustness(self):

        for case in test_cases.CASES:

            net = pf.Network(self.T) # multiperiod
            
            functions = [pf.Function('generation cost',1.,net),
                         pf.Function('phase shift regularization',1.,net),
                         pf.Function('generator powers regularization',1.,net),
                         pf.Function('tap ratio regularization',1.,net),
                         pf.Function('susceptance regularization',1.,net),
                         pf.Function('voltage angle regularization',1.,net),
                         pf.Function('voltage magnitude regularization',1.,net),
                         pf.Function('soft voltage magnitude limits',1.,net),
                         pf.Function('sparse controls penalty',1.,net),
                         pf.Function('consumption utility',1.,net),
                         pf.Function('net consumption cost',1.,net)]

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
            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Gen powers
            for gen in net.generators:
                gen.P = gen.P + 0.01

            functions = [pf.Function('generation cost',1.,net),
                         pf.Function('phase shift regularization',1.,net),
                         pf.Function('generator powers regularization',1.,net),
                         pf.Function('tap ratio regularization',1.,net),
                         pf.Function('susceptance regularization',1.,net),
                         pf.Function('voltage angle regularization',1.,net),
                         pf.Function('voltage magnitude regularization',1.,net),
                         pf.Function('soft voltage magnitude limits',1.,net),
                         pf.Function('sparse controls penalty',1.,net),
                         pf.Function('consumption utility',1.,net),
                         pf.Function('net consumption cost',1.,net)]

            # Update network
            list(map(lambda f: f.update_network(),functions))

            # After updating network
            list(map(lambda f: f.analyze(),functions))
            list(map(lambda f: f.eval(x0),functions))

            for f in functions:
                self.assertTrue(isinstance(f.gphi,np.ndarray))
                self.assertTrue(isinstance(f.Hphi,coo_matrix))
                self.assertEqual(f.gphi.size,0)
                self.assertEqual(f.Hphi.nnz,0)

            # Add variables
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_vars,
                             (2*net.num_buses +
                              2*net.num_generators +
                              net.num_loads +
                              net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts())*self.T)

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
            
    def test_func_DUMMY(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            
            # Vars
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            
            self.assertEqual(net.num_vars,2*net.get_num_generators()*net.num_periods)
            self.assertGreater(net.num_vars,0)

            x0 = net.get_var_values()+np.random.randn(net.num_vars)*1e-2
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Regular function
            funcREF = pf.Function('generation cost',1.2,net)

            # Custom function written in Python
            func = pf.functions.DummyGenCost(0.3,net) 
            
            self.assertEqual(funcREF.weight,1.2)
            self.assertEqual(func.weight,0.3)

            self.assertEqual(funcREF.name,'generation cost')
            self.assertEqual(func.name,'dummy generation cost')

            self.assertTupleEqual(func.gphi.shape,(0,))
            self.assertTupleEqual(func.Hphi.shape,(0,0))
            self.assertTupleEqual(funcREF.gphi.shape,(0,))
            self.assertTupleEqual(funcREF.Hphi.shape,(0,0))

            self.assertTrue(np.all(func.bus_counted == np.zeros(net.num_buses*net.num_periods,np.dtype('bool'))))
            func.bus_counted[2] = True
            self.assertTrue(func.bus_counted[2])

            self.assertEqual(func.phi,0.)
            self.assertEqual(func.gphi.size,0)
            self.assertEqual(func.Hphi.nnz,0)
           
            funcREF.analyze()
            func.analyze()
            
            self.assertEqual(func.phi,funcREF.phi)
            self.assertEqual(func.phi,0.)
            self.assertTrue(np.all(func.gphi == funcREF.gphi))
            self.assertEqual(func.gphi.size,net.num_vars)
            self.assertEqual(func.Hphi.nnz,net.get_num_generators()*net.num_periods)
            self.assertTrue(np.all(func.Hphi.row == funcREF.Hphi.row))
            self.assertTrue(np.all(func.Hphi.col == funcREF.Hphi.col))
            self.assertTrue(np.all(func.Hphi.data == funcREF.Hphi.data))
            
            funcREF.eval(x0)
            func.eval(x0)
            
            self.assertNotEqual(func.phi,0.)
            self.assertLess(abs(func.phi-funcREF.phi),1e-8)
            self.assertTrue(np.all(func.gphi == funcREF.gphi))

            phi = 0
            for t in range(net.num_periods):
                for gen in net.generators:
                    P = x0[gen.index_P[t]]
                    phi += gen.cost_coeff_Q0+gen.cost_coeff_Q1*P+gen.cost_coeff_Q2*P*P
            self.assertLess(abs(func.phi-phi),1e-8)

    def tearDown(self):

        pass
