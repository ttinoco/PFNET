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
from scipy.sparse import coo_matrix,triu,bmat

NUM_TRIALS = 25
EPS = 2e0 # %
TOL = 1e-4

class TestProblem(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()
        self.p = pf.Problem()

        # Random
        np.random.seed(0)

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
            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_NOT_SLACK,
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          pf.GEN_PROP_SLACK,
                          'active power')
            net.set_flags('generator',
                          'variable',
                          pf.GEN_PROP_REG,
                          'reactive power')
            net.set_flags('branch',
                          'variable',
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          'phase shift')
            net.set_flags('shunt',
                          'variable',
                          pf.SHUNT_PROP_SWITCHED_V,
                          'susceptance')
            
            self.assertEqual(net.num_vars,
                             2*(net.num_buses-net.get_num_slack_buses()) +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers_v() + 
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
                             
            # Fixed
            net.set_flags('bus',
                          'fixed',
                          pf.BUS_PROP_REG_BY_GEN,
                          'voltage magnitude')
            net.set_flags('branch',
                          'fixed',
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          'tap ratio')
            net.set_flags('branch',
                          'fixed',
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          'phase shift')
            net.set_flags('shunt',
                          'fixed',
                          pf.SHUNT_PROP_SWITCHED_V,
                          'susceptance')
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
            self.assertTrue(np.all(x0 == p.x))
            
            # Before
            phi = p.phi
            gphi = p.gphi
            Hphi = p.Hphi

            f = p.f
            b = p.b
            A = p.A
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            p.store_sensitivities(np.zeros(p.A.shape[0]),sens,None,None)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,sens[2*bus.index+offset])
                self.assertEqual(bus.sens_Q_balance,sens[2*bus.index+1+offset])
            self.assertRaises(pf.ProblemError,
                              p.store_sensitivities,
                              np.zeros(p.A.shape[0]),
                              np.zeros(p.f.size+5),
                              None,
                              None)

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
            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_NOT_SLACK,
                          ['voltage magnitude','voltage angle'])
            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_NOT_SLACK|pf.BUS_PROP_REG_BY_GEN,
                          'voltage magnitude deviation')
            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_REG_BY_TRAN,
                          'voltage magnitude violation')
            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_REG_BY_SHUNT,
                          'voltage magnitude violation')
            net.set_flags('generator',
                          'variable',
                          pf.GEN_PROP_SLACK,
                          'active power')
            net.set_flags('generator',
                          'variable',
                          pf.GEN_PROP_REG,
                          'reactive power')
            net.set_flags('branch',
                          'variable',
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          ['tap ratio','tap ratio deviation'])
            net.set_flags('shunt',
                          'variable',
                          pf.SHUNT_PROP_SWITCHED_V,
                          ['susceptance','susceptance deviation'])                          

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
            r = np.random.randn(net.num_vars)
            x0 = p.get_init_point()+r
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            self.assertTrue(np.all(x0 == p.x+r))
            
            # Before
            phi = p.phi
            gphi = p.gphi
            Hphi = p.Hphi

            f = p.f
            b = p.b
            A = p.A
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
                error = 100.*np.linalg.norm(gd_exact-gd_approx)/np.maximum(np.linalg.norm(gd_exact),TOL)
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            p.store_sensitivities(np.zeros(p.A.shape[0]),sens,None,None)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.sens_P_balance,sens[2*bus.index+offset])
                self.assertEqual(bus.sens_Q_balance,sens[2*bus.index+1+offset])
            self.assertRaises(pf.ProblemError,
                              p.store_sensitivities,
                              np.zeros(p.A.shape[0]),
                              np.zeros(p.f.size+5),
                              None,
                              None)

    def test_problem_limits(self):

        p = self.p
        net = self.net

        for case in test_cases.CASES:
            
            p.clear()
            net.load(case)
            p.network = net

            net.set_flags('bus',
                          'variable',
                          pf.BUS_PROP_ANY,
                          'voltage magnitude')
            self.assertEqual(net.num_vars,net.num_buses)

            l = p.get_lower_limits()
            u = p.get_upper_limits()
            self.assertTrue(isinstance(l,np.ndarray))
            self.assertTrue(isinstance(u,np.ndarray))
            self.assertTupleEqual(l.shape,(net.num_buses,))
            self.assertTupleEqual(u.shape,(net.num_buses,))
            for bus in net.buses:
                self.assertEqual(bus.v_max,u[bus.index_v_mag])
                self.assertEqual(bus.v_min,l[bus.index_v_mag])

    def test_problem_Glu_construction(self):

        p = self.p
        net = self.net

        for case in test_cases.CASES:
            
            p.clear()
            net.load(case)
            p.network = net

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_bounded,0)
            
            # flags
            net.set_flags('bus',
                          ['variable','bounded'],
                          pf.BUS_PROP_ANY,
                          ['voltage magnitude','voltage angle'])
            
            self.assertGreater(net.num_buses,0)
            self.assertEqual(net.num_vars,net.num_buses*2)
            self.assertEqual(net.num_bounded,net.num_buses*2)
            
            self.assertEqual(len(p.constraints),0)

            p.add_constraint(pf.CONSTR_TYPE_LBOUND)
            p.add_constraint(pf.CONSTR_TYPE_DC_FLOW_LIM)

            self.assertEqual(len(p.constraints),2)

            constr1 = p.find_constraint(pf.CONSTR_TYPE_LBOUND)
            constr2 = p.find_constraint(pf.CONSTR_TYPE_DC_FLOW_LIM)
            self.assertRaises(pf.ProblemError,p.find_constraint,pf.CONSTR_TYPE_PF)

            p.analyze()

            l1 = constr1.l
            u1 = constr1.u
            G1 = constr1.G

            l2 = constr2.l
            u2 = constr2.u
            G2 = constr2.G

            l = p.l
            u = p.u
            G = p.G

            self.assertTupleEqual(l1.shape,(net.num_vars,))
            self.assertTupleEqual(u1.shape,(net.num_vars,))
            self.assertTupleEqual(G1.shape,(net.num_vars,net.num_vars))
            self.assertTupleEqual(l2.shape,(net.num_branches,))
            self.assertTupleEqual(u2.shape,(net.num_branches,))
            self.assertTupleEqual(G2.shape,(net.num_branches,net.num_vars))
            self.assertTupleEqual(l.shape,(net.num_vars+net.num_branches,))
            self.assertTupleEqual(u.shape,(net.num_vars+net.num_branches,))
            self.assertTupleEqual(G.shape,(net.num_vars+net.num_branches,net.num_vars))
            
            self.assertLess(np.linalg.norm(l-np.hstack((l1,l2)),np.inf),1e-12)
            self.assertLess(np.linalg.norm(u-np.hstack((u1,u2)),np.inf),1e-12)

            self.assertGreater(G.nnz,0)
            self.assertGreater(bmat([[G1],[G2]],format='coo').nnz,0)
            E = G - bmat([[G1],[G2]])
            self.assertEqual(E.nnz,0)
            
    def tearDown(self):
        
        pass
