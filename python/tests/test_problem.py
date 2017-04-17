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
from numpy.linalg import norm
from scipy.sparse import coo_matrix,triu,bmat

NUM_TRIALS = 25
EPS = 3. # %
TOL = 1e-4

class TestProblem(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.p = pf.Problem()

        # Random
        np.random.seed(0)

    def test_problem_LSNR(self):

        # Constants
        h = 1e-9
        
        p = self.p

        for case in test_cases.CASES:
            
            p.clear()
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            p.network = net
            
            # Variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
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
                             2*(net.num_buses-net.get_num_slack_buses()) +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers_v() + 
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
                             
            # Fixed
            net.set_flags('bus',
                          'fixed',
                          'regulated by generator',
                          'voltage magnitude')
            net.set_flags('branch',
                          'fixed',
                          'tap changer - v',
                          'tap ratio')
            net.set_flags('branch',
                          'fixed',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'fixed',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_fixed,
                             net.get_num_buses_reg_by_gen() +
                             net.get_num_tap_changers_v() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts())
                             
            # Constraints
            p.add_constraint(pf.Constraint('AC power balance',net))
            p.add_constraint(pf.Constraint('generator active power participation',net))
            p.add_constraint(pf.Constraint('generator reactive power participation',net))
            p.add_constraint(pf.Constraint('variable fixing',net))
            self.assertEqual(len(p.constraints),4)

            # Check adding redundant constraints
            p.add_constraint(pf.Constraint('generator active power participation',net))
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

            # Numbers
            self.assertEqual(x0.size,p.num_primal_variables)
            self.assertEqual(A.shape[0],p.num_linear_equality_constraints)
            self.assertEqual(f.size,p.num_nonlinear_equality_constraints)
            self.assertEqual(p.num_primal_variables,p.get_num_primal_variables())
            self.assertEqual(p.num_linear_equality_constraints,p.get_num_linear_equality_constraints())
            self.assertEqual(p.num_nonlinear_equality_constraints,p.get_num_nonlinear_equality_constraints())
 
            # phi
            self.assertTrue(type(phi) is float)
            self.assertEqual(phi,0.)
            
            # gphi
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(net.num_vars,))
            self.assertLess(norm(gphi),1e-10)

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
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)
                
            # Check Hcombined
            coeff = np.random.randn(f.shape[0])
            p.eval(x0)
            self.assertRaises(pf.ProblemError,p.combine_H,np.zeros(f.shape[0]+1),False)
            self.assertTrue(p.has_error())
            p.clear_error()
            self.assertFalse(p.has_error())
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
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
                if c.name == 'AC power balance':
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

        for case in test_cases.CASES:
            
            p.clear()
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            p.network = net
            
            # Variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('bus',
                          'variable',
                          ['not slack','regulated by generator'],
                          'voltage magnitude deviation')
            net.set_flags('bus',
                          'variable',
                          'regulated by transformer',
                          'voltage magnitude violation')
            net.set_flags('bus',
                          'variable',
                          'regulated by shunt',
                          'voltage magnitude violation')
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          ['tap ratio','tap ratio deviation'])
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
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
            p.add_constraint(pf.Constraint('AC power balance',net))
            p.add_constraint(pf.Constraint('generator active power participation',net))
            p.add_constraint(pf.Constraint('generator reactive power participation',net))
            p.add_constraint(pf.Constraint('voltage regulation by generators',net))
            p.add_constraint(pf.Constraint('voltage regulation by transformers',net))
            p.add_constraint(pf.Constraint('voltage regulation by shunts',net))
            self.assertEqual(len(p.constraints),6)

            # Check adding redundant constraints
            p.add_constraint(pf.Constraint('AC power balance',net))
            self.assertEqual(len(p.constraints),6)
            
            # Functions
            p.add_function(pf.Function('voltage magnitude regularization',1.,net))
            p.add_function(pf.Function('voltage angle regularization',5.,net))
            p.add_function(pf.Function('generator powers regularization',8.,net))
            p.add_function(pf.Function('tap ratio regularization',3.,net))
            p.add_function(pf.Function('susceptance regularization',1.,net))
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

            # Numbers
            self.assertEqual(x0.size,p.num_primal_variables)
            self.assertEqual(A.shape[0],p.num_linear_equality_constraints)
            self.assertEqual(f.size,p.num_nonlinear_equality_constraints)
                        
            # phi
            self.assertTrue(type(phi) is float)
            self.assertGreater(phi,0.)
            man_phi = sum(f.weight*f.phi for f in p.functions)
            self.assertLess(np.abs(man_phi-phi),1e-10)

            # gphi
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(net.num_vars,))
            man_gphi = sum(f.weight*f.gphi for f in p.functions)
            self.assertLess(norm(man_gphi-gphi),1e-10)

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
                error = 100.*norm(gd_exact-gd_approx)/np.maximum(norm(gd_exact),TOL)
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
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Check Hcombined
            coeff = np.random.randn(f.shape[0])
            p.eval(x0)
            self.assertRaises(pf.ProblemError,p.combine_H,np.zeros(f.shape[0]+1),False)
            self.assertTrue(p.has_error())
            p.clear_error()
            self.assertFalse(p.has_error())
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
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
                if c.name == 'AC power balance':
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

        for case in test_cases.CASES:
            
            p.clear()
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            p.network = net

            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')
            self.assertEqual(net.num_vars,net.num_buses)

            l = p.get_lower_limits()
            u = p.get_upper_limits()
            self.assertTrue(isinstance(l,np.ndarray))
            self.assertTrue(isinstance(u,np.ndarray))
            self.assertTupleEqual(l.shape,(net.num_buses,))
            self.assertTupleEqual(u.shape,(net.num_buses,))
            for bus in net.buses:
                self.assertEqual(bus.v_max_reg,u[bus.index_v_mag])
                self.assertEqual(bus.v_min_reg,l[bus.index_v_mag])

    def test_problem_Glu_construction(self):

        p = self.p

        for case in test_cases.CASES:
            
            p.clear()
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            p.network = net

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_bounded,0)
            
            # flags
            net.set_flags('bus',
                          ['variable','bounded'],
                          'any',
                          ['voltage magnitude','voltage angle'])
            
            self.assertGreater(net.num_buses,0)
            self.assertEqual(net.num_vars,net.num_buses*2)
            self.assertEqual(net.num_bounded,net.num_buses*2)
            
            self.assertEqual(len(p.constraints),0)

            p.add_constraint(pf.Constraint('variable bounds',net))
            p.add_constraint(pf.Constraint('DC branch flow limits',net))

            self.assertEqual(len(p.constraints),2)

            constr1 = p.find_constraint('variable bounds')
            constr2 = p.find_constraint('DC branch flow limits')
            self.assertRaises(pf.ProblemError,p.find_constraint,'AC power balance')

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

            num_dc = len([br for br in net.branches if br.ratingA != 0.])

            self.assertTupleEqual(l1.shape,(net.num_vars,))
            self.assertTupleEqual(u1.shape,(net.num_vars,))
            self.assertTupleEqual(G1.shape,(net.num_vars,net.num_vars))
            self.assertTupleEqual(l2.shape,(num_dc,))
            self.assertTupleEqual(u2.shape,(num_dc,))
            self.assertTupleEqual(G2.shape,(num_dc,net.num_vars))
            self.assertTupleEqual(l.shape,(net.num_vars+num_dc,))
            self.assertTupleEqual(u.shape,(net.num_vars+num_dc,))
            self.assertTupleEqual(G.shape,(net.num_vars+num_dc,net.num_vars))
            
            self.assertLess(norm(l-np.hstack((l1,l2)),np.inf),1e-12)
            self.assertLess(norm(u-np.hstack((u1,u2)),np.inf),1e-12)

            self.assertGreater(G.nnz,0)
            self.assertGreater(bmat([[G1],[G2]],format='coo').nnz,0)
            E = G - bmat([[G1],[G2]])
            self.assertEqual(E.nnz,0)

    def test_problem_ACOPF_with_thermal(self):

        p = self.p

        for case in test_cases.CASES:
            
            p.clear()
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            p.network = net

            for branch in net.branches:
                if branch.ratingA == 0.:
                    branch.ratingA = 100.
            
            # Variables
            net.set_flags('bus',
                          ['variable','bounded'],
                          'any',
                          'voltage magnitude')
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('generator',
                          ['variable','bounded'],
                          'adjustable active power',
                          'active power')
            net.set_flags('generator',
                          ['variable','bounded'],
                          'regulator',
                          'reactive power')
            net.set_flags('branch',
                          ['variable','bounded'],
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          ['variable','bounded'],
                          'phase shifter',
                          'phase shift')
           
            self.assertEqual(net.num_vars,(2*net.num_buses-net.get_num_slack_buses() +
                                           net.get_num_P_adjust_gens() + 
                                           net.get_num_reg_gens()+
                                           net.get_num_tap_changers()+
                                           net.get_num_phase_shifters()))
            self.assertEqual(net.num_bounded,(net.num_buses +
                                              net.get_num_P_adjust_gens() + 
                                              net.get_num_reg_gens()+
                                              net.get_num_tap_changers()+
                                              net.get_num_phase_shifters()))

            p.add_constraint(pf.Constraint('AC power balance',net))
            p.add_constraint(pf.Constraint('AC branch flow limits',net))
            p.add_constraint(pf.Constraint('variable bounds',net))
            p.add_function(pf.Function('generation cost',1.,net))
            p.analyze()

            # Extra vars
            self.assertGreater(p.num_extra_vars,0)
            self.assertEqual(p.num_extra_vars,net.num_branches*2)
            
            # Init point
            x0 = p.get_init_point()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars+p.num_extra_vars,))
            x = net.get_var_values()
            x0_check = np.hstack((x,np.zeros(p.num_extra_vars)))
            self.assertTrue(np.all(x0 == x0_check))
            
            p.analyze()
            r = np.random.randn(p.num_extra_vars)
            x0[net.num_vars:] = r
            p.eval(x0)
            
            phi = p.phi
            gphi = p.gphi.copy()
            Hphi = p.Hphi.copy()

            f = p.f.copy()
            b = p.b.copy()
            A = p.A.copy()
            J = p.J.copy()
            G = p.G.copy()
            l = p.l.copy()
            u = p.u.copy()

            # Numbers
            self.assertEqual(x0.size,p.num_primal_variables)
            self.assertEqual(A.shape[0],p.num_linear_equality_constraints)
            self.assertEqual(f.size,p.num_nonlinear_equality_constraints)
                        
            # phi
            self.assertTrue(type(phi) is float)
            self.assertGreaterEqual(phi,0.)
            man_phi = sum(f.weight*f.phi for f in p.functions)
            self.assertLess(np.abs(man_phi-phi),1e-10)

            # gphi
            self.assertTrue(type(gphi) is np.ndarray)
            self.assertTupleEqual(gphi.shape,(net.num_vars+p.num_extra_vars,))
            man_gphi = sum(f.weight*f.gphi for f in p.functions)
            self.assertLess(norm(man_gphi-gphi[:net.num_vars]),1e-10)
            self.assertTrue(np.all(gphi[net.num_vars:] == 0.))

            # Hphi
            self.assertTrue(type(Hphi) is coo_matrix)
            self.assertTupleEqual(Hphi.shape,(net.num_vars+p.num_extra_vars,net.num_vars+p.num_extra_vars))
            self.assertGreater(Hphi.nnz,0)
    
            # f
            self.assertTrue(type(f) is np.ndarray)
            f_size = sum(c.f.shape[0] for c in p.constraints)
            f_man = np.zeros(0)
            offset = 0
            for c in p.constraints:
                if c.Jbar.shape[0]:
                    f_man = np.hstack((f_man,c.f+c.Jbar*r[offset:offset+c.num_extra_vars]))
                else:
                    f_man = np.hstack((f_man,c.f))
                offset += c.num_extra_vars
            self.assertTupleEqual(f.shape,(f_size,))
            self.assertEqual(f.size,f_man.size)
            self.assertTrue(np.all(f_man == f))
            
            # b
            self.assertTrue(type(b) is np.ndarray)
            b_size = sum(c.b.shape[0] for c in p.constraints)
            self.assertTupleEqual(b.shape,(b_size,))

            # J
            self.assertTrue(type(J) is coo_matrix)
            J_size = sum(c.J.shape[0] for c in p.constraints)
            J_nnz = sum(c.J.nnz+c.Jbar.nnz for c in p.constraints)
            J_man = []
            for c in p.constraints:
                if c.num_extra_vars > 0:
                    J_man.append([c.J,c.Jbar])
                else:
                    J_man.append([c.J,None])
            J_man = bmat(J_man,format='coo')
            self.assertTupleEqual(J.shape,(J_size,net.num_vars+p.num_extra_vars))
            self.assertEqual(J.nnz,J_nnz)
            self.assertTupleEqual(J_man.shape,J.shape)
            self.assertLess(norm((J_man-J).data),1e-10)

            # G, l, u
            self.assertTrue(type(G) is coo_matrix)
            G_size = sum(c.G.shape[0] for c in p.constraints)
            G_nnz = sum(c.G.nnz+c.Gbar.nnz for c in p.constraints)
            G_man = []
            for c in p.constraints:
                if c.num_extra_vars > 0:
                    G_man.append([c.G,c.Gbar])
                else:
                    G_man.append([c.G,None])
            G_man = bmat(G_man,format='coo')
            self.assertTupleEqual(G.shape,(G_size,net.num_vars+p.num_extra_vars))
            self.assertEqual(G.nnz,G_nnz)
            self.assertEqual(l.size,G_size)
            self.assertEqual(u.size,G_size)
            self.assertTupleEqual(G_man.shape,G.shape)
            self.assertLess(norm((G_man-G).data),1e-10)

            # A
            self.assertTrue(type(A) is coo_matrix)
            A_size = sum(c.A.shape[0] for c in p.constraints)
            A_nnz = sum(c.A.nnz for c in p.constraints)
            self.assertTupleEqual(A.shape,(A_size,net.num_vars+p.num_extra_vars))
            self.assertEqual(A.nnz,A_nnz)

            # Check gphi
            h = 1e-9
            phi0 = phi
            gphi0 = gphi.copy()
            self.assertTrue(np.all(gphi0[net.num_vars:] == 0.))
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars+p.num_extra_vars)
    
                x = x0 + h*d
                
                p.eval(x)
                phi1 = p.phi
                
                gd_exact = np.dot(gphi0,d)
                gd_approx = (phi1-phi0)/h
                error = 100.*norm(gd_exact-gd_approx)/np.maximum(norm(gd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Check J
            h = 1e-12
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars+p.num_extra_vars)
    
                x = x0 + h*d
                
                p.eval(x)
                f1 = p.f
                
                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)
                
            # Check Hphi
            h = 1e-9
            gphi0 = gphi.copy()
            Hphi0 = Hphi.copy()
            Hphi0 = Hphi0 + Hphi0.T - triu(Hphi0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars+p.num_extra_vars)
    
                x = x0 + h*d
                
                p.eval(x)

                gphi1 = p.gphi.copy()
                
                Hd_exact = Hphi0*d
                Hd_approx = (gphi1-gphi0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Check Hcombined
            h = 1e-13
            coeff = np.random.randn(f.shape[0])
            p.eval(x0)
            self.assertRaises(pf.ProblemError,p.combine_H,np.zeros(f.shape[0]+1),False)
            self.assertTrue(p.has_error())
            p.clear_error()
            self.assertFalse(p.has_error())
            p.combine_H(coeff,False)
            J0 = p.J.copy()
            g0 = J0.T*coeff
            H0 = p.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars+p.num_extra_vars,net.num_vars+p.num_extra_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):
                
                d = np.random.randn(net.num_vars+p.num_extra_vars)
                
                x = x0 + h*d
                
                p.eval(x)
                
                g1 = p.J.T*coeff
                
                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

    def test_problem_with_DUMMY_func(self):
        
        for case in test_cases.CASES:
            
            p1 = pf.Problem()
            p2 = pf.Problem()

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            p1.network = net
            p2.network = net
            
            # Variables
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')

            self.assertEqual(net.num_vars,net.num_generators)
            self.assertEqual(p1.get_num_primal_variables(),net.num_generators)
            self.assertEqual(p2.get_num_primal_variables(),net.num_generators)

            p1.add_function(pf.Function('generation cost',1.,net))
            p2.add_function(pf.functions.DummyGenCost(1,net))

            self.assertEqual(len(p1.functions),1)
            self.assertEqual(len(p2.functions),1)

            self.assertEqual(p1.functions[0].name,"generation cost")
            self.assertEqual(p2.functions[0].name,"dummy generation cost")

            self.assertTupleEqual(p1.functions[0].Hphi.shape,(0,0))
            self.assertTupleEqual(p2.functions[0].Hphi.shape,(0,0))
 
            p1.analyze()
            p2.analyze()

            self.assertEqual(p1.phi,0.)
            self.assertEqual(p2.phi,0.)
            self.assertEqual(p1.gphi.size,p2.gphi.size)
            self.assertTrue(np.all(p1.gphi == p2.gphi))
            self.assertEqual(p1.Hphi.nnz,p2.Hphi.nnz)
            self.assertTrue(np.all(p1.Hphi.row == p2.Hphi.row))
            self.assertTrue(np.all(p1.Hphi.col == p2.Hphi.col))
            self.assertTrue(np.all(p1.Hphi.data == p2.Hphi.data))
                        
            p1.eval(net.get_var_values())
            p2.eval(net.get_var_values())
            
            self.assertGreaterEqual(p1.phi,0)
            self.assertLess(abs(p1.phi-p2.phi),1e-8)
            self.assertTrue(np.all(p1.gphi == p2.gphi))
            self.assertEqual(p1.Hphi.nnz,p2.Hphi.nnz)
            self.assertTrue(np.all(p1.Hphi.row == p2.Hphi.row))
            self.assertTrue(np.all(p1.Hphi.col == p2.Hphi.col))
            self.assertTrue(np.all(p1.Hphi.data == p2.Hphi.data))

    def test_problem_with_DUMMY_constr(self):
        
        for case in test_cases.CASES:
            
            p1 = pf.Problem()
            p2 = pf.Problem()

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            p1.network = net
            p2.network = net

            # Variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            net.set_flags('battery',
                          'variable',
                          'any',
                          'charging power')
            self.assertEqual(net.num_vars,
                             (net.num_buses-net.get_num_slack_buses() +
                              net.num_generators +
                              net.num_loads +
                              net.num_var_generators +
                              net.get_num_phase_shifters()+
                              2*net.num_batteries)*net.num_periods)

            self.assertEqual(p1.get_num_primal_variables(),net.num_vars)
            self.assertEqual(p2.get_num_primal_variables(),net.num_vars)
            
            p1.add_constraint(pf.Constraint('DC power balance',net))
            p2.add_constraint(pf.constraints.DummyDCPF(net))

            self.assertEqual(len(p1.constraints),1)
            self.assertEqual(len(p2.constraints),1)

            self.assertEqual(p1.constraints[0].name,"DC power balance")
            self.assertEqual(p2.constraints[0].name,"dummy DC power balance")

            self.assertEqual(p1.find_constraint("DC power balance").name,"DC power balance")
            self.assertEqual(p2.find_constraint("dummy DC power balance").name,"dummy DC power balance")
            self.assertRaises(pf.ProblemError,p2.find_constraint,"DC power balance")
            
            p1.analyze()
            p2.analyze()

            self.assertTrue(np.all(p1.b == p2.b))
            self.assertTrue(np.all(p1.A.row == p2.A.row))
            self.assertTrue(np.all(p1.A.col == p2.A.col))
            self.assertTrue(np.all(p1.A.data == p2.A.data))

            p1.eval(net.get_var_values())
            p2.eval(net.get_var_values())

            self.assertTrue(np.all(p1.b == p2.b))
            self.assertTrue(np.all(p1.A.row == p2.A.row))
            self.assertTrue(np.all(p1.A.col == p2.A.col))
            self.assertTrue(np.all(p1.A.data == p2.A.data))

    def tearDown(self):
        
        pass
