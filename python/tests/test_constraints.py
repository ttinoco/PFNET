#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import os
import pfnet as pf
import unittest
from . import test_cases
import numpy as np
from numpy.linalg import norm
from scipy.sparse import coo_matrix,triu,tril,eye

NUM_TRIALS = 25
EPS = 2.5 # %
TOL = 1e-4

class TestConstraints(unittest.TestCase):

    def setUp(self):

        # Network
        self.T = 2

        # Random
        np.random.seed(0)

    def test_constr_FIX(self):
        
        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # add vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = vargen.index*1.5
                vargen.Q = vargen.index*2.5
            self.assertGreater(net.num_var_generators,0)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
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
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             2*net.num_buses +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts() +
                             net.num_var_generators*2+
                             3*net.num_batteries+
                             net.num_loads)

            # Fixed
            net.set_flags('bus',
                          'fixed',
                          'slack',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('bus',
                          'fixed',
                          'regulated by generator',
                          'voltage magnitude')
            net.set_flags('generator',
                          'fixed',
                          'regulator',
                          'reactive power')
            net.set_flags('branch',
                          'fixed',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'fixed',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'fixed',
                          'switching - v',
                          'susceptance')
            net.set_flags('variable generator',
                          'fixed',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'fixed',
                          'any',
                          ['charging power','energy level'])
            net.set_flags('load',
                          'fixed',
                          'any',
                          'active power')
            self.assertGreater(net.num_fixed,0)
            self.assertEqual(net.num_fixed,
                             2*(net.get_num_slack_buses()) +
                             (net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts() +
                             net.num_var_generators*2+
                             3*net.num_batteries+
                             net.num_loads)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            constr = pf.Constraint('variable fixing',net)
            self.assertEqual(constr.name,'variable fixing')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            A_nnz = net.num_fixed+net.get_num_buses_reg_by_gen()
            constr.analyze()
            self.assertEqual(A_nnz,constr.A_nnz)
            constr.eval(x0)
            self.assertEqual(0,constr.A_nnz)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Vargen
            for vargen in net.var_generators:
                ar = np.where(A.col == vargen.index_P)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],vargen.index_P)
                self.assertEqual(b[A.row[ar[0]]],vargen.P)
                self.assertEqual(b[A.row[ar[0]]],vargen.index*1.5)
            for vargen in net.var_generators:
                ar = np.where(A.col == vargen.index_Q)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],vargen.index_Q)
                self.assertEqual(b[A.row[ar[0]]],vargen.Q)
                self.assertEqual(b[A.row[ar[0]]],vargen.index*2.5)

            # Batteries
            for bat in net.batteries:
                ar = np.where(A.col == bat.index_Pc)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],bat.index_Pc)
                self.assertEqual(b[A.row[ar[0]]],max([bat.P,0]))
            for bat in net.batteries:
                ar = np.where(A.col == bat.index_Pd)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],bat.index_Pd)
                self.assertEqual(b[A.row[ar[0]]],max([-bat.P,0]))
            for bat in net.batteries:
                ar = np.where(A.col == bat.index_E)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],bat.index_E)
                self.assertEqual(b[A.row[ar[0]]],bat.E)

            # Load
            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertTrue(load.has_flags('fixed','active power'))
                ar = np.where(A.col == load.index_P)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],load.index_P)
                self.assertEqual(b[A.row[ar[0]]],load.P)

        
        # Multiperiods
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # add vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = np.random.rand(self.T)*10
                vargen.Q = np.random.rand(self.T)*10
                self.assertEqual(vargen.num_periods,self.T)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # Vars
            net.set_flags('bus',
                          ['variable','fixed'],
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          ['variable','fixed'],
                          'slack',
                          'active power')
            net.set_flags('generator',
                          ['variable','fixed'],
                          'regulator',
                          'reactive power')
            net.set_flags('branch',
                          ['variable','fixed'],
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          ['variable','fixed'],
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          ['variable','fixed'],
                          'switching - v',
                          'susceptance')
            net.set_flags('variable generator',
                          ['variable','fixed'],
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          ['variable','fixed'],
                          'any',
                          ['charging power','energy level'])
            net.set_flags('load',
                          ['variable','fixed'],
                          'any',
                          'active power')
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,net.num_vars)
            self.assertEqual(net.num_vars,
                             (2*net.num_buses +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              net.num_var_generators*2+
                              3*net.num_batteries+
                              net.num_loads)*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            constr = pf.Constraint('variable fixing',net)
            self.assertEqual(constr.name,'variable fixing')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            A_nnz = net.num_fixed+net.get_num_buses_reg_by_gen()*self.T
            constr.analyze()
            self.assertEqual(A_nnz,constr.A_nnz)
            constr.eval(x0)
            self.assertEqual(0,constr.A_nnz)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(net.num_fixed,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(net.num_fixed,net.num_vars))
            self.assertEqual(A.nnz,net.num_fixed+net.get_num_buses_reg_by_gen()*self.T)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Time loop
            for t in range(self.T):

                # bus
                for bus in net.buses:
                    ar = np.where(A.col == bus.index_v_mag[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],bus.index_v_mag[t])
                    if bus.is_regulated_by_gen():
                        self.assertEqual(b[A.row[ar[0]]],bus.v_set[t])
                    else:
                        self.assertEqual(b[A.row[ar[0]]],bus.v_mag[t])
                    ar = np.where(A.col == bus.index_v_ang[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],bus.index_v_ang[t])
                    self.assertEqual(b[A.row[ar[0]]],bus.v_ang[t])

                # Gens
                for gen in net.generators:
                    if gen.is_slack():
                        ar = np.where(A.col == gen.index_P[t])[0]
                        self.assertEqual(ar.size,1)
                        self.assertEqual(A.col[ar[0]],gen.index_P[t])
                        self.assertEqual(b[A.row[ar[0]]],gen.P[t])
                    if gen.is_regulator():
                        ar = np.where(A.col == gen.index_Q[t])[0]
                        if gen.index == gen.reg_bus.reg_generators[0].index:
                            self.assertEqual(ar.size,2)
                        else:
                            self.assertEqual(ar.size,1)
                        self.assertEqual(A.col[ar[0]],gen.index_Q[t])
                        for i in range(ar.size):
                            if A.data[ar[i]] == 1.:
                                self.assertEqual(b[A.row[ar[i]]],gen.Q[t])
                            else:
                                self.assertEqual(A.data[ar[i]],0)

                # Shunts
                for shunt in net.shunts:
                    if shunt.is_switched_v():
                        ar = np.where(A.col == shunt.index_b[t])[0]
                        self.assertEqual(ar.size,1)
                        self.assertEqual(A.col[ar[0]],shunt.index_b[t])
                        self.assertEqual(b[A.row[ar[0]]],shunt.b[t])

                # Branch
                for branch in net.branches:
                    if branch.is_tap_changer():
                        ar = np.where(A.col == branch.index_ratio[t])[0]
                        self.assertEqual(ar.size,1)
                        self.assertEqual(A.col[ar[0]],branch.index_ratio[t])
                        self.assertEqual(b[A.row[ar[0]]],branch.ratio[t])
                    if branch.is_phase_shifter():
                        ar = np.where(A.col == branch.index_phase[t])[0]
                        self.assertEqual(ar.size,1)
                        self.assertEqual(A.col[ar[0]],branch.index_phase[t])
                        self.assertEqual(b[A.row[ar[0]]],branch.phase[t])

                # Vargen
                for vargen in net.var_generators:
                    ar = np.where(A.col == vargen.index_P[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],vargen.index_P[t])
                    self.assertEqual(b[A.row[ar[0]]],vargen.P[t])
                    ar = np.where(A.col == vargen.index_Q[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],vargen.index_Q[t])
                    self.assertEqual(b[A.row[ar[0]]],vargen.Q[t])

                # Batteries
                for bat in net.batteries:
                    ar = np.where(A.col == bat.index_Pc[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],bat.index_Pc[t])
                    self.assertEqual(b[A.row[ar[0]]],max([bat.P[t],0]))
                    ar = np.where(A.col == bat.index_Pd[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],bat.index_Pd[t])
                    self.assertEqual(b[A.row[ar[0]]],max([-bat.P[t],0]))
                    ar = np.where(A.col == bat.index_E[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],bat.index_E[t])
                    self.assertEqual(b[A.row[ar[0]]],bat.E[t])

                # Load
                for load in net.loads:
                    self.assertTrue(load.has_flags('variable','active power'))
                    self.assertTrue(load.has_flags('fixed','active power'))
                    ar = np.where(A.col == load.index_P[t])[0]
                    self.assertEqual(ar.size,1)
                    self.assertEqual(A.col[ar[0]],load.index_P[t])
                    self.assertEqual(b[A.row[ar[0]]],load.P[t])

    def test_constr_LBOUND(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # add vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = vargen.index*1.5
                vargen.Q = vargen.index*2.5
                vargen.P_ava = vargen.index*3.
                vargen.P_max = 100.
                vargen.P_min = 0.
                vargen.Q_max = 50.
                vargen.Q_min = -50.
            self.assertGreater(net.num_var_generators,0)

            self.assertEqual(net.num_bounded,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)
                load.P_max = 3.3*(load.index+1)
                load.Q = 3.5*load.index
                
            # Vars
            net.set_flags('bus',
                          'variable',
                          'regulated by generator',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          ['active power','reactive power'])
            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          ['active power','reactive power'])
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
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            num_vars_saved = net.num_vars
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             net.get_num_buses_reg_by_gen()*2 +
                             net.get_num_reg_gens()*2 +
                             2*net.get_num_P_adjust_loads() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters()*1 +
                             net.get_num_switched_shunts() +
                             net.num_var_generators*2+
                             3*net.num_batteries)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            constr = pf.Constraint('variable bounds',net)
            self.assertEqual(constr.name,'variable bounds')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            constr.analyze()
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            constr.eval(x0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(net.num_vars,net.num_vars))
            self.assertEqual(G.nnz,net.num_vars)
            self.assertTrue(np.all(G.row == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.col == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.data == np.ones(net.num_vars)))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(net.num_vars,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(net.num_vars,))

            E = G-eye(net.num_vars)
            self.assertGreater(G.nnz,0)
            self.assertGreater(norm(G.data,np.inf),0.5)
            self.assertEqual(E.nnz,0)

            self.assertTrue(not np.any(np.isinf(l)))
            self.assertTrue(not np.any(np.isnan(l)))
            self.assertTrue(not np.any(np.isinf(u)))
            self.assertTrue(not np.any(np.isnan(u)))
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Bounds
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags('variable',
                                                  ['voltage magnitude',
                                                   'voltage angle']))
                    self.assertEqual(u[bus.index_v_mag],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_v_ang],pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_v_mag],0.)
                    self.assertEqual(l[bus.index_v_ang],-pf.BUS_INF_V_ANG)
                else:
                    self.assertFalse(bus.has_flags('variable',
                                                   ['voltage magnitude',
                                                    'voltage angle']))

            for branch in net.branches:
                if branch.is_tap_changer():
                    self.assertTrue(branch.has_flags('variable','tap ratio'))
                    self.assertEqual(u[branch.index_ratio],pf.BRANCH_INF_RATIO)
                    self.assertEqual(l[branch.index_ratio],0.)
                else:
                    self.assertFalse(branch.has_flags('variable','tap ratio'))
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags('variable','phase shift'))
                    self.assertLess(np.abs(u[branch.index_phase]-np.pi*2.),1e-10)
                    self.assertLess(np.abs(l[branch.index_phase]+np.pi*2.),1e-10)
                else:
                    self.assertFalse(branch.has_flags('variable','phase shift'))

            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags('variable',['active power','reactive power']))
                    self.assertEqual(u[gen.index_P],pf.GEN_INF_P)
                    self.assertEqual(u[gen.index_Q],pf.GEN_INF_Q)
                    self.assertEqual(l[gen.index_P],-pf.GEN_INF_P)
                    self.assertEqual(l[gen.index_Q],-pf.GEN_INF_Q)
                else:
                    self.assertFalse(gen.has_flags('variable',
                                                   ['active power','reactive power']))

            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertTrue(load.has_flags('variable','reactive power'))
                self.assertTrue(load.has_flags('variable',['active power','reactive power']))
                self.assertEqual(u[load.index_P],pf.LOAD_INF_P)
                self.assertEqual(l[load.index_P],-pf.LOAD_INF_P)
                self.assertEqual(u[load.index_Q],pf.LOAD_INF_Q)
                self.assertEqual(l[load.index_Q],-pf.LOAD_INF_Q)

            for vargen in net.var_generators:
                self.assertTrue(vargen.has_flags('variable',
                                                 ['active power','reactive power']))
                self.assertEqual(u[vargen.index_P],pf.VARGEN_INF_P)
                self.assertEqual(u[vargen.index_Q],pf.VARGEN_INF_Q)
                self.assertEqual(l[vargen.index_P],-pf.VARGEN_INF_P)
                self.assertEqual(l[vargen.index_Q],-pf.VARGEN_INF_Q)

            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags('variable','susceptance'))
                    self.assertEqual(u[shunt.index_b],pf.SHUNT_INF_SUSC)
                    self.assertEqual(l[shunt.index_b],-pf.SHUNT_INF_SUSC)
                else:
                    self.assertFalse(shunt.has_flags('variable','susceptance'))

            for bat in net.batteries:
                self.assertTrue(bat.has_flags('variable','charging power'))
                self.assertTrue(bat.has_flags('variable','energy level'))
                self.assertEqual(u[bat.index_Pc],pf.BAT_INF_P)
                self.assertEqual(l[bat.index_Pc],0.)
                self.assertEqual(u[bat.index_Pd],pf.BAT_INF_P)
                self.assertEqual(l[bat.index_Pd],0.)
                self.assertEqual(u[bat.index_E],pf.BAT_INF_E)
                self.assertEqual(l[bat.index_E],0.)

            # Add bounded flags
            net.set_flags('bus',
                          'bounded',
                          'regulated by generator',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'bounded',
                          'regulator',
                          ['active power','reactive power'])
            net.set_flags('load',
                          'bounded',
                          'adjustable active power',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'bounded',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'bounded',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'bounded',
                          'switching - v',
                          'susceptance')
            net.set_flags('variable generator',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'bounded',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,num_vars_saved)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,net.num_vars)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            constr = pf.Constraint('variable bounds',net)
            self.assertEqual(constr.name,'variable bounds')

            constr.analyze()

            G = constr.G
            l = constr.l
            u = constr.u

            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(net.num_vars,net.num_vars))
            self.assertEqual(G.nnz,net.num_vars)
            self.assertTrue(np.all(G.row == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.col == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.data == np.ones(net.num_vars)))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(net.num_vars,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(net.num_vars,))

            E = G-eye(net.num_vars)
            self.assertGreater(G.nnz,0)
            self.assertGreater(norm(G.data,np.inf),0.5)
            self.assertEqual(E.nnz,0)

            # Bounds
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags('bounded',
                                                  ['voltage magnitude',
                                                   'voltage angle']))
                    self.assertTrue(bus.has_flags('variable',
                                                  ['voltage magnitude',
                                                   'voltage angle']))
                    self.assertEqual(u[bus.index_v_mag],bus.v_max)
                    self.assertEqual(u[bus.index_v_ang],pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_v_mag],bus.v_min)
                    self.assertEqual(l[bus.index_v_ang],-pf.BUS_INF_V_ANG)
                else:
                    self.assertFalse(bus.has_flags('bounded',
                                                   ['voltage magnitude',
                                                    'voltage angle']))

            for branch in net.branches:
                if branch.is_tap_changer():
                    self.assertTrue(branch.has_flags('bounded','tap ratio'))
                    self.assertEqual(u[branch.index_ratio],branch.ratio_max)
                    self.assertEqual(l[branch.index_ratio],branch.ratio_min)
                else:
                    self.assertFalse(branch.has_flags('bounded','tap ratio'))
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags('bounded','phase shift'))
                    self.assertEqual(u[branch.index_phase],branch.phase_max)
                    self.assertEqual(l[branch.index_phase],branch.phase_min)
                else:
                    self.assertFalse(branch.has_flags('bounded','phase shift'))

            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags('bounded',['active power','reactive power']))
                    self.assertEqual(u[gen.index_P],gen.P_max)
                    self.assertEqual(u[gen.index_Q],gen.Q_max)
                    self.assertEqual(l[gen.index_P],gen.P_min)
                    self.assertEqual(l[gen.index_Q],gen.Q_min)
                else:
                    self.assertFalse(gen.has_flags('bounded',['active power','reactive power']))

            for load in net.loads:
                self.assertTrue(load.has_flags('bounded','active power'))
                self.assertTrue(load.has_flags('bounded','reactive power'))
                self.assertTrue(load.has_flags('bounded',['active power','reactive power']))
                self.assertEqual(u[load.index_P],load.P_max)
                self.assertEqual(l[load.index_P],load.P_min)
                self.assertEqual(u[load.index_Q],pf.LOAD_INF_Q)
                self.assertEqual(l[load.index_Q],-pf.LOAD_INF_Q)

            for vargen in net.var_generators:
                self.assertTrue(vargen.has_flags('bounded',['active power','reactive power']))
                self.assertEqual(u[vargen.index_P],vargen.P_ava)
                self.assertEqual(u[vargen.index_Q],vargen.Q_max)
                self.assertEqual(l[vargen.index_P],vargen.P_min)
                self.assertEqual(l[vargen.index_Q],vargen.Q_min)

            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags('bounded','susceptance'))
                    self.assertEqual(u[shunt.index_b],shunt.b_max)
                    self.assertEqual(l[shunt.index_b],shunt.b_min)
                else:
                    self.assertFalse(shunt.has_flags('bounded','susceptance'))

            for bat in net.batteries:
                self.assertTrue(bat.has_flags('bounded','charging power'))
                self.assertTrue(bat.has_flags('bounded','energy level'))
                self.assertEqual(u[bat.index_Pc],bat.P_max)
                self.assertEqual(l[bat.index_Pc],0.)
                self.assertEqual(u[bat.index_Pd],-bat.P_min)
                self.assertEqual(l[bat.index_Pd],0.)
                self.assertEqual(u[bat.index_E],bat.E_max)
                self.assertEqual(l[bat.index_E],0.)

            # Sensitivities
            net.clear_sensitivities()
            for bus in net.buses:
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
                self.assertEqual(bus.sens_v_mag_u_bound,0.)
                self.assertEqual(bus.sens_v_mag_l_bound,0.)
                self.assertEqual(bus.sens_v_ang_u_bound,0.)
                self.assertEqual(bus.sens_v_ang_l_bound,0.)
            for gen in net.generators:
                self.assertEqual(gen.sens_P_u_bound,0.)
                self.assertEqual(gen.sens_P_l_bound,0.)
            for load in net.loads:
                self.assertEqual(load.sens_P_u_bound,0.)
                self.assertEqual(load.sens_P_l_bound,0.)

            mu = np.random.randn(net.num_vars)
            pi = np.random.randn(net.num_vars)

            constr.store_sensitivities(None,None,mu,pi)

            for bus in net.buses:
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
                self.assertEqual(bus.sens_v_mag_u_bound,0.)
                self.assertEqual(bus.sens_v_mag_l_bound,0.)
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags('variable','voltage angle'))
                    self.assertNotEqual(bus.sens_v_ang_u_bound,0.)
                    self.assertNotEqual(bus.sens_v_ang_l_bound,0.)
                    self.assertEqual(bus.sens_v_ang_u_bound,mu[bus.index_v_ang])
                    self.assertEqual(bus.sens_v_ang_l_bound,pi[bus.index_v_ang])
                else:
                    self.assertEqual(bus.sens_v_ang_u_bound,0.)
                    self.assertEqual(bus.sens_v_ang_l_bound,0.)
            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags('variable','active power'))
                    self.assertNotEqual(gen.sens_P_u_bound,0.)
                    self.assertNotEqual(gen.sens_P_l_bound,0.)
                    self.assertEqual(gen.sens_P_u_bound,mu[gen.index_P])
                    self.assertEqual(gen.sens_P_l_bound,pi[gen.index_P])
                else:
                    self.assertEqual(gen.sens_P_u_bound,0.)
                    self.assertEqual(gen.sens_P_l_bound,0.)
            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertNotEqual(load.sens_P_u_bound,0.)
                self.assertNotEqual(load.sens_P_l_bound,0.)
                self.assertEqual(load.sens_P_u_bound,mu[load.index_P])
                self.assertEqual(load.sens_P_l_bound,pi[load.index_P])

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # add vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = np.random.rand(self.T)
                vargen.Q = np.random.rand(self.T)
                vargen.P_ava = vargen.P*3.4
                vargen.P_max = 100.
                vargen.P_min = 0.
                vargen.Q_max = 50.
                vargen.Q_min = -50.
                self.assertEqual(vargen.num_periods,self.T)
                for t in range(self.T):
                    self.assertEqual(vargen.P_ava[t],vargen.P[t]*3.4)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_bounded,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)*np.array(range(net.num_periods))
                load.P_max = 3.3*(load.index+1)*np.array(range(net.num_periods))
                load.Q = 3.5*load.index*np.array(range(net.num_periods))

            # Vars
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
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          ['tap ratio'])
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          ['susceptance'])
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             (net.num_buses*2 +
                              net.num_generators*2 +
                              2*net.num_loads +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              net.num_var_generators*2+
                              3*net.num_batteries)*self.T)

            x0 = net.get_var_values()
            constr = pf.Constraint('variable bounds',net)
            self.assertEqual(constr.name,'variable bounds')
            constr.analyze()
            constr.eval(x0)

            G = constr.G
            l = constr.l
            u = constr.u

            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(net.num_vars,net.num_vars))
            self.assertEqual(G.nnz,net.num_vars)
            self.assertTrue(np.all(G.row == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.col == np.array(range(net.num_vars))))
            self.assertTrue(np.all(G.data == np.ones(net.num_vars)))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(net.num_vars,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(net.num_vars,))

            for t in range(self.T):
                for bus in net.buses:
                    self.assertEqual(u[bus.index_v_mag[t]],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_v_ang[t]],pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_v_mag[t]],0)
                    self.assertEqual(l[bus.index_v_ang[t]],-pf.BUS_INF_V_ANG)
                for gen in net.generators:
                    self.assertEqual(u[gen.index_P[t]],pf.GEN_INF_P)
                    self.assertEqual(u[gen.index_Q[t]],pf.GEN_INF_Q)
                    self.assertEqual(l[gen.index_P[t]],-pf.GEN_INF_P)
                    self.assertEqual(l[gen.index_Q[t]],-pf.GEN_INF_Q)
                for branch in net.branches:
                    if branch.is_tap_changer():
                        self.assertEqual(u[branch.index_ratio[t]],pf.BRANCH_INF_RATIO)
                        self.assertEqual(l[branch.index_ratio[t]],0.)
                    if branch.is_phase_shifter():
                        self.assertLess(np.abs(u[branch.index_phase[t]]-np.pi*2.),1e-10)
                        self.assertLess(np.abs(l[branch.index_phase[t]]+np.pi*2.),1e-10)
                for vargen in net.var_generators:
                    self.assertEqual(u[vargen.index_P[t]],pf.VARGEN_INF_P)
                    self.assertEqual(u[vargen.index_Q[t]],pf.VARGEN_INF_Q)
                    self.assertEqual(l[vargen.index_P[t]],-pf.VARGEN_INF_P)
                    self.assertEqual(l[vargen.index_Q[t]],-pf.VARGEN_INF_Q)
                for load in net.loads:
                    self.assertEqual(u[load.index_P[t]],pf.LOAD_INF_P)
                    self.assertEqual(l[load.index_P[t]],-pf.LOAD_INF_P)
                    self.assertEqual(u[load.index_Q[t]],pf.LOAD_INF_Q)
                    self.assertEqual(l[load.index_Q[t]],-pf.LOAD_INF_Q)
                for shunt in net.shunts:
                    if shunt.is_switched_v():
                        self.assertEqual(u[shunt.index_b[t]],pf.SHUNT_INF_SUSC)
                        self.assertEqual(l[shunt.index_b[t]],-pf.SHUNT_INF_SUSC)

            # Bounded
            net.set_flags('bus',
                          'bounded',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('load',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'bounded',
                          'tap changer',
                          ['tap ratio'])
            net.set_flags('branch',
                          'bounded',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'bounded',
                          'switching - v',
                          ['susceptance'])
            net.set_flags('variable generator',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'bounded',
                          'any',
                          ['charging power','energy level'])
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_bounded,net.num_vars)

            x0 = net.get_var_values()
            constr = pf.Constraint('variable bounds',net)
            self.assertEqual(constr.name,'variable bounds')
            constr.analyze()
            constr.eval(x0)

            G = constr.G
            l = constr.l
            u = constr.u

            for t in range(self.T):
                for bus in net.buses:
                    self.assertEqual(u[bus.index_v_mag[t]],bus.v_max)
                    self.assertEqual(u[bus.index_v_ang[t]],pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_v_mag[t]],bus.v_min)
                    self.assertEqual(l[bus.index_v_ang[t]],-pf.BUS_INF_V_ANG)
                for gen in net.generators:
                    self.assertEqual(u[gen.index_P[t]],gen.P_max)
                    self.assertEqual(u[gen.index_Q[t]],gen.Q_max)
                    self.assertEqual(l[gen.index_P[t]],gen.P_min)
                    self.assertEqual(l[gen.index_Q[t]],gen.Q_min)
                for branch in net.branches:
                    if branch.is_tap_changer():
                        self.assertEqual(u[branch.index_ratio[t]],branch.ratio_max)
                        self.assertEqual(l[branch.index_ratio[t]],branch.ratio_min)
                    if branch.is_phase_shifter():
                        self.assertEqual(u[branch.index_phase[t]],branch.phase_max)
                        self.assertEqual(l[branch.index_phase[t]],branch.phase_min)
                for vargen in net.var_generators:
                    self.assertEqual(u[vargen.index_P[t]],vargen.P_ava[t])
                    self.assertEqual(u[vargen.index_Q[t]],vargen.Q_max)
                    self.assertEqual(l[vargen.index_P[t]],vargen.P_min)
                    self.assertEqual(l[vargen.index_Q[t]],vargen.Q_min)
                for load in net.loads:
                    self.assertEqual(u[load.index_P[t]],load.P_max[t])
                    self.assertEqual(l[load.index_P[t]],load.P_min[t])
                    self.assertEqual(u[load.index_P[t]],3.3*(load.index+1)*t)
                    self.assertEqual(l[load.index_P[t]],-2.4*(load.index+1)*t)
                    self.assertEqual(u[load.index_Q[t]],pf.LOAD_INF_Q)
                    self.assertEqual(l[load.index_Q[t]],-pf.LOAD_INF_Q)
                for shunt in net.shunts:
                    if shunt.is_switched_v():
                        self.assertEqual(u[shunt.index_b[t]],shunt.b_max)
                        self.assertEqual(l[shunt.index_b[t]],shunt.b_min)

    def test_constr_PAR_GEN_P(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            self.assertEqual(net.num_vars,0)

            # Vars
            net.set_flags('generator',
                          'variable',
                          'slack',
                          ['active power','reactive power'])
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_vars,(net.get_num_slack_gens()+net.get_num_reg_gens())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('generator active power participation',net)
            self.assertEqual(constr.name,'generator active power participation')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)

            # Manual count
            nnz = 0
            num_constr = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_slack():
                    num_constr += len(bus.generators)-1 # P participation
                    nnz += 2*(len(bus.generators)-1)

            constr.analyze()
            self.assertEqual(nnz*self.T,constr.A_nnz)
            constr.eval(x0)
            self.assertEqual(0,constr.A_nnz)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b

            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(num_constr*self.T,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(num_constr*self.T,net.num_vars))
            self.assertEqual(A.nnz,nnz*self.T)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Detailed check
            Ai = A.row
            Aj = A.col
            Ad = A.data
            self.assertEqual(Ai.size,nnz*self.T)
            self.assertEqual(Aj.size,nnz*self.T)
            self.assertEqual(Ad.size,nnz*self.T)
            i = 0
            row = 0
            counted = {}
            for t in range(self.T):
                for k in range(net.num_branches):
                    br = net.get_branch(k)
                    for bus in [br.bus_k,br.bus_m]:
                        if (bus.number,t) in counted:
                            continue
                        counted[(bus.number,t)] = True
                        if bus.is_slack():
                            gens = bus.generators
                            self.assertGreater(len(gens),0)
                            g1 = gens[0]
                            for g2 in gens[1:]:
                                self.assertEqual(b[row],0.)
                                self.assertEqual(Ai[i],row)
                                self.assertEqual(Aj[i],g1.index_P[t])
                                self.assertEqual(Ad[i],1.)
                                i += 1
                                self.assertEqual(Ai[i],row)
                                self.assertEqual(Aj[i],g2.index_P[t])
                                self.assertEqual(Ad[i],-1.)
                                i += 1
                                row += 1
            self.assertEqual(i,nnz*self.T)

            # Last check
            x = np.zeros(net.num_vars)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    if bus.is_slack():
                        self.assertGreater(len(bus.generators),0)
                        for g in bus.generators:
                            self.assertTrue(g.has_flags('variable','active power'))
                            x[g.index_P[t]] = 10.
            self.assertGreater(norm(x),0)
            self.assertTrue(norm(A*x-b) < 1e-10)

    def test_constr_PAR_GEN_Q(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            self.assertEqual(net.num_vars,0)

            # Vars
            net.set_flags('generator',
                          'variable',
                          'slack',
                          ['active power','reactive power'])
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_vars,(net.get_num_slack_gens()+net.get_num_reg_gens())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('generator reactive power participation',net)
            self.assertEqual(constr.name,'generator reactive power participation')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)

            # Manual count
            nnz = 0
            num_constr = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen():
                    num_constr += len(bus.reg_generators)-1 # Q participation
                    nnz += 2*(len(bus.reg_generators)-1)

            constr.analyze()
            self.assertEqual(nnz*self.T,constr.A_nnz)
            constr.eval(x0)
            self.assertEqual(0,constr.A_nnz)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b

            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(num_constr*self.T,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(num_constr*self.T,net.num_vars))
            self.assertEqual(A.nnz,nnz*self.T)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Detailed check
            Ai = A.row
            Aj = A.col
            Ad = A.data
            self.assertEqual(Ai.size,nnz*self.T)
            self.assertEqual(Aj.size,nnz*self.T)
            self.assertEqual(Ad.size,nnz*self.T)
            i = 0
            row = 0
            counted = {}
            for t in range(self.T):
                for k in range(net.num_branches):
                    br = net.get_branch(k)
                    for bus in [br.bus_k,br.bus_m]:
                        if (bus.number,t) in counted:
                            continue
                        counted[(bus.number,t)] = True
                        if bus.is_regulated_by_gen():
                            reg_gens = bus.reg_generators
                            self.assertGreater(len(reg_gens),0)
                            g1 = reg_gens[0]
                            self.assertGreater(g1.Q_max,g1.Q_min)
                            for g2 in reg_gens[1:]:
                                self.assertTrue(np.abs(b[row]-(g1.Q_min/(g1.Q_max-g1.Q_min)-g2.Q_min/(g2.Q_max-g2.Q_min))) < 1e-10)
                                self.assertGreater(g2.Q_max,g2.Q_min)
                                self.assertEqual(Ai[i],row)
                                self.assertEqual(Aj[i],g1.index_Q[t])
                                self.assertTrue(np.abs(Ad[i]-1./(g1.Q_max-g1.Q_min)) < 1e-10)
                                i += 1
                                self.assertEqual(Ai[i],row)
                                self.assertEqual(Aj[i],g2.index_Q[t])
                                self.assertTrue(np.abs(Ad[i]+1./(g2.Q_max-g2.Q_min)) < 1e-10)
                                i += 1
                                row += 1
            self.assertEqual(i,nnz*self.T)

            # Last check
            x = np.zeros(net.num_vars)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    if bus.is_regulated_by_gen():
                        self.assertGreater(len(bus.reg_generators),0)
                        for g in bus.reg_generators:
                            self.assertTrue(g.has_flags('variable','reactive power'))
                            x[g.index_Q[t]] = (g.Q_max+g.Q_min)/2.
            self.assertTrue(norm(A*x-b) < 1e-10)

    def test_constr_ACPF(self):

        # Constants
        h = 1e-10

        # Multiperiods
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            
            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))
            for vargen in net.var_generators:
                vargen.P = np.random.rand(net.num_periods)
                vargen.Q = np.random.randn(net.num_periods)

            # Add batteries
            gen_buses = net.get_generator_buses()
            net.add_batteries(gen_buses,20.,40.,0.8,0.9)
            self.assertGreater(net.num_batteries,0)
            self.assertEqual(net.num_batteries,len(gen_buses))
            for bat in net.batteries:
                bat.P = np.random.randn(net.num_periods)

            # No vars
            self.assertEqual(net.num_vars,0)
            
            # Constraint
            constr = pf.Constraint('AC power balance',net)
            self.assertEqual(constr.name,'AC power balance')

            x0 = net.get_var_values()
            self.assertEqual(x0.size,0)
            constr.analyze()
            constr.eval(x0)
            
            f = constr.f
            self.assertEqual(f.size,2*net.num_buses*net.num_periods)

            # Check mismatches (no vars)
            for t in range(net.num_periods):
                for bus in net.buses:
                    P_mis = 0.
                    Q_mis = 0.
                    for branch in bus.branches_k:
                        P_mis -= branch.get_P_km()[t]
                        Q_mis -= branch.get_Q_km()[t]
                    for branch in bus.branches_m:
                        P_mis -= branch.get_P_mk()[t]
                        Q_mis -= branch.get_Q_mk()[t]
                    for gen in bus.generators:
                        P_mis += gen.P[t]
                        Q_mis += gen.Q[t]
                    for vargen in bus.var_generators:
                        P_mis += vargen.P[t]
                        Q_mis += vargen.Q[t]
                    for load in bus.loads:
                        P_mis -= load.P[t]
                        Q_mis -= load.Q[t]
                    for bat in bus.batteries:
                        P_mis -= bat.P[t]
                    for shunt in bus.shunts:
                        P_mis -= shunt.g*(bus.v_mag[t]**2.)
                        Q_mis -= -shunt.b[t]*(bus.v_mag[t]**2.)
                    self.assertAlmostEqual(P_mis,f[bus.index_P+t*2*net.num_buses])
                    self.assertAlmostEqual(Q_mis,f[bus.index_Q+t*2*net.num_buses])

            # Cross check mismatches with net properties (no vars)
            net.update_properties()
            dP_list = dict([(t,list()) for t in range(self.T)])
            dQ_list = dict([(t,list()) for t in range(self.T)])
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    dP = f[bus.index_P+t*2*net.num_buses]
                    dQ = f[bus.index_Q+t*2*net.num_buses]
                    dP_list[t].append(dP)
                    dQ_list[t].append(dQ)
                    self.assertAlmostEqual(dP,bus.P_mis[t])
                    self.assertAlmostEqual(dQ,bus.Q_mis[t])
            self.assertAlmostEqual(net.bus_P_mis[t],np.max(np.abs(dP_list[t]))*net.base_power)
            self.assertAlmostEqual(net.bus_Q_mis[t],np.max(np.abs(dQ_list[t]))*net.base_power)

            # Vars
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
                          ['active power','reactive power'])
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
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,
                             (2*net.num_buses +
                              2*net.num_generators +
                              2*net.num_loads +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              3*net.num_batteries +
                              net.num_var_generators*2)*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint('AC power balance',net)
            self.assertEqual(constr.name,'AC power balance')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            num_Jnnz = (net.num_buses*4 +
                        net.num_branches*8 +
                        net.get_num_tap_changers()*4 +
                        net.get_num_phase_shifters()*4 +
                        net.get_num_switched_shunts() +
                        net.num_generators*2 +
                        net.num_loads*2 +
                        net.num_batteries*2 +
                        net.num_var_generators*2)*self.T

            constr.analyze()
            self.assertEqual(num_Jnnz,constr.J_nnz)
            constr.eval(x0)
            self.assertEqual(num_Jnnz,constr.J_nnz)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            constr.combine_H(np.ones(f.size),False)
            Hcomb = constr.H_combined

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(2*net.num_buses*self.T,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(2*net.num_buses*self.T,net.num_vars))
            self.assertEqual(J.nnz,num_Jnnz)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertTupleEqual(Hcomb.shape,(net.num_vars,net.num_vars))
            self.assertEqual(Hcomb.nnz,2*(net.get_num_buses()*3 +
                                          net.get_num_branches()*12 +
                                          net.get_num_tap_changers()*9 +
                                          net.get_num_phase_shifters()*10 +
                                          net.get_num_switched_shunts())*self.T)

            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))

            # Check mismatches
            x1 = x0+np.random.randn(net.num_vars)
            constr.eval(x1)
            for t in range(net.num_periods):
                for bus in net.buses:
                    P_mis = 0.
                    Q_mis = 0.
                    for branch in bus.branches_k:
                        P_mis -= branch.get_P_km(x1)[t]
                        Q_mis -= branch.get_Q_km(x1)[t]
                    for branch in bus.branches_m:
                        P_mis -= branch.get_P_mk(x1)[t]
                        Q_mis -= branch.get_Q_mk(x1)[t]
                    for gen in bus.generators:
                        P_mis += x1[gen.index_P[t]]
                        Q_mis += x1[gen.index_Q[t]]
                    for vargen in bus.var_generators:
                        P_mis += x1[vargen.index_P[t]]
                        Q_mis += x1[vargen.index_Q[t]]
                    for load in bus.loads:
                        P_mis -= x1[load.index_P[t]]
                        Q_mis -= x1[load.index_Q[t]]
                    for bat in bus.batteries:
                        P_mis -= x1[bat.index_Pc[t]]-x1[bat.index_Pd[t]]
                    for shunt in bus.shunts:
                        if shunt.has_flags('variable','susceptance'):
                            b = x1[shunt.index_b[t]]
                        else:
                            b = shunt.b[t]
                        if bus.has_flags('variable','voltage magnitude'):
                            v = x1[bus.index_v_mag[t]]
                        else:
                            v = bus.v_mag[t]
                        P_mis -= shunt.g*v*v
                        Q_mis -= -b*v*v
                    self.assertAlmostEqual(P_mis,f[bus.index_P+t*2*net.num_buses])
                    self.assertAlmostEqual(Q_mis,f[bus.index_Q+t*2*net.num_buses])

            # Cross check mismatches with net properties
            constr.eval(x1)
            net.update_properties(x1)
            dP_list = dict([(t,list()) for t in range(self.T)])
            dQ_list = dict([(t,list()) for t in range(self.T)])
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    dP = f[bus.index_P+t*2*net.num_buses]
                    dQ = f[bus.index_Q+t*2*net.num_buses]
                    dP_list[t].append(dP)
                    dQ_list[t].append(dQ)
                    self.assertAlmostEqual(dP,bus.P_mis[t])
                    self.assertAlmostEqual(dQ,bus.Q_mis[t])
            self.assertAlmostEqual(net.bus_P_mis[t],np.max(np.abs(dP_list[t]))*net.base_power)
            self.assertAlmostEqual(net.bus_Q_mis[t],np.max(np.abs(dQ_list[t]))*net.base_power)

            # Check mismatches across time
            for vargen in net.var_generators:
                vargen.P = np.ones(net.num_periods)*0.2 # static
                vargen.Q = np.ones(net.num_periods)*0.1 # static
            for bat in net.batteries:
                bat.P = np.ones(net.num_periods)*0.1    # static
            x0 = net.get_var_values()
            constr.eval(x0)
            f = constr.f
            J = constr.J
            P_list = []
            for t in range(self.T):
                P_list.append(net.get_var_projection('all','all',t_start=t,t_end=t))
            f_list = [f[t*2*net.num_buses:(t+1)*2*net.num_buses] for t in range(self.T)]
            for t in range(self.T-1):
                self.assertLess(norm(f_list[t]-f_list[t+1]),1e-12*norm(f_list[t]))
            Jx = J*x0
            Jx_list = [Jx[t*2*net.num_buses:(t+1)*2*net.num_buses] for t in range(self.T)]
            for t in range(self.T-1):
                self.assertLess(norm(Jx_list[t]-Jx_list[t+1]),1e-12*norm(Jx_list[t]))
            for i in range(10):
                H_list = []
                j = np.random.randint(0,2*net.num_buses)
                for t in range(self.T):
                    H_list.append(coo_matrix(P_list[t]*constr.get_H_single(t*2*net.num_buses+j)*P_list[t].T))
                for t in range(self.T-1):
                    self.assertTrue(np.all(H_list[t].row == H_list[t+1].row))
                    self.assertTrue(np.all(H_list[t].col == H_list[t+1].col))
                    self.assertLess(norm(H_list[t].data-H_list[t+1].data),1e-12*norm(H_list[t].data))

            # Jacobian check
            x0 = x0.copy()
            constr.eval(x0)
            f0 = constr.f.copy()
            J0 = constr.J.copy()
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars)

                x = x0 + h*d

                constr.eval(x)
                f1 = constr.f

                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_P_balance[t],0.)
                    self.assertEqual(bus.sens_Q_balance[t],0.)
            sens = np.zeros(2*net.num_buses*self.T)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.index_P,2*bus.index)
                    self.assertEqual(bus.index_Q,2*bus.index+1)
                    sens[bus.index_P+t*2*net.num_buses] = 3.5*bus.index_P+0.33+t*2*net.num_buses
                    sens[bus.index_Q+t*2*net.num_buses] = 3.4*bus.index_Q+0.32+t*2*net.num_buses
            constr.store_sensitivities(None,sens,None,None)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_P_balance[t],3.5*bus.index_P+0.33+t*2*net.num_buses)
                    self.assertEqual(bus.sens_Q_balance[t],3.4*bus.index_Q+0.32+t*2*net.num_buses)

    def test_constr_REG_GEN(self):

        # Constants
        h = 1e-8

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
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
            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-net.get_num_slack_buses()) +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('voltage regulation by generators',net)
            self.assertEqual(constr.name,'voltage regulation by generators')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.num_extra_vars,0)

            Jnnz = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    Jnnz += 2 + 2*len(bus.reg_generators)

            Annz = 3*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())

            rowsJ = 2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            rowsA = net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()

            constr.analyze()
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,Annz*self.T)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,rowsA*self.T)
            self.assertEqual(constr.num_extra_vars,rowsJ*self.T)
            
            y0 = np.random.rand(constr.num_extra_vars)
            constr.eval(x0,y0)
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ*self.T,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA*self.T,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(J.nnz,Jnnz*self.T)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(A.nnz,Annz*self.T)

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            # Ax=b check
            self.assertEqual(norm(A.data,1),rowsA*3*self.T)
            self.assertEqual(np.sum(A.data),(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())*self.T)
            for k in range(J.shape[0]//2):
                index1 = np.where(A.col == net.num_vars+2*k)[0]
                index2 = np.where(A.col == net.num_vars+2*k+1)[0]
                self.assertEqual(index1.size,1)
                self.assertEqual(index2.size,1)
                self.assertEqual(A.row[index1[0]],A.row[index2[0]])
                index3 = np.where(A.row == A.row[index1[0]])[0]
                self.assertEqual(index3.size,3)
                for i in index3:
                    if A.col[i] == net.num_vars+2*k: # y
                        self.assertEqual(A.data[i],-1.)
                    elif A.col[i] == net.num_vars+2*k+1:
                        self.assertEqual(A.data[i],1.) # z
                    else:
                        self.assertEqual(A.data[i],1.) # v

            # f check
            flags = {}
            eps = 1e-8
            for t in range(self.T):
                for bus in net.buses:
                    flags[(t,bus.index)] = False
            J_row = 0
            for t in range(self.T):
                for branch in net.branches:
                    for bus in [branch.bus_k,branch.bus_m]:
                        if not flags[(t,bus.index)]:
                            if bus.is_regulated_by_gen() and not bus.is_slack():
                                y = y0[J_row]
                                z = y0[J_row+1]
                                Q = sum([g.Q[t] for g in bus.reg_generators])
                                Qmax = sum([g.Q_max for g in bus.reg_generators])
                                Qmin = sum([g.Q_min for g in bus.reg_generators])
                                CompY = (Q-Qmin)+y-np.sqrt((Q-Qmin)**2.+y**2.+2*eps)
                                CompZ = (Qmax-Q)+z-np.sqrt((Qmax-Q)**2.+z**2.+2*eps)
                                self.assertAlmostEqual(CompY,f[J_row])
                                self.assertAlmostEqual(CompZ,f[J_row+1])
                                J_row += 2
                        flags[(t,bus.index)] = True            

            # Jacobian check
            f0 = f.copy()
            J0 = J.copy()
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                f1 = constr.f

                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if f.shape[0] == 0:
                    break

                j = np.random.randint(0,f.shape[0])

                constr.eval(x0,y0)

                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular

                H0 = (H0 + H0.T) - triu(H0)

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0,y0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,
                                  (net.num_vars+constr.num_extra_vars,
                                   net.num_vars+constr.num_extra_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                
                g1 = constr.J.T*coeff

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_v_reg_by_gen[t],0.)
            sens = np.zeros(constr.f.size)
            self.assertEqual(sens.size,rowsJ*self.T)
            Ji = constr.J.row
            Jj = constr.J.col
            
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    if bus.is_regulated_by_gen() and not bus.is_slack():
                        indices = Ji[np.where(Jj == bus.reg_generators[0].index_Q[t])[0]]
                        self.assertEqual(indices.size,2)
                        sens[indices[0]] = -bus.index-10
                        sens[indices[1]] = bus.index+11*(bus.index % 2)
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    if bus.is_regulated_by_gen() and not bus.is_slack():
                        if bus.index % 2 == 1:
                            self.assertEqual(bus.sens_v_reg_by_gen[t],bus.index+11)
                        else:
                            self.assertEqual(bus.sens_v_reg_by_gen[t],-bus.index-10)

    def test_constr_NBOUND(self):

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
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
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
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts())*self.T)

            # Bound constraints
            net.set_flags('bus',
                          'bounded',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'bounded',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'bounded',
                          'phase shifter',
                          'phase shift')
            net.set_flags('shunt',
                          'bounded',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_bounded,
                             (2*net.num_buses +
                              2*net.num_generators +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('variable nonlinear bounds',net)
            self.assertEqual(constr.name,'variable nonlinear bounds')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)

            J_nnz = 2*net.num_vars
            constr.analyze()
            self.assertEqual(J_nnz,constr.J_nnz)
            constr.eval(x0)
            self.assertEqual(J_nnz,constr.J_nnz)
            constr.store_sensitivities(None,np.zeros(constr.J.shape[0]),None,None)
            self.assertEqual(J_nnz,constr.J_nnz)
            constr.eval(x0)
            self.assertEqual(J_nnz,constr.J_nnz)

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
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
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
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_v_mag_u_bound[t],0.)
                    self.assertEqual(bus.sens_v_mag_l_bound[t],0.)
            sens = np.zeros(constr.f.size)
            Ji = constr.J.row
            Jj = constr.J.col
            for t in range(self.T):
                for i in range(net.num_buses): # buses
                    bus = net.get_bus(i)
                    indices = Ji[np.where(Jj == bus.index_v_mag[t])[0]]
                    self.assertEqual(indices.size,2)
                    sens[indices[0]] = bus.index*10.
                    sens[indices[1]] = -bus.index*10.
            constr.store_sensitivities(None,sens,None,None)
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_v_mag_u_bound[t],bus.index*10.)
                    self.assertEqual(bus.sens_v_mag_l_bound[t],-bus.index*10.)

    def test_constr_REG_TRAN(self):

        # Constants
        h = 1e-8
        normal = 1e0
        eta = 1e-8

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'regulated by transformer',
                          'voltage magnitude')
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            self.assertEqual(net.num_vars,
                             (net.get_num_buses_reg_by_tran() +
                              net.get_num_tap_changers_v())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('voltage regulation by transformers',net)
            self.assertEqual(constr.name,'voltage regulation by transformers')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)

            Jnnz = 10*net.get_num_tap_changers_v()
            Annz = 3*net.get_num_tap_changers_v()
            self.assertGreaterEqual(Jnnz,0)
            self.assertGreaterEqual(Annz,0)

            rowsJ = 4*net.get_num_tap_changers_v()
            rowsA = net.get_num_tap_changers_v()
            self.assertGreaterEqual(rowsJ,0)
            self.assertGreaterEqual(rowsA,0)

            constr.analyze()
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,Annz*self.T)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,rowsA*self.T)
            constr.eval(x0)
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            
            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ*self.T,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA*self.T,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(J.nnz,Jnnz*self.T)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(A.nnz,Annz*self.T)
            self.assertEqual(constr.num_extra_vars,rowsJ*self.T)

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            y0 = np.random.rand(constr.num_extra_vars)

            # Ax=b check
            self.assertEqual(norm(A.data,1),rowsA*3*self.T)
            self.assertEqual(np.sum(A.data),net.get_num_tap_changers_v()*self.T)

            # f check
            index = 0
            for t in range(self.T):
                for i in range(net.num_branches):
                    br = net.get_branch(i)
                    if br.is_tap_changer_v():
                        self.assertTrue(br.has_flags('variable','tap ratio'))
                        bus = br.reg_bus
                        fvmin = ((bus.v_mag[t]-bus.v_min_reg) - np.sqrt((bus.v_mag[t]-bus.v_min_reg)**2. + 2*eta))*normal
                        fvmax = ((bus.v_max_reg-bus.v_mag[t]) - np.sqrt((bus.v_max_reg-bus.v_mag[t])**2. + 2*eta))*normal
                        ftmax = ((br.ratio_max-br.ratio[t]) - np.sqrt((br.ratio_max-br.ratio[t])**2. + 2*eta))*normal
                        ftmin = ((br.ratio[t]-br.ratio_min) - np.sqrt((br.ratio[t]-br.ratio_min)**2. + 2*eta))*normal
                        self.assertLess(np.abs(fvmin-f[index]),1e-10*(1+np.abs(fvmin)))
                        self.assertLess(np.abs(fvmax-f[index+1]),1e-10*(1+np.abs(fvmax)))
                        self.assertLess(np.abs(ftmax-f[index+2]),1e-10*(1+np.abs(ftmax)))
                        self.assertLess(np.abs(ftmin-f[index+3]),1e-10*(1+np.abs(ftmin)))
                        index += 4

            # Jacobian check
            constr.eval(x0,y0)
            f0 = constr.f.copy()
            J0 = constr.J.copy()            
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                f1 = constr.f

                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if not f.size:
                    continue

                j = np.random.randint(0,f.shape[0])

                constr.eval(x0,y0)

                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular

                H0 = (H0 + H0.T) - triu(H0)

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0,y0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,
                                  (net.num_vars+constr.num_extra_vars,
                                   net.num_vars+constr.num_extra_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.T*coeff

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_v_reg_by_tran[t],0.)
            sens = np.zeros(constr.f.size)
            counter = 0
            for t in range(self.T):
                for branch in net.branches:
                    if branch.is_tap_changer():
                        sens[counter:counter+4] = branch.reg_bus.index*t
                        counter += 4
            self.assertEqual(counter,constr.f.size)
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
            for t in range(self.T):
                for branch in net.branches:
                    if branch.is_tap_changer():
                        self.assertEqual(branch.reg_bus.sens_v_reg_by_tran[t],branch.reg_bus.index*t)

    def test_constr_REG_SHUNT(self):

        # Constants
        h = 1e-8
        normal = 1e0
        eta = 1e-8

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'regulated by shunt',
                          'voltage magnitude')
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_vars,
                             (net.get_num_buses_reg_by_shunt() +
                              net.get_num_switched_shunts())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('voltage regulation by shunts',net)
            self.assertEqual(constr.name,'voltage regulation by shunts')

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

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)

            Jnnz = 10*net.get_num_switched_shunts()
            Annz = 3*net.get_num_switched_shunts()
            self.assertGreaterEqual(Jnnz,0)
            self.assertGreaterEqual(Annz,0)

            rowsJ = 4*net.get_num_switched_shunts()
            rowsA = net.get_num_switched_shunts()
            self.assertGreaterEqual(rowsJ,0)
            self.assertGreaterEqual(rowsA,0)

            constr.analyze()
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,Annz*self.T)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,rowsA*self.T)
            constr.eval(x0)
            self.assertEqual(constr.J_nnz,Jnnz*self.T)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,rowsJ*self.T)
            self.assertEqual(constr.A_row,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(rowsJ*self.T,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(rowsA*self.T,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(rowsJ*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(J.nnz,Jnnz*self.T)
            self.assertTrue(np.all(J.row <= rowsJ*self.T-1))
            self.assertTrue(np.all(J.col <= net.num_vars+constr.num_extra_vars-1))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(rowsA*self.T,net.num_vars+constr.num_extra_vars))
            self.assertEqual(A.nnz,Annz*self.T)
            self.assertTrue(np.all(A.row <= rowsA*self.T-1))
            self.assertTrue(np.all(A.col <= net.num_vars+constr.num_extra_vars-1))
            self.assertEqual(constr.num_extra_vars,rowsJ*self.T)

            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            self.assertTrue(not np.any(np.isinf(J.data)))
            self.assertTrue(not np.any(np.isnan(J.data)))
            self.assertTrue(not np.any(np.isinf(A.data)))
            self.assertTrue(not np.any(np.isnan(A.data)))

            y0 = np.random.rand(constr.num_extra_vars)

            # Ax=b check
            self.assertEqual(norm(A.data,1),rowsA*3*self.T)
            self.assertEqual(np.sum(A.data),net.get_num_switched_shunts()*self.T)

            # f check
            index = 0
            counted = {}
            for t in range(self.T):
                for i in range(net.num_branches):
                    br = net.get_branch(i)
                    for bus in [br.bus_k,br.bus_m]:
                        if (bus.index,t) not in counted:
                            for s in bus.reg_shunts:
                                self.assertEqual(bus.number,s.reg_bus.number)
                                self.assertTrue(bus.has_flags('variable','voltage magnitude'))
                                self.assertTrue(s.has_flags('variable','susceptance'))
                                fvmin = ((bus.v_mag[t]-bus.v_min_reg) - np.sqrt((bus.v_mag[t]-bus.v_min_reg)**2. + 2.*eta))*normal
                                fvmax = ((bus.v_max_reg-bus.v_mag[t]) - np.sqrt((bus.v_max_reg-bus.v_mag[t])**2. + 2.*eta))*normal
                                fbmax = ((s.b_max-s.b[t]) - np.sqrt((s.b_max-s.b[t])**2. + 2*eta))*normal
                                fbmin = ((s.b[t]-s.b_min) - np.sqrt((s.b[t]-s.b_min)**2. + 2*eta))*normal
                                self.assertLess(np.abs(fvmin-f[index]),1e-10*(1+np.abs(fvmin)))
                                self.assertLess(np.abs(fvmax-f[index+1]),1e-10*(1+np.abs(fvmax)))
                                self.assertLess(np.abs(fbmax-f[index+2]),1e-10*(1+np.abs(fbmax)))
                                self.assertLess(np.abs(fbmin-f[index+3]),1e-10*(1+np.abs(fbmin)))
                                index += 4
                        counted[(bus.index,t)] = True

            # Jacobian check
            constr.eval(x0,y0)
            f0 = constr.f.copy()
            J0 = constr.J.copy()
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                f1 = constr.f

                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if not f.size:
                    continue

                j = np.random.randint(0,f.shape[0])

                constr.eval(x0,y0)

                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular

                H0 = (H0 + H0.T) - triu(H0)

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.tocsr()[j,:].toarray().flatten()

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0,y0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,
                                  (net.num_vars+constr.num_extra_vars,
                                   net.num_vars+constr.num_extra_vars))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+constr.num_extra_vars)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.T*coeff

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sensitivities
            net.clear_sensitivities()
            for t in range(self.T):
                for i in range(net.num_buses):
                    bus = net.get_bus(i)
                    self.assertEqual(bus.sens_v_reg_by_shunt[t],0.)
            sens = np.zeros(constr.f.size)
            counter = 0
            flags = dict([((t,b.index),False) for t in range(self.T) for b in net.buses])
            for t in range(self.T):
                for branch in net.branches:
                    for bus in [branch.bus_k,branch.bus_m]:
                        if not flags[(t,bus.index)]:
                            for shunt in bus.reg_shunts:
                                sens[counter:counter+4] = bus.index*t
                                counter += 4
                        flags[(t,bus.index)] = True
            self.assertEqual(counter,constr.f.size)
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
            flags = dict([((t,b.index),False) for t in range(self.T) for b in net.buses])
            for t in range(self.T):
                for branch in net.branches:
                    for bus in [branch.bus_k,branch.bus_m]:
                        if not flags[(t,bus.index)]:
                            for shunt in bus.reg_shunts:
                                self.assertEqual(bus.sens_v_reg_by_shunt[t],bus.index*t)
                        flags[(t,bus.index)] = True

    def test_robustness(self):

        for case in test_cases.CASES:

            net = pf.Network(self.T)

            constraints = [pf.Constraint('variable nonlinear bounds',net),
                           pf.Constraint('variable fixing',net),
                           pf.Constraint('generator active power participation',net),
                           pf.Constraint('generator reactive power participation',net),
                           pf.Constraint('AC power balance',net),
                           pf.Constraint('DC power balance',net),
                           pf.Constraint('voltage regulation by generators',net),
                           pf.Constraint('voltage regulation by transformers',net),
                           pf.Constraint('voltage regulation by shunts',net),
                           pf.Constraint('AC branch flow limits',net)]

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

            list(map(lambda c: c.eval(x0),constraints))
            list(map(lambda c: c.analyze(),constraints))
            list(map(lambda c: c.eval(x0),constraints))

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
            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            constraints = [pf.Constraint('variable nonlinear bounds',net),
                           pf.Constraint('variable fixing',net),
                           pf.Constraint('generator active power participation',net),
                           pf.Constraint('generator reactive power participation',net),
                           pf.Constraint('AC power balance',net),
                           pf.Constraint('DC power balance',net),
                           pf.Constraint('voltage regulation by generators',net),
                           pf.Constraint('voltage regulation by transformers',net),
                           pf.Constraint('voltage regulation by shunts',net),
                           pf.Constraint('AC branch flow limits',net)]

            # Update network
            list(map(lambda c: c.update_network(),constraints))

            # After updating network
            list(map(lambda c: c.analyze(),constraints))
            list(map(lambda c: c.eval(x0),constraints))

            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))

            # Add variables
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
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
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,
                             (2*net.num_buses +
                              2*net.num_generators +
                              net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts()+
                              3*net.num_batteries)*self.T)

            x0 = net.get_var_values()

            # Before analyzing
            list(map(lambda c: c.clear_error(),constraints))
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            list(map(lambda c: c.clear_error(),constraints))

            # Do it right
            list(map(lambda c: c.analyze(),constraints))
            list(map(lambda c: c.eval(x0),constraints))
            for c in constraints:
                self.assertTrue(isinstance(c.b,np.ndarray))
                self.assertTrue(isinstance(c.A,coo_matrix))
                self.assertTrue(isinstance(c.f,np.ndarray))
                self.assertTrue(isinstance(c.J,coo_matrix))
                self.assertEqual(c.A.shape[1],net.num_vars+c.num_extra_vars)
                self.assertEqual(c.J.shape[1],net.num_vars+c.num_extra_vars)
                if c.f.size:
                    self.assertTupleEqual(c.get_H_single(0).shape,
                                          (net.num_vars+c.num_extra_vars,net.num_vars+c.num_extra_vars))
                else:
                    self.assertTupleEqual(c.get_H_single(0).shape,(0,0))

    def test_constr_DCPF(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len([b for b in net.buses if b.loads]))
            for b in net.buses:
                if b.loads:
                    self.assertGreater(len(b.var_generators),0)
                    for vargen in b.var_generators:
                        self.assertEqual(vargen.bus,b)

            # batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P *= -1.

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
                              2*net.num_batteries))

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('DC power balance',net)
            self.assertEqual(constr.name,'DC power balance')
            
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
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)

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
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.A_nnz,
                             (net.num_generators +
                              net.num_loads +
                              net.num_var_generators +
                              4*net.num_branches -
                              2*r +
                              2*net.get_num_phase_shifters()+
                              2*net.num_batteries))
            self.assertTupleEqual(b.shape,(net.num_buses,))
            self.assertTupleEqual(f.shape,(0,))
            self.assertTupleEqual(A.shape,(net.num_buses,net.num_vars))
            self.assertEqual(A.nnz,constr.A_nnz)
            self.assertTupleEqual(J.shape,(0,net.num_vars))

            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(A.nnz,
                             (net.num_generators +
                              net.num_loads +
                              net.num_var_generators +
                              4*net.num_branches -
                              2*r +
                              2*net.get_num_phase_shifters()+
                              2*net.num_batteries))

            # Extract pieces
            P1 = net.get_var_projection('bus','voltage angle')
            P2 = net.get_var_projection('generator','active power')
            P3 = net.get_var_projection('variable generator','active power')
            P4 = net.get_var_projection('branch','phase shift')
            P5 = net.get_var_projection('load','active power')
            P6 = net.get_var_projection('battery','charging power')

            G = A*P2.T
            R = A*P3.T
            Atheta = -A*P1.T
            Aphi = -A*P4.T
            L = -A*P5.T
            B = -A*P6.T
            x = np.random.randn(net.num_vars)
            p = P2*x
            r = P3*x
            theta = P1*x
            phi = P4*x
            l = P5*x
            Pb = P6*x
            self.assertLess(norm((G*p+R*r-Atheta*theta-Aphi*phi-L*l-B*Pb)-A*x),1e-10)

            # Sensitivities
            for bus in net.buses:
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
            new_sens = np.random.randn(net.num_buses)
            constr.store_sensitivities(new_sens,None,None,None)
            for bus in net.buses:
                self.assertNotEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
                self.assertEqual(bus.sens_P_balance,new_sens[bus.index])

            # mismatches
            mismatches = A*x0-b
            for bus in net.buses:
                mis = 0
                for gen in bus.generators:
                    mis += gen.P
                for vargen in bus.var_generators:
                    mis += vargen.P
                for load in bus.loads:
                    mis -= load.P
                for bat in bus.batteries:
                    mis -= bat.P
                for br in bus.branches_k:
                    mis -= br.P_km_DC
                for br in bus.branches_m:
                    mis += br.P_km_DC
                self.assertLess(np.abs(mismatches[bus.index]-mis),1e-8)

            # No variables
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            constr.del_matvec()
            constr.analyze()
            f1 = constr.f
            J1 = constr.J
            A1 = constr.A
            b1 = constr.b
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertTupleEqual(b1.shape,(net.num_buses,))
            self.assertTupleEqual(f1.shape,(0,))
            self.assertTupleEqual(A1.shape,(net.num_buses,net.num_vars))
            self.assertEqual(A1.nnz,constr.A_nnz)
            self.assertTupleEqual(J1.shape,(0,net.num_vars))
            x1 = net.get_var_values()
            self.assertTrue(type(x1) is np.ndarray)
            self.assertTupleEqual(x1.shape,(net.num_vars,))

            mismatches1 = A1*x1-b1
            for bus in net.buses:
                mis = 0
                for gen in bus.generators:
                    mis += gen.P
                for vargen in bus.var_generators:
                    mis += vargen.P
                for load in bus.loads:
                    mis -= load.P
                for bat in bus.batteries:
                    mis -= bat.P
                for br in bus.branches_k:
                    mis -= br.P_km_DC
                for br in bus.branches_m:
                    mis -= br.P_mk_DC
                self.assertLess(np.abs(mismatches1[bus.index]-mis),1e-8)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)

            # batteries
            for bat in net.batteries:
                bat.P = np.random.randn(self.T)*10

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
                              2*net.num_batteries)*self.T)
            x0 = net.get_var_values()

            # Count something
            r = 0
            for b in net.buses:
                if b.is_slack():
                    r += len(b.branches)

            # Constraint
            constr = pf.Constraint('DC power balance',net)
            self.assertEqual(constr.name,'DC power balance')

            # Analyze
            constr.analyze()
            A = constr.A
            b = constr.b
            self.assertEqual(constr.A_nnz,
                             (net.num_generators +
                              net.num_loads +
                              net.num_var_generators +
                              4*net.num_branches -
                              2*r +
                              2*net.get_num_phase_shifters()+
                              2*net.num_batteries)*self.T)
            self.assertTupleEqual(b.shape,(net.num_buses*self.T,))
            self.assertTupleEqual(A.shape,(net.num_buses*self.T,net.num_vars))
            self.assertEqual(A.nnz,constr.A_nnz)

            # Eval
            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(A.nnz,
                             (net.num_generators +
                              net.num_loads +
                              net.num_var_generators +
                              4*net.num_branches -
                              2*r +
                              2*net.get_num_phase_shifters()+
                              2*net.num_batteries)*self.T)

            # Mismatches
            mismatches = A*x0-b
            for t in range(self.T):
                for bus in net.buses:
                    mis = 0
                    for gen in bus.generators:
                        mis += gen.P[t]
                    for vargen in bus.var_generators:
                        mis += vargen.P[t]
                    for load in bus.loads:
                        mis -= load.P[t]
                    for bat in bus.batteries:
                        mis -= bat.P[t]
                    for br in bus.branches_k:
                        mis -= br.P_km_DC[t]
                    for br in bus.branches_m:
                        mis -= br.P_mk_DC[t]
                    self.assertLess(np.abs(mismatches[bus.index+t*net.num_buses]-mis),1e-8)

            # No variables
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            constr.del_matvec()
            constr.analyze()
            A1 = constr.A
            b1 = constr.b
            x1 = net.get_var_values()
            self.assertTupleEqual(x1.shape,(0,))

            mismatches1 = A1*x1-b1
            for t in range(self.T):
                for bus in net.buses:
                    mis = 0
                    for gen in bus.generators:
                        mis += gen.P[t]
                    for vargen in bus.var_generators:
                        mis += vargen.P[t]
                    for load in bus.loads:
                        mis -= load.P[t]
                    for bat in bus.batteries:
                        mis -= bat.P[t]
                    for br in bus.branches_k:
                        mis -= br.P_km_DC[t]
                    for br in bus.branches_m:
                        mis -= br.P_mk_DC[t]
                    self.assertLess(np.abs(mismatches1[bus.index+t*net.num_buses]-mis),1e-8)

    def test_constr_DC_FLOW_LIM(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)

            # Variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            self.assertEqual(net.num_vars,net.num_buses-net.get_num_slack_buses())

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('DC branch flow limits',net)
            self.assertEqual(constr.name,'DC branch flow limits')

            # Num constr
            num_constr = len([br for br in net.branches if br.ratingA != 0.])

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

            # Before
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,0))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,0))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.G_row,0)

            # Analyze
            constr.analyze()
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            u = constr.u
            G = constr.G
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.G_row,num_constr)

            self.assertTupleEqual(b.shape,(0,))
            self.assertTupleEqual(f.shape,(0,))
            self.assertTupleEqual(l.shape,(num_constr,))
            self.assertTupleEqual(u.shape,(num_constr,))

            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertTupleEqual(G.shape,(num_constr,net.num_vars))
            self.assertEqual(G.nnz,constr.G_nnz)

            self.assertTrue(np.all(l <= u))

            num = 0
            for br in net.branches:
                if br.ratingA == 0.:
                    continue
                if not br.bus_k.is_slack():
                    num += 1
                if not br.bus_m.is_slack():
                    num += 1
            self.assertEqual(num,constr.G_nnz)

            counter = 0
            index = 0
            for br in net.branches:
                if br.ratingA == 0.:
                    continue
                off = 0
                if br.bus_k.is_slack():
                    off = br.b*br.bus_k.v_ang
                else:
                    self.assertEqual(G.row[counter],index)
                    self.assertEqual(G.col[counter],br.bus_k.index_v_ang)
                    self.assertEqual(G.data[counter],-br.b)
                    counter += 1
                if br.bus_m.is_slack():
                    off = -br.b*br.bus_m.v_ang
                else:
                    self.assertEqual(G.row[counter],index)
                    self.assertEqual(G.col[counter],br.bus_m.index_v_ang)
                    self.assertEqual(G.data[counter],br.b)
                    counter += 1
                rating = br.ratingA
                self.assertEqual(l[index],-rating+off-br.b*br.phase)
                self.assertEqual(u[index],rating+off-br.b*br.phase)
                index += 1
            self.assertEqual(counter,G.nnz)
            self.assertEqual(index,G.shape[0])

            # Flow
            Gx0 = constr.G*x0
            self.assertTupleEqual(Gx0.shape,(num_constr,))
            index = 0
            for branch in net.branches:
                if branch.ratingA == 0.:
                    continue
                bus1 = branch.bus_k
                bus2 = branch.bus_m
                if bus1.is_slack():
                    flow = Gx0[index]-branch.b*(bus1.v_ang-branch.phase)
                elif bus2.is_slack():
                    flow = Gx0[index]-branch.b*(-bus2.v_ang-branch.phase)
                else:
                    flow = Gx0[index]-branch.b*(-branch.phase)
                self.assertLess(np.abs(branch.P_km_DC-flow),1e-10)
                index += 1

            # Sensitivities
            index = 0
            for branch in net.branches:
                self.assertEqual(branch.sens_P_u_bound,0.)
                self.assertEqual(branch.sens_P_l_bound,0.)
            mu = np.random.randn(num_constr)
            pi = np.random.randn(num_constr)
            self.assertEqual(constr.G.shape[0],num_constr)
            constr.store_sensitivities(None,None,mu,pi)
            for branch in net.branches:
                if branch.ratingA == 0.:
                    continue
                self.assertEqual(branch.sens_P_u_bound,mu[index])
                self.assertEqual(branch.sens_P_l_bound,pi[index])
                index += 1
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.G_row,num_constr)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            self.assertEqual(net.num_vars,0)

            # Nonzero angles
            for bus in net.buses:
                bus.v_ang = np.random.randn()*np.ones(self.T)

            # Variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            self.assertEqual(net.num_vars,(net.num_buses-net.get_num_slack_buses())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Num constr
            num_constr = len([br for br in net.branches if br.ratingA != 0.])

            # Constraint
            constr = pf.Constraint('DC branch flow limits',net)
            self.assertEqual(constr.name,'DC branch flow limits')
            constr.analyze()
            G = constr.G
            l = constr.l
            u = constr.u
            self.assertTupleEqual(l.shape,(num_constr*self.T,))
            self.assertTupleEqual(u.shape,(num_constr*self.T,))
            self.assertTupleEqual(G.shape,(num_constr*self.T,net.num_vars))
            Projs = []
            for t in range(self.T):
                Projs.append(net.get_var_projection('all','all',t,t))
            Gs = [G*P.T for P in Projs]
            x0s = [P*x0 for P in Projs]
            Gx0s = [(Gs[t]*x0s[t])[t*num_constr:(t+1)*num_constr] for t in range(self.T)]
            ls = [l[t*num_constr:(t+1)*num_constr] for t in range(self.T)]
            us = [u[t*num_constr:(t+1)*num_constr] for t in range(self.T)]
            for t in range(self.T):
                self.assertLessEqual(norm(Gx0s[t]-Gx0s[0]),1e-10*norm(Gx0s[0]))
                self.assertLessEqual(norm(ls[t]-ls[0]),1e-10*norm(ls[0]))
                self.assertLessEqual(norm(us[t]-us[0]),1e-10*norm(us[0]))

    def test_constr_LINPF(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # load
            if sum([l.P[0] for l in net.loads]) < 0:
                lmin = np.min([l.P for l in net.loads])
                for l in net.loads:
                    l.P = l.P + np.abs(lmin)

            # add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))
            for vargen in net.var_generators:
                vargen.Q = np.abs(vargen.P)
                for t in range(self.T):
                    self.assertGreater(vargen.Q[t],0.)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
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
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            self.assertEqual(net.num_vars,
                             (2*net.get_num_buses() +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              net.num_var_generators*2)*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('linearized AC power balance',net)
            self.assertEqual(constr.name,'linearized AC power balance')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G

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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)

            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            num_Annz = (net.get_num_buses()*4 +
                        net.get_num_branches()*8 +
                        net.get_num_tap_changers()*4 +
                        net.get_num_phase_shifters()*4 +
                        net.get_num_switched_shunts() +
                        net.get_num_slack_gens() +
                        net.get_num_reg_gens()+
                        net.num_var_generators*2)

            constr.analyze()
            self.assertEqual(constr.A_nnz,0)
            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G

            # After
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(2*net.num_buses*self.T,))
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(2*net.num_buses*self.T,net.num_vars))
            self.assertEqual(A.nnz,num_Annz*self.T)
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertTrue(not np.any(np.isinf(b)))
            self.assertTrue(not np.any(np.isnan(b)))

            # Check with ACPF
            constrPF = pf.Constraint('AC power balance',net)
            self.assertEqual(constrPF.name,'AC power balance')
            constrPF.analyze()
            constrPF.eval(x0)
            self.assertEqual(A.nnz,constrPF.J.nnz)
            self.assertTrue(np.all(A.row == constrPF.J.row))
            self.assertTrue(np.all(A.col == constrPF.J.col))
            self.assertTrue(np.all(A.data == constrPF.J.data))
            self.assertGreater(norm(A.row),0)
            self.assertGreater(norm(A.col),0)
            self.assertGreater(norm(A.data),0)
            self.assertGreater(norm(b),0)
            self.assertLess(norm(b-(constrPF.J*x0-constrPF.f)),1e-10*(norm(b)+1))

            # After eval
            constr.eval(np.zeros(x0.size))
            self.assertEqual(constr.A.nnz,constrPF.J.nnz)
            self.assertTrue(np.all(constr.A.row == constrPF.J.row))
            self.assertTrue(np.all(constr.A.col == constrPF.J.col))
            self.assertTrue(np.all(constr.A.data == constrPF.J.data))
            self.assertGreater(norm(constr.A.row),0)
            self.assertGreater(norm(constr.A.col),0)
            self.assertGreater(norm(constr.A.data),0)
            self.assertGreater(norm(constr.b),0)
            self.assertLess(norm(constr.b-(constrPF.J*x0-constrPF.f)),1e-10*(norm(b)+1))

    def test_constr_GEN_RAMP(self):

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            self.assertEqual(net.num_vars,0)

            # Gens
            for gen in net.generators:
                gen.dP_max = np.random.rand()*100.
                gen.P_prev = np.random.rand()*10.
                gen.P = np.random.rand()*20

            # Vars
            net.set_flags('generator',
                          'variable',
                          'not slack',
                          'active power')
            num = net.num_generators-net.get_num_slack_gens()
            self.assertEqual(net.num_vars,num*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('generator ramp limits',net)
            self.assertEqual(constr.name,'generator ramp limits')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

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
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            constr.analyze()
            self.assertEqual(constr.A_nnz,0)
            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(num*self.T,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(num*self.T,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(num*self.T,net.num_vars))
            self.assertEqual(G.nnz,num*(1 + (self.T-1)*2))
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            for t in range(self.T):
                for gen in net.generators:
                    if not gen.is_slack():
                        ac = np.where(G.col == gen.index_P[t])[0]

                        # Last time
                        if t == self.T-1:
                            self.assertEqual(ac.size,1)
                            i = G.row[ac[0]]
                            self.assertEqual(G.data[ac[0]],1.)
                            self.assertEqual(l[i],-gen.dP_max)
                            self.assertEqual(u[i],gen.dP_max)
                            ar = np.where(G.row == i)[0]
                            self.assertEqual(ar.size,2)
                            for j in ar:
                                if G.col[j] == gen.index_P[t]:
                                    pass
                                else:
                                    self.assertEqual(G.col[j],gen.index_P[t-1])
                                    self.assertEqual(G.data[j],-1.)

                        # Not last time
                        else:
                            self.assertEqual(ac.size,2)
                            for i in ac:
                                self.assertEqual(G.col[i],gen.index_P[t])

                                # added
                                if G.data[i] == -1.:
                                    self.assertEqual(l[G.row[i]],-gen.dP_max)
                                    self.assertEqual(u[G.row[i]],gen.dP_max)
                                    ar = np.where(G.row == G.row[i])[0]
                                    self.assertEqual(ar.size,2)
                                    for j in ar:
                                        if G.col[j] == gen.index_P[t]:
                                            pass
                                        else:
                                            self.assertEqual(G.col[j],gen.index_P[t+1])
                                            self.assertEqual(G.data[j],1.)

                                # subtracted
                                else:
                                    if t == 0:
                                        self.assertEqual(l[G.row[i]],-gen.dP_max+gen.P_prev)
                                        self.assertEqual(u[G.row[i]],gen.dP_max+gen.P_prev)
                                    else:
                                        self.assertEqual(l[G.row[i]],-gen.dP_max)
                                        self.assertEqual(u[G.row[i]],gen.dP_max)
                                        ar = np.where(G.row == G.row[i])[0]
                                        self.assertEqual(ar.size,2)
                                        for j in ar:
                                            if G.col[j] == gen.index_P[t]:
                                                pass
                                            else:
                                                self.assertEqual(G.col[j],gen.index_P[t-1])
                                                self.assertEqual(G.data[j],-1.)

    def test_constr_AC_FLOW_LIM(self):

        # Constants
        h = 1e-11
        tol = 1e-2
        eps = 1.1 # %
        param = 1e-6

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            self.assertEqual(net.num_vars,
                             (2*net.get_num_buses() +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters())*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constr
            constr = pf.Constraint('AC branch flow limits',net)
            self.assertEqual(constr.name,'AC branch flow limits')
            constr.analyze()
            num_constr = len([br for br in net.branches if br.ratingA != 0.])*2*net.num_periods
            self.assertTupleEqual(constr.f.shape,(num_constr,))
            self.assertEqual(constr.J_row,num_constr)

            # zero ratings
            for br in net.branches:
                if br.ratingA == 0.:
                    br.ratingA = 100.

            # Constraint
            constr = pf.Constraint('AC branch flow limits',net)
            self.assertEqual(constr.name,'AC branch flow limits')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

            # Before
            self.assertEqual(constr.num_extra_vars,0)
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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.G_row,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            self.assertEqual(constr.num_extra_vars,0)

            num_constr = net.get_num_branches()*2*self.T
            num_Jnnz = (net.get_num_branches()*8 +
                        net.get_num_tap_changers()*2 +
                        net.get_num_phase_shifters()*2)*self.T+num_constr

            constr.analyze()
            self.assertEqual(num_Jnnz,constr.J_nnz)
            self.assertEqual(0,constr.G_nnz)
            self.assertEqual(num_constr,constr.J_row)
           
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u
            
            # After analyze
            self.assertEqual(constr.num_extra_vars,num_constr)
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(num_constr,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(num_constr,net.num_vars+num_constr))
            self.assertEqual(J.nnz,num_Jnnz)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars+num_constr))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(num_constr,net.num_vars+num_constr))
            self.assertEqual(G.nnz,num_constr)
            self.assertTrue(np.all(G.row == np.array(range(num_constr))))
            self.assertTrue(np.all(G.col == np.array(range(net.num_vars,net.num_vars+num_constr))))
            self.assertTrue(np.all(G.row == G.col-net.num_vars))
            self.assertTrue(np.all(G.data == 1.))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(num_constr,))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(num_constr,))
            self.assertTrue(np.all(l == 0.))
            for t in range(net.num_periods):
                for branch in net.branches:
                    i = t*net.num_branches*2+2*branch.index
                    self.assertEqual(u[i],branch.ratingA)
                    self.assertEqual(u[i+1],branch.ratingA)

            # Hessian structure
            for i in range(constr.J.shape[0]):
                H = constr.get_H_single(i)
                self.assertTupleEqual(H.shape,(net.num_vars+num_constr,net.num_vars+num_constr))
                self.assertTrue(np.all(H.row >= H.col))
            Hcomb = constr.H_combined
            H_comb_nnz = 2*(net.num_branches*10 +
                            net.get_num_tap_changers()*5+
                            net.get_num_phase_shifters()*5)*self.T
            self.assertTupleEqual(Hcomb.shape,(net.num_vars+num_constr,net.num_vars+num_constr))
            self.assertTrue(np.all(Hcomb.row >= Hcomb.col))
            self.assertEqual(Hcomb.nnz,H_comb_nnz)

            y0 = np.random.randn(num_constr)

            constr.eval(x0,y0)
            self.assertEqual(num_constr,constr.J_row)
            self.assertEqual(0,constr.G_nnz)
            self.assertEqual(num_Jnnz,constr.J_nnz)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u
            constr.combine_H(np.ones(f.size),False)
            Hcomb = constr.H_combined

            # After eval
            self.assertTrue(not np.any(np.isinf(f)))
            self.assertTrue(not np.any(np.isnan(f)))
            
            # Cross check current magnitudes
            J_row = 0
            for t in range(net.num_periods):
                for branch in net.branches:
                    i = t*net.num_branches*2+2*branch.index
                    Pkm = branch.get_P_km()[t]
                    Qkm = branch.get_Q_km()[t]
                    Pmk = branch.get_P_mk()[t]
                    Qmk = branch.get_Q_mk()[t]
                    vk = branch.bus_k.v_mag[t]
                    vm = branch.bus_m.v_mag[t]
                    ikmmag = np.sqrt(np.abs((Pkm/vk) + 1j*(Qkm/vk))**2.+param)
                    imkmag = np.sqrt(np.abs((Pmk/vm) + 1j*(Qmk/vm))**2.+param)
                    error_km = 100.*np.abs(ikmmag-f[i]-y0[J_row])/max([ikmmag,tol])
                    error_mk = 100.*np.abs(imkmag-f[i+1]-y0[J_row+1])/max([imkmag,tol])
                    self.assertLess(error_km,eps)
                    self.assertLess(error_mk,eps)
                    J_row += 2

            # Jacobian check 1
            constr.eval(x0,y0)
            f0 = constr.f.copy()
            J0 = constr.J.copy()
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+num_constr)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                f1 = constr.f

                Jd_exact = J0*d
                Jd_approx = (f1-f0)/h
                error = 100.*norm(Jd_exact-Jd_approx)/np.maximum(norm(Jd_exact),tol)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                j = np.random.randint(0,f.shape[0])

                constr.eval(x0,y0)

                g0 = constr.J.tocsr()[j,:].toarray().flatten()
                H0 = constr.get_H_single(j)

                self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
                
                H0 = (H0 + H0.T) - triu(H0)

                d = np.random.randn(net.num_vars+num_constr)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)
                
                g1 = constr.J.tocsr()[j,:].toarray().flatten()

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),tol)
                self.assertLessEqual(error,EPS)

            # Combined Hessian check
            h = 1e-12
            coeff = np.random.randn(f0.shape[0])
            constr.eval(x0,y0)
            constr.combine_H(coeff,False)
            J0 = constr.J
            g0 = J0.T*coeff
            H0 = constr.H_combined.copy()
            self.assertTrue(type(H0) is coo_matrix)
            self.assertTupleEqual(H0.shape,(net.num_vars+num_constr,net.num_vars+num_constr))
            self.assertTrue(np.all(H0.row >= H0.col)) # lower triangular
            H0 = (H0 + H0.T) - triu(H0)
            for i in range(NUM_TRIALS):

                d = np.random.randn(net.num_vars+num_constr)

                x = x0 + h*d[:net.num_vars]
                y = y0 + h*d[net.num_vars:]

                constr.eval(x,y)

                g1 = constr.J.T*coeff

                Hd_exact = H0*d
                Hd_approx = (g1-g0)/h
                error = 100.*norm(Hd_exact-Hd_approx)/np.maximum(norm(Hd_exact),tol)
                self.assertLessEqual(error,EPS)

            # Jacobian check 2
            if os.environ.get('TEST_SLOW'):
                h = 1e-12
                tol = 1e0
                eps = 5.0 # %
                constr.eval(x0,y0)
                f0 = constr.f.copy()
                J0 = constr.J.copy()
                self.assertEqual(constr.num_extra_vars,num_constr)
                def check_J(index,nnz=None):
                    e = np.zeros(net.num_vars+constr.num_extra_vars)
                    e[index] = 1.
                    x = x0 + h*e[:net.num_vars]
                    y = y0 + h*e[net.num_vars:]
                    constr.eval(x,y)
                    f1 = constr.f
                    Je_exact = J0*e
                    Je_approx = (f1-f0)/h
                    error = 100.*norm(Je_exact-Je_approx)/np.maximum(norm(Je_exact),tol)
                    self.assertLessEqual(np.sum((f1-f0) != 0.),nnz)
                    try:
                        self.assertLessEqual(error,eps)
                        return True
                    except AssertionError:
                        return False
                for t in range(net.num_periods):
                    for branch in net.branches:
                        if branch.is_tap_changer():
                            self.assertTrue(check_J(branch.index_ratio[t],2))
                    for branch in net.branches:
                        if branch.is_phase_shifter():
                            self.assertTrue(check_J(branch.index_phase[t],2))
                    num_bad = 0
                    for bus in net.buses:
                        if not check_J(bus.index_v_mag[t],2*bus.degree):
                            num_bad += 1
                    self.assertLess((100.*num_bad)/net.num_buses,1.) # less then 1 %
                    num_bad = 0
                    for bus in net.buses:
                        if not check_J(bus.index_v_ang[t],2*bus.degree):
                            num_bad += 1
                    self.assertLess((100.*num_bad)/net.num_buses,1.) # less then 1 %

    def test_constr_DUMMY(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len([b for b in net.buses if b.loads]))
            for b in net.buses:
                if b.loads:
                    self.assertGreater(len(b.var_generators),0)
                    for vargen in b.var_generators:
                        self.assertEqual(vargen.bus,b)

            # batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P *= -1.

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

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Ref constraint
            constrREF = pf.Constraint('DC power balance',net)
            self.assertEqual(constrREF.name,'DC power balance')

            # Dummy constraint
            constr = pf.constraints.DummyDCPF(net)
            self.assertEqual(constr.name,'dummy DC power balance')
            
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.A_row,constrREF.A_row)
            self.assertEqual(constr.A_nnz,constrREF.A_nnz)
            
            self.assertEqual(constr.b.size,0)
            self.assertEqual(constr.A.shape[0],0)
            self.assertEqual(constr.A.shape[1],0)
            self.assertEqual(constr.A.nnz,0)
           
            constrREF.analyze()
            constr.analyze()

            self.assertEqual(constr.A_row,0)
            self.assertGreater(constr.A_nnz,0)
            self.assertEqual(constr.A_row,constrREF.A_row)
            self.assertEqual(constr.A_nnz,constrREF.A_nnz)
            
            self.assertTrue(np.all(constr.b == constrREF.b))
            self.assertTrue(np.all(constr.A.row == constrREF.A.row))
            self.assertTrue(np.all(constr.A.col == constrREF.A.col))
            self.assertTrue(np.all(constr.A.data == constrREF.A.data))
            
            self.assertTupleEqual(constr.l.shape,(0,))
            self.assertTupleEqual(constr.u.shape,(0,))
            self.assertTupleEqual(constr.f.shape,(0,))
            self.assertTupleEqual(constr.G.shape,(0,net.num_vars))
            self.assertTupleEqual(constr.J.shape,(0,net.num_vars))

            constrREF.eval(net.get_var_values())
            constr.eval(net.get_var_values())

            self.assertTrue(np.all(constr.b == constrREF.b))
            self.assertTrue(np.all(constr.A.row == constrREF.A.row))
            self.assertTrue(np.all(constr.A.col == constrREF.A.col))
            self.assertTrue(np.all(constr.A.data == constrREF.A.data))

    def test_constr_BAT_DYN(self):

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,5)
            self.assertEqual(net.num_periods,5)
            self.assertEqual(net.num_vars,0)

            # Add battries
            gen_buses = net.get_generator_buses()
            net.add_batteries(gen_buses,20.,40.,0.8,0.7)
            self.assertEqual(net.num_batteries,len(gen_buses))
            self.assertGreater(net.num_batteries,0)

            # Vars
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,5*3*net.num_batteries)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('battery dynamics',net)
            self.assertEqual(constr.name,'battery dynamics')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

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
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            constr.analyze()
            self.assertEqual(constr.A_row,(5+1)*net.num_batteries)
            self.assertEqual(constr.A_nnz,5*4*net.num_batteries)
            self.assertEqual(constr.G_nnz,0)
            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(6*net.num_batteries,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(6*net.num_batteries,net.num_vars))
            self.assertEqual(A.nnz,5*4*net.num_batteries)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            
            for t in range(5):
                for bat in net.batteries:
                    self.assertTrue(bat.has_flags('variable',['charging power','energy level']))
                    
                    aPc = np.where(A.col == bat.index_Pc[t])[0]
                    aPd = np.where(A.col == bat.index_Pd[t])[0]
                    aE = np.where(A.col == bat.index_E[t])[0]
                    if t < 5-1:
                        aEE = np.where(A.col == bat.index_E[t+1])[0]
                    
                    self.assertEqual(aPc.size,1)
                    self.assertEqual(aPd.size,1)
                   
                    eq_row = A.row[aPc[0]]
                    self.assertEqual(eq_row,A.row[aPd[0]])
                    self.assertEqual(A.data[aPc[0]],-bat.eta_c)
                    self.assertEqual(A.data[aPd[0]],1./bat.eta_d)
 
                    if t == 0:
                        self.assertEqual(aE.size,2)
                        
                        # init eq
                        j = aE[0] 
                        self.assertEqual(A.data[j],1.)
                        self.assertEqual(b[A.row[j]],bat.E_init)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,1)
                        
                        # update eq E_{t+1} - E_t - eta_c Pc_t + (1/eta_d) Pd_t = 0
                        j = aE[1]
                        self.assertEqual(A.data[j],-1.)
                        self.assertEqual(b[A.row[j]],0.)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,4)
                        self.assertEqual(A.row[j],eq_row)
                        self.assertEqual(A.row[j],A.row[aEE[0]])
                        
                    elif t < 5-1:
                        self.assertEqual(aE.size,2)

                        # update eq E_t - E_{t-1} - eta_c Pc_{t-1} + (1/eta_d) Pd_{t-1} = 0
                        j = aE[0]
                        self.assertEqual(A.data[j],1.)
                        self.assertEqual(b[A.row[j]],0.)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,4)
                        self.assertNotEqual(A.row[j],eq_row)
                        self.assertNotEqual(A.row[j],A.row[aEE[0]])

                        # update eq E_{t+1} - E_t - eta_c Pc_t + (1/eta_d) Pd_t = 0
                        j = aE[1]
                        self.assertEqual(A.data[j],-1.)
                        self.assertEqual(b[A.row[j]],0.)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,4)
                        self.assertEqual(A.row[j],eq_row)
                        self.assertEqual(A.row[j],A.row[aEE[0]])

                    else:
                        self.assertEqual(aE.size,2)

                        # update eq E_t - E_{t-1} - eta_c Pc_{t-1} + (1/eta_d) Pd_{t-1} = 0
                        j = aE[0]
                        self.assertEqual(A.data[j],1.)
                        self.assertEqual(b[A.row[j]],0.)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,4)
                        self.assertNotEqual(A.row[j],eq_row)

                        # update eq - E_t - eta_c Pc_t + (1/eta_d) Pd_t = -E_final
                        j = aE[1]
                        self.assertEqual(A.data[j],-1.)
                        self.assertEqual(b[A.row[j]],-bat.E_final)
                        self.assertEqual(np.where(A.row == A.row[j])[0].size,3)
                        self.assertEqual(A.row[j],eq_row)

    def test_constr_LOAD_PF(self):

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            self.assertEqual(net.num_vars,0)

            # Powers
            for load in net.loads:
                load.P = np.random.rand(net.num_periods)

            # Target power factors
            for load in net.loads:
                load.target_power_factor = np.random.rand()
                self.assertTrue(0 < load.target_power_factor < 1.)

            # Vars
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            self.assertEqual(net.num_vars,2*net.num_loads*self.T)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('load constant power factor',net)
            self.assertEqual(constr.name,'load constant power factor')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

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
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            constr.analyze()
            self.assertEqual(constr.A_nnz,2*net.num_loads*self.T)
            self.assertEqual(constr.A_row,net.num_loads*self.T)
            constr.eval(x0)
            self.assertEqual(constr.A_nnz,0)

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            G = constr.G
            u = constr.u

            # After
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(net.num_loads*self.T,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(net.num_loads*self.T,net.num_vars))
            self.assertEqual(A.nnz,2*net.num_loads*self.T)
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)

            for load in net.loads:
                for t in range(net.num_periods):
                    indices = np.where(A.col == load.index_P[t])[0]
                    self.assertEqual(indices.size,1)
                    row = A.row[indices[0]]
                    indices = np.where(A.row == row)[0]
                    self.assertEqual(indices.size,2)
                    for i in indices:
                        if A.col[i] == load.index_P[t]:
                            gamma = load.target_power_factor
                            factor = np.sqrt(1.-gamma**2.)/gamma
                            load.Q = load.P*factor
                            self.assertAlmostEqual(A.data[i],-factor)
                        else:
                            self.assertEqual(A.col[i],load.index_Q[t])
                            self.assertEqual(A.data[i],1.)
            
            x = net.get_var_values()
            self.assertLess(np.linalg.norm(constr.A*x-constr.b),1e-10)
            
            for load in net.loads:
                for t in range(net.num_periods):
                    self.assertAlmostEqual(load.power_factor[t],load.target_power_factor)

    def test_constr_AC_LIN_FLOW_LIM(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            self.assertEqual(net.num_vars,
                             (2*net.get_num_buses()-net.get_num_slack_buses() +
                              net.get_num_tap_changers() +
                              net.get_num_phase_shifters())*self.T)

            # Zero ratings
            for br in net.branches:
                if br.ratingA == 0.:
                    br.ratingA = 100.

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))

            # Constraint
            constr = pf.Constraint('linearized AC branch flow limits',net)
            self.assertEqual(constr.name,'linearized AC branch flow limits')

            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u

            # Before
            self.assertEqual(constr.num_extra_vars,0)
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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,0))
            self.assertEqual(G.nnz,0)
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(0,))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(0,))
            self.assertEqual(constr.J_row,0)
            self.assertEqual(constr.A_row,0)
            self.assertEqual(constr.G_row,0)
            self.assertEqual(constr.J_nnz,0)
            self.assertEqual(constr.A_nnz,0)
            self.assertEqual(constr.G_nnz,0)
            self.assertEqual(constr.num_extra_vars,0)
      
            # Tap ratios and phase shifts
            if net.get_num_tap_changers()+net.get_num_phase_shifters() > 0:
                self.assertRaises(pf.ConstraintError,constr.analyze)
                constr.clear_error()
                continue

            # No voltage magnitude bounds
            self.assertRaises(pf.ConstraintError,constr.analyze)
            self.assertRaisesRegexp(pf.ConstraintError,
                                    "AC_LIN_FLOW_LIM constraint requires variable voltage magnitudes to be bounded",
                                    constr.analyze)
            constr.clear_error()

            net.set_flags('bus',
                          'bounded',
                          'any',
                          'voltage magnitude')

            self.assertEqual(net.num_bounded,net.num_buses*self.T)

            # Library unavailable
            if not pf.info['line_flow']:
                self.assertRaisesRegexp(pf.ConstraintError,
                                        "LINE_FLOW library not available",
                                        constr.analyze)
                constr.clear_error()
                continue
        
            constr.analyze()

            self.assertGreaterEqual(constr.G_nnz,constr.G_row)
           
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
            l = constr.l
            u = constr.u
            
            # After analyze
            self.assertEqual(constr.num_extra_vars,0)
            self.assertTrue(type(f) is np.ndarray)
            self.assertTupleEqual(f.shape,(0,))
            self.assertTrue(type(b) is np.ndarray)
            self.assertTupleEqual(b.shape,(0,))
            self.assertTrue(type(J) is coo_matrix)
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertEqual(J.nnz,0)
            self.assertTrue(type(A) is coo_matrix)
            self.assertTupleEqual(A.shape,(0,net.num_vars))
            self.assertEqual(A.nnz,0)
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(constr.G_row,net.num_vars))
            self.assertFalse(np.any(np.isnan(G.data)))
            self.assertTrue(type(u) is np.ndarray)
            self.assertTupleEqual(u.shape,(constr.G_row,))
            self.assertFalse(np.any(np.isnan(u)))
            self.assertTrue(type(l) is np.ndarray)
            self.assertTupleEqual(l.shape,(constr.G_row,))
            self.assertTrue(np.all(l == -1e8))

    def tearDown(self):

        pass
