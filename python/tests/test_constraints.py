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
from scipy.sparse import coo_matrix,triu,tril,eye

NUM_TRIALS = 25
EPS = 2e0 # %
TOL = 1e-4

class TestConstraints(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

        # Random
        np.random.seed(0)

    def test_constr_FIX(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            # add vargens
            net.add_vargens(net.get_load_buses(),50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = vargen.index*1.5
                vargen.Q = vargen.index*2.5
            self.assertGreater(net.num_vargens,0)

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
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             2*net.num_buses +
                             net.get_num_slack_gens() + 
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts() +
                             net.num_vargens*2+
                             2*net.num_bats)
            
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
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_FIXED,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_FIXED,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            self.assertGreater(net.num_fixed,0)
            self.assertEqual(net.num_fixed,
                             2*(net.get_num_slack_buses()) +
                             (net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts() +
                             net.num_vargens*2+
                             2*net.num_bats)
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            constr = pf.Constraint(pf.CONSTR_TYPE_FIX,net)

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
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)

            Acounter = net.num_fixed+net.get_num_buses_reg_by_gen()
            constr.analyze()
            self.assertEqual(Acounter,constr.Acounter)
            constr.eval(x0)
            self.assertEqual(0,constr.Acounter)

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
                ar = np.where(A.col == bat.index_P)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],bat.index_P)
                self.assertEqual(b[A.row[ar[0]]],bat.P)
            for bat in net.batteries:
                ar = np.where(A.col == bat.index_E)[0]
                self.assertEqual(ar.size,1)
                self.assertEqual(A.col[ar[0]],bat.index_E)
                self.assertEqual(b[A.row[ar[0]]],bat.E)

    def test_constr_LBOUND(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            # add vargens
            net.add_vargens(net.get_load_buses(),50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = vargen.index*1.5
                vargen.Q = vargen.index*2.5
            self.assertGreater(net.num_vargens,0)
            
            self.assertEqual(net.num_bounded,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)
                load.P_max = 3.3*(load.index+1)

            # Vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_GEN,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG,pf.BUS_VAR_VDEV,pf.BUS_VAR_VVIO])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          [pf.GEN_VAR_P,pf.GEN_VAR_Q])
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_P_ADJUST,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          [pf.BRANCH_VAR_RATIO,pf.BRANCH_VAR_RATIO_DEV])
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          [pf.SHUNT_VAR_SUSC,pf.SHUNT_VAR_SUSC_DEV])
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            num_vars_saved = net.num_vars
            self.assertGreater(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_vars,
                             net.get_num_buses_reg_by_gen()*6 +
                             net.get_num_reg_gens()*2 +
                             net.get_num_P_adjust_loads() + 
                             net.get_num_tap_changers()*3 +
                             net.get_num_phase_shifters()*1 +
                             net.get_num_switched_shunts()*3 +
                             net.num_vargens*2+
                             2*net.num_bats)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            constr = pf.Constraint(pf.CONSTR_TYPE_LBOUND,net)

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
            
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)

            constr.analyze()
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)
            constr.eval(x0)
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)

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
            self.assertGreater(np.linalg.norm(G.data,np.inf),0.5)
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
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,
                                                  pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO))
                    self.assertEqual(u[bus.index_v_mag],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_v_ang],pf.BUS_INF_V_ANG)
                    self.assertEqual(u[bus.index_y],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_z],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_vl],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_vh],pf.BUS_INF_V_MAG)
                    self.assertEqual(l[bus.index_v_mag],0.)
                    self.assertEqual(l[bus.index_v_ang],-pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_y],0.)
                    self.assertEqual(l[bus.index_z],0.)
                    self.assertEqual(l[bus.index_vl],0.)
                    self.assertEqual(l[bus.index_vh],0.)
                else:
                    self.assertFalse(bus.has_flags(pf.FLAG_VARS,
                                                   pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO))
            
            for branch in net.branches:
                if branch.is_tap_changer():
                    self.assertTrue(branch.has_flags(pf.FLAG_VARS,
                                                     pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV))
                    self.assertEqual(u[branch.index_ratio],pf.BRANCH_INF_RATIO)
                    self.assertEqual(u[branch.index_ratio_y],pf.BRANCH_INF_RATIO)
                    self.assertEqual(u[branch.index_ratio_z],pf.BRANCH_INF_RATIO)
                    self.assertEqual(l[branch.index_ratio],0.)
                    self.assertEqual(l[branch.index_ratio_y],0.)
                    self.assertEqual(l[branch.index_ratio_z],0.)
                else:
                    self.assertFalse(branch.has_flags(pf.FLAG_VARS,
                                                      pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV))
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags(pf.FLAG_VARS,
                                                     pf.BRANCH_VAR_PHASE))
                    self.assertLess(np.abs(u[branch.index_phase]-np.pi*2.),1e-10)
                    self.assertLess(np.abs(l[branch.index_phase]+np.pi*2.),1e-10)
                else:
                    self.assertFalse(branch.has_flags(pf.FLAG_VARS,
                                                      pf.BRANCH_VAR_PHASE))
            
            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags(pf.FLAG_VARS,
                                                  pf.GEN_VAR_P|pf.GEN_VAR_Q))
                    self.assertEqual(u[gen.index_P],pf.GEN_INF_P)
                    self.assertEqual(u[gen.index_Q],pf.GEN_INF_Q)
                    self.assertEqual(l[gen.index_P],-pf.GEN_INF_P)
                    self.assertEqual(l[gen.index_Q],-pf.GEN_INF_Q)
                else:
                    self.assertFalse(gen.has_flags(pf.FLAG_VARS,
                                                   pf.GEN_VAR_P|pf.GEN_VAR_Q))

            for load in net.loads:
                self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                self.assertEqual(u[load.index_P],pf.LOAD_INF_P)
                self.assertEqual(l[load.index_P],-pf.LOAD_INF_P)

            for vargen in net.var_generators:
                self.assertTrue(vargen.has_flags(pf.FLAG_VARS,
                                                 pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q))
                self.assertEqual(u[vargen.index_P],pf.VARGEN_INF_P)
                self.assertEqual(u[vargen.index_Q],pf.VARGEN_INF_Q)
                self.assertEqual(l[vargen.index_P],-pf.VARGEN_INF_P)
                self.assertEqual(l[vargen.index_Q],-pf.VARGEN_INF_Q)

            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags(pf.FLAG_VARS,
                                                     pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV))
                    self.assertEqual(u[shunt.index_b],pf.SHUNT_INF_SUSC)
                    self.assertEqual(u[shunt.index_y],pf.SHUNT_INF_SUSC)
                    self.assertEqual(u[shunt.index_z],pf.SHUNT_INF_SUSC)
                    self.assertEqual(l[shunt.index_b],-pf.SHUNT_INF_SUSC)
                    self.assertEqual(l[shunt.index_y],0.)
                    self.assertEqual(l[shunt.index_z],0.)
                else:
                    self.assertFalse(shunt.has_flags(pf.FLAG_VARS,
                                                     pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV))              

            for bat in net.batteries:
                self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P))
                self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_E))
                self.assertEqual(u[bat.index_P],pf.BAT_INF_P)
                self.assertEqual(l[bat.index_P],-pf.BAT_INF_P)
                self.assertEqual(u[bat.index_E],pf.BAT_INF_E)
                self.assertEqual(l[bat.index_E],0.)
                    
            # Add bounded flags
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_REG_BY_GEN,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG,pf.BUS_VAR_VDEV,pf.BUS_VAR_VVIO])
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_BOUNDED,
                          pf.GEN_PROP_REG,
                          [pf.GEN_VAR_P,pf.GEN_VAR_Q])
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_BOUNDED,
                          pf.LOAD_PROP_P_ADJUST,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          [pf.BRANCH_VAR_RATIO,pf.BRANCH_VAR_RATIO_DEV])
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_BOUNDED,
                          pf.SHUNT_PROP_SWITCHED_V,
                          [pf.SHUNT_VAR_SUSC,pf.SHUNT_VAR_SUSC_DEV])
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_BOUNDED,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_BOUNDED,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            self.assertEqual(net.num_vars,num_vars_saved)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,net.num_vars)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            constr = pf.Constraint(pf.CONSTR_TYPE_LBOUND,net)
                        
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
            self.assertGreater(np.linalg.norm(G.data,np.inf),0.5)
            self.assertEqual(E.nnz,0)

            # Bounds
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(bus.has_flags(pf.FLAG_BOUNDED,
                                                  pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO))
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,
                                                  pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO))
                    self.assertEqual(u[bus.index_v_mag],bus.v_max)
                    self.assertEqual(u[bus.index_v_ang],pf.BUS_INF_V_ANG)
                    self.assertEqual(u[bus.index_y],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_z],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_vl],pf.BUS_INF_V_MAG)
                    self.assertEqual(u[bus.index_vh],pf.BUS_INF_V_MAG)
                    self.assertEqual(l[bus.index_v_mag],bus.v_min)
                    self.assertEqual(l[bus.index_v_ang],-pf.BUS_INF_V_ANG)
                    self.assertEqual(l[bus.index_y],0.)
                    self.assertEqual(l[bus.index_z],0.)
                    self.assertEqual(l[bus.index_vl],0.)
                    self.assertEqual(l[bus.index_vh],0.)
                else:
                    self.assertFalse(bus.has_flags(pf.FLAG_BOUNDED,
                                                   pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO))
            
            for branch in net.branches:
                if branch.is_tap_changer():
                    self.assertTrue(branch.has_flags(pf.FLAG_BOUNDED,
                                                     pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV))
                    self.assertEqual(u[branch.index_ratio],branch.ratio_max)
                    self.assertEqual(u[branch.index_ratio_y],pf.BRANCH_INF_RATIO)
                    self.assertEqual(u[branch.index_ratio_z],pf.BRANCH_INF_RATIO)
                    self.assertEqual(l[branch.index_ratio],branch.ratio_min)
                    self.assertEqual(l[branch.index_ratio_y],0.)
                    self.assertEqual(l[branch.index_ratio_z],0.)
                else:
                    self.assertFalse(branch.has_flags(pf.FLAG_BOUNDED,
                                                      pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV))
                if branch.is_phase_shifter():
                    self.assertTrue(branch.has_flags(pf.FLAG_BOUNDED,
                                                     pf.BRANCH_VAR_PHASE))
                    self.assertEqual(u[branch.index_phase],branch.phase_max)
                    self.assertEqual(l[branch.index_phase],branch.phase_min)
                else:
                    self.assertFalse(branch.has_flags(pf.FLAG_BOUNDED,
                                                      pf.BRANCH_VAR_PHASE))
            
            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags(pf.FLAG_BOUNDED,
                                                  pf.GEN_VAR_P|pf.GEN_VAR_Q))
                    self.assertEqual(u[gen.index_P],gen.P_max)
                    self.assertEqual(u[gen.index_Q],gen.Q_max)
                    self.assertEqual(l[gen.index_P],gen.P_min)
                    self.assertEqual(l[gen.index_Q],gen.Q_min)
                else:
                    self.assertFalse(gen.has_flags(pf.FLAG_BOUNDED,
                                                   pf.GEN_VAR_P|pf.GEN_VAR_Q))

            for load in net.loads:
                self.assertTrue(load.has_flags(pf.FLAG_BOUNDED,pf.LOAD_VAR_P))
                self.assertEqual(u[load.index_P],load.P_max)
                self.assertEqual(l[load.index_P],load.P_min)

            for vargen in net.var_generators:
                self.assertTrue(vargen.has_flags(pf.FLAG_BOUNDED,
                                                 pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q))
                self.assertEqual(u[vargen.index_P],vargen.P_max)
                self.assertEqual(u[vargen.index_Q],vargen.Q_max)
                self.assertEqual(l[vargen.index_P],vargen.P_min)
                self.assertEqual(l[vargen.index_Q],vargen.Q_min)

            for shunt in net.shunts:
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags(pf.FLAG_BOUNDED,
                                                     pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV))
                    self.assertEqual(u[shunt.index_b],shunt.b_max)
                    self.assertEqual(u[shunt.index_y],pf.SHUNT_INF_SUSC)
                    self.assertEqual(u[shunt.index_z],pf.SHUNT_INF_SUSC)
                    self.assertEqual(l[shunt.index_b],shunt.b_min)
                    self.assertEqual(l[shunt.index_y],0.)
                    self.assertEqual(l[shunt.index_z],0.)
                else:
                    self.assertFalse(shunt.has_flags(pf.FLAG_BOUNDED,
                                                     pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV))

            for bat in net.batteries:
                self.assertTrue(bat.has_flags(pf.FLAG_BOUNDED,pf.BAT_VAR_P))
                self.assertTrue(bat.has_flags(pf.FLAG_BOUNDED,pf.BAT_VAR_E))
                self.assertEqual(u[bat.index_P],bat.P_max)
                self.assertEqual(l[bat.index_P],bat.P_min)
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
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG))
                    self.assertNotEqual(bus.sens_v_ang_u_bound,0.)
                    self.assertNotEqual(bus.sens_v_ang_l_bound,0.)
                    self.assertEqual(bus.sens_v_ang_u_bound,mu[bus.index_v_ang])
                    self.assertEqual(bus.sens_v_ang_l_bound,pi[bus.index_v_ang])
                else:
                    self.assertEqual(bus.sens_v_ang_u_bound,0.)
                    self.assertEqual(bus.sens_v_ang_l_bound,0.)
            for gen in net.generators:
                if gen.is_regulator():
                    self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                    self.assertNotEqual(gen.sens_P_u_bound,0.)
                    self.assertNotEqual(gen.sens_P_l_bound,0.)
                    self.assertEqual(gen.sens_P_u_bound,mu[gen.index_P])
                    self.assertEqual(gen.sens_P_l_bound,pi[gen.index_P])
                else:
                    self.assertEqual(gen.sens_P_u_bound,0.)
                    self.assertEqual(gen.sens_P_l_bound,0.)
            for load in net.loads:
                self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                self.assertNotEqual(load.sens_P_u_bound,0.)
                self.assertNotEqual(load.sens_P_l_bound,0.)
                self.assertEqual(load.sens_P_u_bound,mu[load.index_P])
                self.assertEqual(load.sens_P_l_bound,pi[load.index_P])

    def test_constr_PAR_GEN_P(self):
        
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
            constr = pf.Constraint(pf.CONSTR_TYPE_PAR_GEN_P,net)

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
                    if bus.number in counted:
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
            self.assertGreater(np.linalg.norm(x),0)
            self.assertTrue(np.linalg.norm(A*x-b) < 1e-10)

    def test_constr_PAR_GEN_Q(self):
        
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
            constr = pf.Constraint(pf.CONSTR_TYPE_PAR_GEN_Q,net)

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
                    if bus.number in counted:
                        continue
                    counted[bus.number] = True
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
                if bus.is_regulated_by_gen():
                    self.assertGreater(len(bus.reg_gens),0)
                    for g in bus.reg_gens:
                        self.assertTrue(g.has_flags(pf.FLAG_VARS,pf.GEN_VAR_Q))
                        x[g.index_Q] = (g.Q_max+g.Q_min)/2.
            self.assertTrue(np.linalg.norm(A*x-b) < 1e-10)

    def test_constr_PF(self):
        
        # Constants
        h = 1e-10
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            # add vargens
            load_buses = net.get_load_buses()
            net.add_vargens(load_buses,50.,30.,5,0.05)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len(load_buses))
            for vargen in net.var_generators:
                vargen.Q = np.abs(vargen.P)
                self.assertGreater(vargen.Q,0.)

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
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            self.assertEqual(net.num_vars,
                             2*net.get_num_buses() +
                             net.get_num_slack_gens() +
                             net.get_num_reg_gens() +
                             net.get_num_tap_changers() +
                             net.get_num_phase_shifters() +
                             net.get_num_switched_shunts() +
                             net.num_vargens*2)

            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_PF,net)

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

            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)

            num_Jnnz = (net.get_num_buses()*4 +
                        net.get_num_branches()*8 +
                        net.get_num_tap_changers()*4 +
                        net.get_num_phase_shifters()*4 +
                        net.get_num_switched_shunts() +
                        net.get_num_slack_gens() +
                        net.get_num_reg_gens()+
                        net.num_vargens*2)
            
            constr.analyze()
            self.assertEqual(num_Jnnz,constr.Jcounter)
            constr.eval(x0)
            self.assertEqual(num_Jnnz,constr.Jcounter)
            
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            G = constr.G
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
            self.assertTrue(type(G) is coo_matrix)
            self.assertTupleEqual(G.shape,(0,net.num_vars))
            self.assertEqual(G.nnz,0)
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
            
            # Remove vargen injections
            fsaved = f.copy()
            x0saved = x0.copy()
            for vargen in net.var_generators:
                x0[vargen.index_P] = 0.
                x0[vargen.index_Q] = 0.
            constr.eval(x0)
            for vargen in net.var_generators:
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       constr.f[vargen.bus.index_P]-vargen.P),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       constr.f[vargen.bus.index_Q]-vargen.Q),1e-10)
            for vargen in net.var_generators:
                vargen.bus.loads[0].P = vargen.bus.loads[0].P-vargen.P
                vargen.bus.loads[0].Q = vargen.bus.loads[0].Q-vargen.Q
                pass
            constr.eval(x0)
            for vargen in net.var_generators:
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       constr.f[vargen.bus.index_P]),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       constr.f[vargen.bus.index_Q]),1e-10)

            # Jacobian check
            x0 = x0saved.copy()
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            constr.store_sensitivities(None,sens,None,None)
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
            
            Jnnz = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_gen() and not bus.is_slack():
                    Jnnz += 2 + 2*len(bus.reg_gens)
                    
            Annz = 3*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            
            rowsJ = 2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            rowsA = net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()
                        
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
                self.assertLessEqual(error,EPS)

            # Sigle Hessian check
            for i in range(NUM_TRIALS):

                if f.shape[0] == 0:
                    break

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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
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
            constr.store_sensitivities(None,np.zeros(constr.J.shape[0]),None,None)
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            constr.store_sensitivities(None,sens,None,None)
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
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
                error = 100.*np.linalg.norm(Jd_exact-Jd_approx)/np.maximum(np.linalg.norm(Jd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
                error = 100.*np.linalg.norm(Hd_exact-Hd_approx)/np.maximum(np.linalg.norm(Hd_exact),TOL)
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
            constr.store_sensitivities(np.zeros(constr.A.shape[0]),sens,None,None)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_shunt():
                    self.assertEqual(bus.sens_v_reg_by_shunt,-bus.index-max([s.index for s in bus.reg_shunts]))

    def test_robustness(self):

        for case in test_cases.CASES:

            net = pf.Network()

            constraints = [pf.Constraint(pf.CONSTR_TYPE_BOUND,net),
                           pf.Constraint(pf.CONSTR_TYPE_FIX,net),
                           pf.Constraint(pf.CONSTR_TYPE_PAR_GEN_P,net),
                           pf.Constraint(pf.CONSTR_TYPE_PAR_GEN_Q,net),
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
            net.load(case)

            # Before updating network
            list(map(lambda c: c.clear_error(),constraints))
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            list(map(lambda c: c.clear_error(),constraints))
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            list(map(lambda c: c.clear_error(),constraints))
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.analyze)
            list(map(lambda c: c.clear_error(),constraints))
            for c in constraints:
                self.assertRaises(pf.ConstraintError,c.eval,x0)
            list(map(lambda c: c.clear_error(),constraints))

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
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            self.assertEqual(net.num_vars,
                             (2*net.num_buses + 
                              2*net.num_gens +
                              net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts()+
                              2*net.num_bats))

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
                self.assertEqual(c.A.shape[1],net.num_vars)
                self.assertEqual(c.J.shape[1],net.num_vars)
                if c.f.size:
                    self.assertTupleEqual(c.get_H_single(0).shape,(net.num_vars,net.num_vars))
                else:
                    self.assertTupleEqual(c.get_H_single(0).shape,(0,0))

    def test_constr_DCPF(self):
                
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertEqual(net.num_vars,0)
            
            # Add vargens
            load_buses = net.get_load_buses()
            net.add_vargens(load_buses,50.,30.,5,0.05)
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
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P)
            self.assertEqual(net.num_vars,
                             (net.num_buses-net.get_num_slack_buses() +
                              net.num_gens +
                              net.num_loads + 
                              net.num_vargens +
                              net.get_num_phase_shifters()+
                              net.num_bats))
            
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
                              net.num_loads + 
                              net.num_vargens +
                              4*net.num_branches - 
                              2*r +
                              2*net.get_num_phase_shifters()+
                              net.num_bats))
            self.assertTupleEqual(b.shape,(net.num_buses,))
            self.assertTupleEqual(f.shape,(0,))
            self.assertTupleEqual(A.shape,(net.num_buses,net.num_vars))
            self.assertEqual(A.nnz,constr.Acounter)
            self.assertTupleEqual(J.shape,(0,net.num_vars)) 
            
            constr.eval(x0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(A.nnz,
                             (net.num_gens +
                              net.num_loads + 
                              net.num_vargens +
                              4*net.num_branches - 
                              2*r +
                              2*net.get_num_phase_shifters()+
                              net.num_bats))
            
            # Extract pieces
            P1 = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)
            P2 = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P)
            P3 = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)
            P4 = net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_PHASE)
            P5 = net.get_var_projection(pf.OBJ_LOAD,pf.LOAD_VAR_P)
            P6 = net.get_var_projection(pf.OBJ_BAT,pf.BAT_VAR_P)

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
            self.assertLess(np.linalg.norm((G*p+R*r-Atheta*theta-Aphi*phi-L*l-B*Pb)-A*x),1e-10)

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
                for gen in bus.gens:
                    mis += gen.P
                for vargen in bus.vargens:
                    mis += vargen.P
                for load in bus.loads:
                    mis -= load.P
                for bat in bus.bats:
                    mis -= bat.P
                for br in bus.branches_from:
                    mis -= br.P_flow_DC
                for br in bus.branches_to:
                    mis += br.P_flow_DC
                self.assertLess(np.abs(mismatches[bus.index]-mis),1e-8)

            # no variables
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            constr.del_matvec()
            constr.analyze()
            f1 = constr.f
            J1 = constr.J
            A1 = constr.A
            b1 = constr.b
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            self.assertEqual(constr.Acounter,0)
            self.assertTupleEqual(b1.shape,(net.num_buses,))
            self.assertTupleEqual(f1.shape,(0,))
            self.assertTupleEqual(A1.shape,(net.num_buses,net.num_vars))
            self.assertEqual(A1.nnz,constr.Acounter)
            self.assertTupleEqual(J1.shape,(0,net.num_vars))
            x1 = net.get_var_values()
            self.assertTrue(type(x1) is np.ndarray)
            self.assertTupleEqual(x1.shape,(net.num_vars,))
            
            mismatches1 = A1*x1-b1
            for bus in net.buses:
                mis = 0
                for gen in bus.gens:
                    mis += gen.P
                for vargen in bus.vargens:
                    mis += vargen.P
                for load in bus.loads:
                    mis -= load.P
                for bat in bus.bats:
                    mis -= bat.P
                for br in bus.branches_from:
                    mis -= br.P_flow_DC
                for br in bus.branches_to:
                    mis += br.P_flow_DC
                self.assertLess(np.abs(mismatches1[bus.index]-mis),1e-8)

    def test_constr_DC_FLOW_LIM(self):
                
        net = self.net
        
        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_vars,0)
            
            # Variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VANG)
            self.assertEqual(net.num_vars,net.num_buses-net.get_num_slack_buses())
            
            x0 = net.get_var_values()
            self.assertTrue(type(x0) is np.ndarray)
            self.assertTupleEqual(x0.shape,(net.num_vars,))
            
            # Constraint
            constr = pf.Constraint(pf.CONSTR_TYPE_DC_FLOW_LIM,net)
            
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
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            self.assertEqual(constr.Gconstr_index,0)
            
            # Analyze
            constr.analyze()
            f = constr.f
            J = constr.J
            A = constr.A
            b = constr.b
            l = constr.l
            u = constr.u
            G = constr.G
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            self.assertEqual(constr.Gconstr_index,net.num_branches)

            self.assertTupleEqual(b.shape,(0,))
            self.assertTupleEqual(f.shape,(0,))
            self.assertTupleEqual(l.shape,(net.num_branches,))
            self.assertTupleEqual(u.shape,(net.num_branches,))
            
            self.assertTupleEqual(A.shape,(0,net.num_vars)) 
            self.assertTupleEqual(J.shape,(0,net.num_vars))
            self.assertTupleEqual(G.shape,(net.num_branches,net.num_vars))
            self.assertEqual(G.nnz,constr.Gcounter)
            
            self.assertTrue(np.all(l <= u))

            num = 0
            for br in net.branches:
                if not br.bus_from.is_slack():
                    num += 1
                if not br.bus_to.is_slack():
                    num += 1
            self.assertEqual(num,constr.Gcounter)

            counter = 0
            index = 0
            for br in net.branches:
                off = 0
                if br.bus_from.is_slack():
                    off = br.b*br.bus_from.v_ang
                else:
                    self.assertEqual(G.row[counter],index)
                    self.assertEqual(G.col[counter],br.bus_from.index_v_ang)
                    self.assertEqual(G.data[counter],-br.b)
                    counter += 1
                if br.bus_to.is_slack():
                    off = -br.b*br.bus_to.v_ang
                else:
                    self.assertEqual(G.row[counter],index)
                    self.assertEqual(G.col[counter],br.bus_to.index_v_ang)
                    self.assertEqual(G.data[counter],br.b)
                    counter += 1
                rating = br.ratingA if br.ratingA > 0 else pf.BRANCH_INF_FLOW
                self.assertEqual(l[index],-rating+off-br.b*br.phase)
                self.assertEqual(u[index],rating+off-br.b*br.phase)
                index += 1
            self.assertEqual(counter,G.nnz)
            self.assertEqual(index,G.shape[0])

            # Flow
            Gx0 = constr.G*x0
            self.assertTupleEqual(Gx0.shape,(net.num_branches,))
            index = 0
            for branch in net.branches:
                bus1 = branch.bus_from
                bus2 = branch.bus_to
                if bus1.is_slack():
                    flow = Gx0[index]-branch.b*(bus1.v_ang-branch.phase)
                elif bus2.is_slack():
                    flow = Gx0[index]-branch.b*(-bus2.v_ang-branch.phase)
                else:
                    flow = Gx0[index]-branch.b*(-branch.phase)
                self.assertLess(np.abs(branch.P_flow_DC-flow),1e-10)
                index += 1

            # Sensitivities
            index = 0
            for branch in net.branches:
                self.assertEqual(branch.sens_P_u_bound,0.)
                self.assertEqual(branch.sens_P_l_bound,0.)
            mu = np.random.randn(net.num_branches)
            pi = np.random.randn(net.num_branches)
            self.assertEqual(constr.G.shape[0],net.num_branches)
            constr.store_sensitivities(None,None,mu,pi)
            for branch in net.branches:
                self.assertEqual(index,branch.index)
                self.assertEqual(branch.sens_P_u_bound,mu[index])
                self.assertEqual(branch.sens_P_l_bound,pi[index])
                index += 1
            self.assertEqual(constr.Jcounter,0)
            self.assertEqual(constr.Acounter,0)
            self.assertEqual(constr.Gcounter,0)
            self.assertEqual(constr.Jconstr_index,0)
            self.assertEqual(constr.Aconstr_index,0)
            self.assertEqual(constr.Gconstr_index,net.num_branches)

    def tearDown(self):
        
        pass
