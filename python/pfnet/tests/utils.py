#***************************************************#
# This file is part of PFNET1.                      #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import numpy as np
from numpy.linalg import norm

norminf = lambda x: norm(x,np.inf) if isinstance(x,np.ndarray) else np.abs(x)

def compare_buses(test, bus1, bus2, check_internals=False, check_indices=True, eps=1e-10):
    """
    Method for checking if two |Bus| objects are similar.
    
    Parameters
    ----------
    test : unittest.TestCase
    bus1 : |Bus|
    bus2 : |Bus|
    check_internals : |TrueFalse|
    check_indices : |TrueFalse|
    eps : float
    """

    test.assertTrue(bus1 is not bus2)
    test.assertEqual(bus1.number, bus2.number)
    test.assertEqual(bus1.num_periods, bus2.num_periods)
    test.assertEqual(bus1.name.upper().strip(), bus2.name.upper().strip())
    test.assertLess(norminf(bus1.v_base-bus2.v_base), eps)
    test.assertLess(norminf(bus1.v_mag-bus2.v_mag), eps)
    test.assertLess(norminf(bus1.v_ang-bus2.v_ang), eps)
    test.assertLess(norminf(bus1.v_set-bus2.v_set), eps)
    test.assertLess(norminf(bus1.v_max_reg-bus2.v_max_reg), eps)
    test.assertLess(norminf(bus1.v_min_reg-bus2.v_min_reg), eps)
    test.assertLess(norminf(bus1.v_max_norm-bus2.v_max_norm), eps)
    test.assertLess(norminf(bus1.v_min_norm-bus2.v_min_norm), eps)
    test.assertLess(norminf(bus1.v_max_emer-bus2.v_max_emer), eps)
    test.assertLess(norminf(bus1.v_min_emer-bus2.v_min_emer), eps)
    test.assertEqual(bus1.is_slack(), bus2.is_slack())
    test.assertEqual(bus1.is_regulated_by_gen(),bus2.is_regulated_by_gen())
    test.assertEqual(bus1.is_regulated_by_tran(),bus2.is_regulated_by_tran())
    test.assertEqual(bus1.is_regulated_by_shunt(),bus2.is_regulated_by_shunt())
    test.assertLess(norminf(bus1.price-bus2.price), eps)
    test.assertEqual(len(bus1.generators),len(bus2.generators))
    test.assertEqual(len(bus1.reg_generators),len(bus2.reg_generators))
    test.assertEqual(len(bus1.loads),len(bus2.loads))
    test.assertEqual(len(bus1.shunts),len(bus2.shunts))
    test.assertEqual(len(bus1.branches_k),len(bus2.branches_k))
    test.assertEqual(len(bus1.branches_m),len(bus2.branches_m))
    test.assertEqual(len(bus1.batteries),len(bus2.batteries))
    test.assertEqual(len(bus1.var_generators),len(bus2.var_generators))
    test.assertEqual(len(bus1.reg_trans),len(bus2.reg_trans))
    test.assertEqual(len(bus1.reg_shunts),len(bus2.reg_shunts))
    if check_indices:
        test.assertEqual(set([o.index for o in bus1.generators]),
                         set([o.index for o in bus2.generators]))
        test.assertEqual(set([o.index for o in bus1.reg_generators]),
                         set([o.index for o in bus2.reg_generators]))
        test.assertEqual(set([o.index for o in bus1.loads]),
                         set([o.index for o in bus2.loads]))
        test.assertEqual(set([o.index for o in bus1.shunts]),
                         set([o.index for o in bus2.shunts]))
        test.assertEqual(set([o.index for o in bus1.reg_shunts]),
                         set([o.index for o in bus2.reg_shunts]))
        test.assertEqual(set([o.index for o in bus1.branches_k]),
                         set([o.index for o in bus2.branches_k]))
        test.assertEqual(set([o.index for o in bus1.branches_m]),
                         set([o.index for o in bus2.branches_m]))
        test.assertEqual(set([o.index for o in bus1.reg_trans]),
                         set([o.index for o in bus2.reg_trans]))
        test.assertEqual(set([o.index for o in bus1.var_generators]),
                         set([o.index for o in bus2.var_generators]))
        test.assertEqual(set([o.index for o in bus1.batteries]),
                         set([o.index for o in bus2.batteries]))
    test.assertLess(norminf(bus1.sens_P_balance-bus2.sens_P_balance), eps)
    test.assertLess(norminf(bus1.sens_Q_balance-bus2.sens_Q_balance), eps)
    test.assertLess(norminf(bus1.sens_v_mag_u_bound-bus2.sens_v_mag_u_bound), eps)
    test.assertLess(norminf(bus1.sens_v_mag_l_bound-bus2.sens_v_mag_l_bound), eps)
    test.assertLess(norminf(bus1.sens_v_ang_u_bound-bus2.sens_v_ang_u_bound), eps)
    test.assertLess(norminf(bus1.sens_v_ang_l_bound-bus2.sens_v_ang_l_bound), eps)
    test.assertLess(norminf(bus1.sens_v_reg_by_gen-bus2.sens_v_reg_by_gen), eps)
    test.assertLess(norminf(bus1.sens_v_reg_by_tran-bus2.sens_v_reg_by_tran), eps)
    test.assertLess(norminf(bus1.sens_v_reg_by_shunt-bus2.sens_v_reg_by_shunt), eps)
    if check_internals:
        test.assertLess(norminf(bus1.index_v_mag-bus2.index_v_mag),eps)
        test.assertLess(norminf(bus1.index_v_ang-bus2.index_v_ang),eps)
        test.assertEqual(bus1.flags_vars,bus2.flags_vars)
        test.assertEqual(bus1.flags_fixed,bus2.flags_fixed)
        test.assertEqual(bus1.flags_bounded,bus2.flags_bounded)
        test.assertEqual(bus1.flags_sparse,bus2.flags_sparse)

def compare_generators(test, gen1, gen2, check_internals=False, eps=1e-10):
    """
    Method for checking if two |Generator| objects are similar.
    
    Parameters
    ----------
    test : unittest.TestCase
    gen1 : |Generator|
    gen2 : |Generator|
    check_configuration : |TrueFalse|
    eps : float
    """
    
    test.assertTrue(gen1 is not gen2)
    test.assertEqual(gen1.name, gen2.name)
    test.assertEqual(gen1.num_periods, gen2.num_periods)
    test.assertEqual(gen1.bus.index, gen2.bus.index)
    test.assertEqual(gen1.is_on_outage(), gen2.is_on_outage())
    test.assertEqual(gen1.is_slack(), gen2.is_slack())
    test.assertEqual(gen1.is_regulator(), gen2.is_regulator())
    test.assertEqual(gen1.is_P_adjustable(), gen2.is_P_adjustable())
    if gen1.is_regulator():
        test.assertEqual(gen1.reg_bus.index, gen2.reg_bus.index)
    test.assertLess(norminf(gen1.P-gen2.P), eps)
    test.assertLess(norminf(gen1.P_max-gen2.P_max), eps)
    test.assertLess(norminf(gen1.P_min-gen2.P_min), eps)
    test.assertLess(norminf(gen1.dP_max-gen2.dP_max), eps)
    test.assertLess(norminf(gen1.P_prev-gen2.P_prev), eps)
    test.assertLess(norminf(gen1.Q-gen2.Q), eps)
    test.assertLess(norminf(gen1.Q_max-gen2.Q_max), eps)
    test.assertLess(norminf(gen1.Q_min-gen2.Q_min), eps)
    test.assertLess(norminf(gen1.cost_coeff_Q0-gen2.cost_coeff_Q0), eps)
    test.assertLess(norminf(gen1.cost_coeff_Q1-gen2.cost_coeff_Q1), eps)
    test.assertLess(norminf(gen1.cost_coeff_Q2-gen2.cost_coeff_Q2), eps)
    if check_internals:
        test.assertLess(norminf(gen1.index_P-gen2.index_P),eps)
        test.assertLess(norminf(gen1.index_Q-gen2.index_Q),eps)
        test.assertEqual(gen1.flags_vars,gen2.flags_vars)
        test.assertEqual(gen1.flags_fixed,gen2.flags_fixed)
        test.assertEqual(gen1.flags_bounded,gen2.flags_bounded)
        test.assertEqual(gen1.flags_sparse,gen2.flags_sparse)
        test.assertLess(norminf(gen1.sens_P_u_bound-gen2.sens_P_u_bound), eps)
        test.assertLess(norminf(gen1.sens_P_l_bound-gen2.sens_P_l_bound), eps)
        test.assertLess(norminf(gen1.sens_Q_u_bound-gen2.sens_Q_u_bound), eps)
        test.assertLess(norminf(gen1.sens_Q_l_bound-gen2.sens_Q_l_bound), eps)

def compare_loads(test, load1, load2, check_internals=False, eps=1e-10):
    """
    Method for checking if two |Load| objects are similar.
    
    Parameters
    ----------
    test : unittest.TestCase
    load1 : |Load|
    load2 : |Load|
    check_configuration : |TrueFalse|
    eps : float
    """

    test.assertTrue(load1 is not load2)
    test.assertEqual(load1.name, load2.name)
    test.assertEqual(load1.num_periods, load2.num_periods)
    test.assertEqual(load1.bus.index, load2.bus.index)
    test.assertLess(norminf(load1.P-load2.P), eps)
    test.assertLess(norminf(load1.P_max-load2.P_max), eps)
    test.assertLess(norminf(load1.P_min-load2.P_min), eps)
    test.assertLess(norminf(load1.Q-load2.Q), eps)
    test.assertLess(norminf(load1.target_power_factor-load2.target_power_factor), eps)
    test.assertLess(norminf(load1.util_coeff_Q0-load2.util_coeff_Q0), eps)
    test.assertLess(norminf(load1.util_coeff_Q1-load2.util_coeff_Q1), eps)
    test.assertLess(norminf(load1.util_coeff_Q2-load2.util_coeff_Q2), eps)
    if check_internals:
        test.assertLess(norminf(load1.index_P-load2.index_P),eps)
        test.assertLess(norminf(load1.index_Q-load2.index_Q),eps)
        test.assertEqual(load1.flags_vars,load2.flags_vars)
        test.assertEqual(load1.flags_fixed,load2.flags_fixed)
        test.assertEqual(load1.flags_bounded,load2.flags_bounded)
        test.assertEqual(load1.flags_sparse,load2.flags_sparse)
        test.assertLess(norminf(load1.sens_P_u_bound-load2.sens_P_u_bound), eps)
        test.assertLess(norminf(load1.sens_P_l_bound-load2.sens_P_l_bound), eps)

def compare_shunts(test, shunt1, shunt2, check_internals=False, eps=1e-10):
    """
    Method for checking if two |Shunt| objects are similar.
    
    Parameters
    ----------
    test : unittest.TestCase
    shunt1 : |Shunt|
    shunt2 : |Shunt|
    check_configuration : |TrueFalse|
    eps : float
    """

    test.assertTrue(shunt1 is not shunt2)
    test.assertEqual(shunt1.name, shunt2.name)
    test.assertEqual(shunt1.num_periods, shunt2.num_periods)
    test.assertEqual(shunt1.bus.index, shunt2.bus.index)
    test.assertEqual(shunt1.is_fixed(), shunt2.is_fixed())
    test.assertEqual(shunt1.is_switched_v(), shunt2.is_switched_v())
    if shunt1.is_switched_v():
        test.assertEqual(shunt1.reg_bus.index, shunt2.reg_bus.index)
    test.assertLess(norminf(shunt1.g-shunt2.g),eps*(1+norminf(shunt1.g)))
    test.assertLess(norminf(shunt1.b-shunt2.b),eps*(1+norminf(shunt1.b)))
    test.assertLess(norminf(shunt1.b_max-shunt2.b_max), eps*(1+norminf(shunt1.b_max)))
    test.assertLess(norminf(shunt1.b_min-shunt2.b_min), eps*(1+norminf(shunt1.b_min)))
    if check_internals:
        test.assertLess(norminf(shunt1.index_b-shunt2.index_b),eps)
        test.assertEqual(shunt1.flags_vars,shunt2.flags_vars)
        test.assertEqual(shunt1.flags_fixed,shunt2.flags_fixed)
        test.assertEqual(shunt1.flags_bounded,shunt2.flags_bounded)
        test.assertEqual(shunt1.flags_sparse,shunt2.flags_sparse)
        test.assertLess(norminf(shunt1.sens_b_u_bound-shunt2.sens_b_u_bound), eps)
        test.assertLess(norminf(shunt1.sens_b_l_bound-shunt2.sens_b_l_bound), eps)

def compare_branches(test, branch1, branch2, check_internals=False, eps=1e-10):
    """
    Method for checking if two |Branch| objects are similar.
    
    Parameters
    ----------
    test : unittest.TestCase
    branch1 : |Branch|
    branch2 : |Branch|
    check_configuration : |TrueFalse|
    eps : float
    """

    test.assertTrue(branch1 is not branch2)
    test.assertEqual(branch1.name, branch2.name)
    test.assertEqual(branch1.num_periods, branch2.num_periods)
    test.assertEqual(branch1.bus_k.index, branch2.bus_k.index)
    test.assertEqual(branch1.bus_m.index, branch2.bus_m.index)
    test.assertEqual(branch1.is_fixed_tran(), branch2.is_fixed_tran())
    test.assertEqual(branch1.is_line(), branch2.is_line())
    test.assertEqual(branch1.is_phase_shifter(), branch2.is_phase_shifter())
    test.assertEqual(branch1.is_tap_changer(), branch2.is_tap_changer())
    test.assertEqual(branch1.is_tap_changer_v(), branch2.is_tap_changer_v())
    test.assertEqual(branch1.is_tap_changer_Q(), branch2.is_tap_changer_Q())
    if branch1.is_tap_changer_v():
        test.assertEqual(branch1.reg_bus.index, branch2.reg_bus.index)
    test.assertLess(norminf(branch1.g-branch2.g), eps*(1+norminf(branch1.g)))
    test.assertLess(norminf(branch1.g_k-branch2.g_k), eps*(1+norminf(branch1.g_k)))
    test.assertLess(norminf(branch1.g_m-branch2.g_m), eps*(1+norminf(branch1.g_m)))
    test.assertLess(norminf(branch1.b-branch2.b), eps*(1+norminf(branch1.b)))
    test.assertLess(norminf(branch1.b_k-branch2.b_k), eps*(1+norminf(branch1.b_k)))
    test.assertLess(norminf(branch1.b_m-branch2.b_m), eps*(1+norminf(branch1.b_m)))
    test.assertLess(norminf(branch1.ratio-branch2.ratio), eps)
    test.assertLess(norminf(branch1.ratio_max-branch2.ratio_max), eps)
    test.assertLess(norminf(branch1.ratio_min-branch2.ratio_min), eps)
    test.assertLess(norminf(branch1.phase-branch2.phase), eps)
    test.assertLess(norminf(branch1.phase_max-branch2.phase_max), eps)
    test.assertLess(norminf(branch1.phase_min-branch2.phase_min), eps)
    test.assertLess(norminf(branch1.ratingA-branch2.ratingA), eps)
    test.assertLess(norminf(branch1.ratingB-branch2.ratingB), eps)
    test.assertLess(norminf(branch1.ratingC-branch2.ratingC), eps)
    test.assertEqual(branch1.is_on_outage(), branch2.is_on_outage())
    test.assertEqual(branch1.has_pos_ratio_v_sens(), branch2.has_pos_ratio_v_sens())
    if check_internals:
        test.assertLess(norminf(branch1.index_ratio-branch2.index_ratio),eps)
        test.assertLess(norminf(branch1.index_phase-branch2.index_phase),eps)
        test.assertEqual(branch1.flags_vars,branch2.flags_vars)
        test.assertEqual(branch1.flags_fixed,branch2.flags_fixed)
        test.assertEqual(branch1.flags_bounded,branch2.flags_bounded)
        test.assertEqual(branch1.flags_sparse,branch2.flags_sparse)
        test.assertLess(norminf(branch1.sens_P_u_bound-branch2.sens_P_u_bound), eps)
        test.assertLess(norminf(branch1.sens_P_l_bound-branch2.sens_P_l_bound), eps)
        test.assertLess(norminf(branch1.sens_i_mag_u_bound-branch2.sens_i_mag_u_bound), eps)
        test.assertLess(norminf(branch1.sens_ratio_u_bound-branch2.sens_ratio_u_bound), eps)
        test.assertLess(norminf(branch1.sens_ratio_l_bound-branch2.sens_ratio_l_bound), eps)
        test.assertLess(norminf(branch1.sens_phase_u_bound-branch2.sens_phase_u_bound), eps)
        test.assertLess(norminf(branch1.sens_phase_l_bound-branch2.sens_phase_l_bound), eps)
        
def compare_networks(test, net1, net2, check_internals=False):
    """
    Method for checking if two |Network| objects are held in different
    memory locations but are otherwise identical.

    Parameters
    ----------
    test : unittest.TestCase
    net1 : |Network|
    net2 : |Network|
    check_configuration : |TrueFalse|
    """

    eps = 1e-10

    # Network
    test.assertTrue(net1 is not net2)
    test.assertFalse(net1.has_same_ptr(net2))
    test.assertEqual(net1.num_periods,net2.num_periods)
    test.assertEqual(net1.base_power,net2.base_power)
    if check_internals:
        test.assertEqual(net1.has_error(),net2.has_error())
        test.assertEqual(net1.error_string,net2.error_string)
        test.assertEqual(net1.num_vars,net2.num_vars)
        test.assertEqual(net1.num_fixed,net2.num_fixed)
        test.assertEqual(net1.num_bounded,net2.num_bounded)
        test.assertEqual(net1.num_sparse,net2.num_sparse)
    
    # Buses
    test.assertEqual(net1.num_buses, net2.num_buses)
    for i in range(net1.num_buses):
        bus1 = net1.buses[i]
        bus2 = net2.buses[i]
        compare_buses(test, bus1, bus2, check_internals=check_internals, eps=eps) 

    # Branches
    test.assertEqual(net1.num_branches, net2.num_branches)
    for i in range(net1.num_branches):
        branch1 = net1.branches[i]
        branch2 = net2.branches[i]
        compare_branches(test, branch1, branch2, check_internals=check_internals, eps=eps)

    # Generators
    test.assertEqual(net1.num_generators, net2.num_generators)
    for i in range(net1.num_generators):
        gen1 = net1.generators[i]
        gen2 = net2.generators[i]
        compare_generators(test, gen1, gen2, check_internals=check_internals, eps=eps)
        
    # Var generators
    test.assertEqual(net1.num_var_generators, net2.num_var_generators)
    for i in range(net1.num_var_generators):
        vargen1 = net1.var_generators[i]
        vargen2 = net2.var_generators[i]
        test.assertTrue(vargen1 is not vargen2)
        test.assertEqual(vargen1.num_periods, vargen2.num_periods)
        test.assertEqual(vargen1.bus.index, vargen2.bus.index)
        test.assertEqual(vargen1.name, vargen2.name)
        test.assertLess(norminf(vargen1.P-vargen2.P), eps)
        test.assertLess(norminf(vargen1.P_ava-vargen2.P_ava), eps)
        test.assertLess(norminf(vargen1.P_max-vargen2.P_max), eps)
        test.assertLess(norminf(vargen1.P_min-vargen2.P_min), eps)
        test.assertLess(norminf(vargen1.P_std-vargen2.P_std), eps)
        test.assertLess(norminf(vargen1.Q-vargen2.Q), eps)
        test.assertLess(norminf(vargen1.Q_max-vargen2.Q_max), eps)
        test.assertLess(norminf(vargen1.Q_min-vargen2.Q_min), eps)
        if check_internals:
            test.assertLess(norminf(vargen1.index_P-vargen2.index_P),eps)
            test.assertLess(norminf(vargen1.index_Q-vargen2.index_Q),eps)
            test.assertEqual(vargen1.flags_vars,vargen2.flags_vars)
            test.assertEqual(vargen1.flags_fixed,vargen2.flags_fixed)
            test.assertEqual(vargen1.flags_bounded,vargen2.flags_bounded)
            test.assertEqual(vargen1.flags_sparse,vargen2.flags_sparse)

    # Shunts
    test.assertEqual(net1.num_shunts, net2.num_shunts)
    for i in range(net1.num_shunts):
        shunt1 = net1.shunts[i]
        shunt2 = net2.shunts[i]
        compare_shunts(test, shunt1, shunt2, check_internals=check_internals, eps=eps) 

    # Loads
    test.assertEqual(net1.num_loads, net2.num_loads)
    for i in range(net1.num_loads):
        load1 = net1.loads[i]
        load2 = net2.loads[i]
        compare_loads(test, load1, load2, check_internals=check_internals, eps=eps)
        
    # Batteries
    test.assertEqual(net1.num_batteries, net2.num_batteries)
    for i in range(net1.num_batteries):
        bat1 = net1.batteries[i]
        bat2 = net2.batteries[i]
        test.assertTrue(bat1 is not bat2)
        test.assertEqual(bat1.name, bat2.name)
        test.assertEqual(bat1.num_periods, bat2.num_periods)
        test.assertEqual(bat1.bus.index, bat2.bus.index)
        test.assertLess(norminf(bat1.P-bat2.P), eps)
        test.assertLess(norminf(bat1.P_max-bat2.P_max), eps)
        test.assertLess(norminf(bat1.P_min-bat2.P_min), eps)
        test.assertLess(norminf(bat1.eta_c-bat2.eta_c), eps)
        test.assertLess(norminf(bat1.eta_d-bat2.eta_d), eps)
        test.assertLess(norminf(bat1.E-bat2.E), eps)
        test.assertLess(norminf(bat1.E_init-bat2.E_init), eps)
        test.assertLess(norminf(bat1.E_final-bat2.E_final), eps)
        test.assertLess(norminf(bat1.E_max-bat2.E_max), eps)
        if check_internals:
            test.assertLess(norminf(bat1.index_Pc-bat2.index_Pc),eps)
            test.assertLess(norminf(bat1.index_Pd-bat2.index_Pd),eps)
            test.assertLess(norminf(bat1.index_E-bat2.index_E),eps)
            test.assertEqual(bat1.flags_vars,bat2.flags_vars)
            test.assertEqual(bat1.flags_fixed,bat2.flags_fixed)
            test.assertEqual(bat1.flags_bounded,bat2.flags_bounded)
            test.assertEqual(bat1.flags_sparse,bat2.flags_sparse)

    # Hashes                 
    for bus in net1.buses:
        test.assertEqual(bus.number,net2.get_bus_from_number(bus.number).number)
        test.assertEqual(bus.name,net2.get_bus_from_number(bus.number).name)

def check_network(test, net):
    """
    Checks validity of network.

    Parameters
    ----------
    test :  unittest.TestCase
    net : |Network|
    """

    test.assertTrue(isinstance(net, pf.Network))

    test.assertEqual(len(net.buses), net.num_buses)
    test.assertEqual(len(net.branches), net.num_branches)
    test.assertEqual(len(net.generators), net.num_generators)
    test.assertEqual(len(net.loads), net.num_loads)
    test.assertEqual(len(net.shunts), net.num_shunts)
    test.assertEqual(len(net.batteries), net.num_batteries)
    test.assertEqual(len(net.var_generators), net.num_var_generators)
    
    for bus in net.buses:
        for gen in bus.generators:
            test.assertTrue(gen.bus.is_equal(bus))
            test.assertTrue(gen.is_equal(net.get_generator(gen.index)))
        for reg_gen in bus.reg_generators:
            test.assertTrue(reg_gen.reg_bus.is_equal(bus))
            test.assertTrue(reg_gen.is_equal(net.get_generator(reg_gen.index)))
        for brk in bus.branches_k:
            test.assertTrue(brk.bus_k.is_equal(bus))
            test.assertTrue(brk.is_equal(net.get_branch(brk.index)))
        for brm in bus.branches_m:
            test.assertTrue(brm.bus_m.is_equal(bus))
            test.assertTrue(brm.is_equal(net.get_branch(brm.index)))
        for br in bus.reg_trans:
            test.assertTrue(br.reg_bus.is_equal(bus))
            test.assertTrue(br.is_equal(net.get_branch(br.index)))
        for br in bus.branches:
            test.assertTrue(br.bus_k.is_equal(bus) or br.bus_m.is_equal(bus))
            test.assertTrue(br.is_equal(net.get_branch(br.index)))
        for load in bus.loads:
            test.assertTrue(load.bus.is_equal(bus))
            test.assertTrue(load.is_equal(net.get_load(load.index)))
        for shunt in bus.shunts:
            test.assertTrue(shunt.bus.is_equal(bus))
            test.assertTrue(shunt.is_equal(net.get_shunt(shunt.index)))
        for reg_shunt in bus.reg_shunts:
            test.assertTrue(reg_shunt.reg_bus.is_equal(bus))
            test.assertTrue(reg_shunt.is_equal(net.get_shunt(reg_shunt.index)))
        for bat in bus.batteries:
            test.assertTrue(bat.bus.is_equal(bus))
            test.assertTrue(bat.is_equal(net.get_battery(bat.index)))
        for vargen in bus.var_generators:
            test.assertTrue(vargen.bus.is_equal(bus))
            test.assertTrue(vargen.is_equal(net.get_var_generator(vargen.index)))

    for gen in net.generators:
        test.assertTrue(gen.index in [x.index for x in gen.bus.generators])
        test.assertTrue(gen.bus.is_equal(net.get_bus(gen.bus.index)))
        if gen.is_regulator():
            test.assertTrue(gen.index in [x.index for x in gen.reg_bus.reg_generators])
    for br in net.branches:
        test.assertTrue(br.bus_k.is_equal(net.get_bus(br.bus_k.index)))
        test.assertTrue(br.bus_m.is_equal(net.get_bus(br.bus_m.index)))
        test.assertTrue(br.index in [x.index for x in br.bus_k.branches_k])
        test.assertTrue(br.index in [x.index for x in br.bus_m.branches_m])
        if br.is_tap_changer_v():
            test.assertTrue(br.index in [x.index for x in br.reg_bus.reg_trans])
    for load in net.loads:
        test.assertTrue(load.bus.is_equal(net.get_bus(load.bus.index)))
        test.assertTrue(load.index in [x.index for x in load.bus.loads])
    for shunt in net.shunts:
        test.assertTrue(shunt.bus.is_equal(net.get_bus(shunt.bus.index)))
        test.assertTrue(shunt.index in [x.index for x in shunt.bus.shunts])
        if shunt.is_switched_v():
            test.assertTrue(shunt.index in [x.index for x in shunt.reg_bus.reg_shunts])
    for bat in net.batteries:
        test.assertTrue(bat.bus.is_equal(net.get_bus(bat.bus.index)))
        test.assertTrue(bat.index in [x.index for x in bat.bus.batteries])
    for gen in net.var_generators:
        test.assertTrue(gen.bus.is_equal(net.get_bus(gen.bus.index)))
        test.assertTrue(gen.index in [x.index for x in gen.bus.var_generators])

    for bus in net.buses:
        test.assertEqual(bus.number, net.get_bus_from_number(bus.number).number)
        test.assertEqual(bus.name, net.get_bus_from_name(bus.name).name)
        test.assertTrue(bus.is_equal(net.get_bus_from_number(bus.number)))
            
