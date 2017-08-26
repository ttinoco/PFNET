#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import numpy as np

def compare_two_networks(unittest, net, new_net):
    """
    Method for checking if two :class:`Network <pfnet.Network>` objects are held in different 
    memory locations but are otherwise identical

    Parameters
    ----------
    unittest : unittest.TestCase
    net : :class:`Network <pfnet.Network>`
    new_net : :class:`Network <pfnet.Network>`
    """

    norminf = lambda x: norm(x,np.inf) if not np.isscalar(x) else np.abs(x)
    eps = 1e-10

    # Network
    unittest.assertTrue(net is not new_net)
    unittest.assertFalse(net.has_same_data(new_net))
    unittest.assertEqual(net.num_periods,new_net.num_periods)
    unittest.assertEqual(net.base_power,new_net.base_power)

    # Buses
    unittest.assertEqual(net.num_buses, new_net.num_buses)
    for i in range(net.num_buses):
        bus = net.buses[i]
        new_bus = new_net.buses[i]
        unittest.assertTrue(bus is not new_bus)
        unittest.assertEqual(bus.number, new_bus.number)
        unittest.assertEqual(bus.name, new_bus.name)
        unittest.assertLess(norminf(bus.v_base-new_bus.v_base), eps)
        unittest.assertLess(norminf(bus.v_mag-new_bus.v_mag), eps)
        unittest.assertLess(norminf(bus.v_ang-new_bus.v_ang), eps)
        unittest.assertLess(norminf(bus.v_set-new_bus.v_set), eps)
        unittest.assertLess(norminf(bus.v_mag-new_bus.v_mag), eps)
        unittest.assertLess(norminf(bus.v_max_reg-new_bus.v_max_reg), eps)
        unittest.assertLess(norminf(bus.v_min_reg-new_bus.v_min_reg), eps)
        unittest.assertLess(norminf(bus.v_max_norm-new_bus.v_max_norm), eps)
        unittest.assertLess(norminf(bus.v_min_norm-new_bus.v_min_norm), eps)
        unittest.assertLess(norminf(bus.v_max_emer-new_bus.v_max_emer), eps)
        unittest.assertLess(norminf(bus.v_min_emer-new_bus.v_min_emer), eps)
        unittest.assertEqual(bus.is_slack(), new_bus.is_slack())
        unittest.assertEqual(bus.is_regulated_by_gen(),new_bus.is_regulated_by_gen())
        unittest.assertEqual(bus.is_regulated_by_tran(),new_bus.is_regulated_by_tran())
        unittest.assertEqual(bus.is_regulated_by_shunt(),new_bus.is_regulated_by_shunt())
        unittest.assertLess(norminf(bus.price-new_bus.price), eps)
        unittest.assertEqual(set([o.index for o in bus.generators]),
                             set([o.index for o in new_bus.generators]))
        unittest.assertEqual(set([o.index for o in bus.reg_generators]),
                             set([o.index for o in new_bus.reg_generators]))
        unittest.assertEqual(set([o.index for o in bus.loads]),
                             set([o.index for o in new_bus.loads]))
        unittest.assertEqual(set([o.index for o in bus.shunts]),
                             set([o.index for o in new_bus.shunts]))
        unittest.assertEqual(set([o.index for o in bus.reg_shunts]),
                             set([o.index for o in new_bus.reg_shunts]))
        unittest.assertEqual(set([o.index for o in bus.branches_k]),
                             set([o.index for o in new_bus.branches_k]))
        unittest.assertEqual(set([o.index for o in bus.branches_m]),
                             set([o.index for o in new_bus.branches_m]))
        unittest.assertEqual(set([o.index for o in bus.reg_trans]),
                             set([o.index for o in new_bus.reg_trans]))
        unittest.assertEqual(set([o.index for o in bus.var_generators]),
                             set([o.index for o in new_bus.var_generators]))
        unittest.assertEqual(set([o.index for o in bus.batteries]),
                             set([o.index for o in new_bus.batteries]))

    # Branches
    unittest.assertEqual(net.num_branches, new_net.num_branches)
    for i in range(net.num_branches):
        branch = net.branches[i]
        new_branch = new_net.branches[i]
        unittest.assertTrue(branch is not new_branch)
        unittest.assertEqual(branch.num_periods, new_branch.num_periods)
        unittest.assertEqual(branch.bus_k.index, new_branch.bus_k.index)
        unittest.assertEqual(branch.bus_m.index, new_branch.bus_m.index)
        unittest.assertEqual(branch.is_fixed_tran(), new_branch.is_fixed_tran())
        unittest.assertEqual(branch.is_line(), new_branch.is_line())
        unittest.assertEqual(branch.is_phase_shifter(), new_branch.is_phase_shifter())
        unittest.assertEqual(branch.is_tap_changer(), new_branch.is_tap_changer())
        unittest.assertEqual(branch.is_tap_changer_v(), new_branch.is_tap_changer_v())
        unittest.assertEqual(branch.is_tap_changer_Q(), new_branch.is_tap_changer_Q())
        if branch.is_tap_changer_v():
            unittest.assertEqual(branch.reg_bus.index, new_branch.reg_bus.index)
        unittest.assertLess(norminf(branch.g-new_branch.g), eps*(1+norminf(branch.g)))
        unittest.assertLess(norminf(branch.g_k-new_branch.g_k), eps*(1+norminf(branch.g_k)))
        unittest.assertLess(norminf(branch.g_m-new_branch.g_m), eps*(1+norminf(branch.g_m)))
        unittest.assertLess(norminf(branch.b-new_branch.b), eps*(1+norminf(branch.b)))
        unittest.assertLess(norminf(branch.b_k-new_branch.b_k), eps*(1+norminf(branch.b_k)))
        unittest.assertLess(norminf(branch.b_m-new_branch.b_m), eps*(1+norminf(branch.b_m)))
        unittest.assertLess(norminf(branch.ratio-new_branch.ratio), eps)
        unittest.assertLess(norminf(branch.ratio_max-new_branch.ratio_max), eps)
        unittest.assertLess(norminf(branch.ratio_min-new_branch.ratio_min), eps)
        unittest.assertLess(norminf(branch.phase-new_branch.phase), eps)
        unittest.assertLess(norminf(branch.phase_max-new_branch.phase_max), eps)
        unittest.assertLess(norminf(branch.phase_min-new_branch.phase_min), eps)
        unittest.assertLess(norminf(branch.ratingA-new_branch.ratingA), eps)
        unittest.assertLess(norminf(branch.ratingB-new_branch.ratingB), eps)
        unittest.assertLess(norminf(branch.ratingC-new_branch.ratingC), eps)
        unittest.assertEqual(branch.is_on_outage(), new_branch.is_on_outage())
        unittest.assertEqual(branch.has_pos_ratio_v_sens(), new_branch.has_pos_ratio_v_sens())

    # Generators
    unittest.assertEqual(net.num_generators, new_net.num_generators)
    for i in range(net.num_generators):
        gen = net.generators[i]
        new_gen = new_net.generators[i]
        unittest.assertTrue(gen is not new_gen)
        unittest.assertEqual(gen.num_periods, new_gen.num_periods)
        unittest.assertEqual(gen.bus.index, new_gen.bus.index)
        unittest.assertEqual(gen.is_on_outage(), new_gen.is_on_outage())
        unittest.assertEqual(gen.is_slack(), new_gen.is_slack())
        unittest.assertEqual(gen.is_regulator(), new_gen.is_regulator())
        unittest.assertEqual(gen.is_P_adjustable(), new_gen.is_P_adjustable())
        if gen.is_regulator():
            unittest.assertEqual(gen.reg_bus.index, new_gen.reg_bus.index)
        unittest.assertLess(norminf(gen.P-new_gen.P), eps)
        unittest.assertLess(norminf(gen.P_max-new_gen.P_max), eps)
        unittest.assertLess(norminf(gen.P_min-new_gen.P_min), eps)
        unittest.assertLess(norminf(gen.dP_max-new_gen.dP_max), eps)
        unittest.assertLess(norminf(gen.P_prev-new_gen.P_prev), eps)
        unittest.assertLess(norminf(gen.Q-new_gen.Q), eps)
        unittest.assertLess(norminf(gen.Q_max-new_gen.Q_max), eps)
        unittest.assertLess(norminf(gen.Q_min-new_gen.Q_min), eps)
        unittest.assertLess(norminf(gen.cost_coeff_Q0-new_gen.cost_coeff_Q0), eps)
        unittest.assertLess(norminf(gen.cost_coeff_Q1-new_gen.cost_coeff_Q1), eps)
        unittest.assertLess(norminf(gen.cost_coeff_Q2-new_gen.cost_coeff_Q2), eps)

    # Var generators
    unittest.assertEqual(net.num_var_generators, new_net.num_var_generators)
    for i in range(net.num_var_generators):
        vargen = net.var_generators[i]
        new_vargen = new_net.var_generators[i]
        unittest.assertTrue(vargen is not new_vargen)
        unittest.assertEqual(vargen.num_periods, new_vargen.num_periods)
        unittest.assertEqual(vargen.bus.index, new_vargen.bus.index)
        unittest.assertEqual(vargen.name, new_vargen.name)
        unittest.assertLess(norminf(vargen.P-new_vargen.P), eps)
        unittest.assertLess(norminf(vargen.P_ava-new_vargen.P_ava), eps)
        unittest.assertLess(norminf(vargen.P_max-new_vargen.P_max), eps)
        unittest.assertLess(norminf(vargen.P_min-new_vargen.P_min), eps)
        unittest.assertLess(norminf(vargen.P_std-new_vargen.P_std), eps)
        unittest.assertLess(norminf(vargen.Q-new_vargen.Q), eps)
        unittest.assertLess(norminf(vargen.Q_max-new_vargen.Q_max), eps)
        unittest.assertLess(norminf(vargen.Q_min-new_vargen.Q_min), eps)

    # Shunts
    unittest.assertEqual(net.num_shunts, new_net.num_shunts)
    for i in range(net.num_shunts):
        shunt = net.shunts[i]
        new_shunt = new_net.shunts[i]
        unittest.assertTrue(shunt is not new_shunt)
        unittest.assertEqual(shunt.num_periods, new_shunt.num_periods)
        unittest.assertEqual(shunt.bus.index, new_shunt.bus.index)
        unittest.assertEqual(shunt.is_fixed(), new_shunt.is_fixed())
        unittest.assertEqual(shunt.is_switched_v(), new_shunt.is_switched_v())
        if shunt.is_switched_v():
            unittest.assertEqual(shunt.reg_bus.index, new_shunt.reg_bus.index)
        unittest.assertLess(norminf(shunt.g-new_shunt.g),eps*(1+norminf(shunt.g)))
        unittest.assertLess(norminf(shunt.b-new_shunt.b),eps*(1+norminf(shunt.b)))
        unittest.assertLess(norminf(shunt.b_max-new_shunt.b_max), eps*(1+norminf(shunt.b_max)))
        unittest.assertLess(norminf(shunt.b_min-new_shunt.b_min), eps*(1+norminf(shunt.b_min)))

    # Loads
    unittest.assertEqual(net.num_loads, new_net.num_loads)
    for i in range(net.num_loads):
        load = net.loads[i]
        new_load = new_net.loads[i]
        unittest.assertTrue(load is not new_load)
        unittest.assertEqual(load.num_periods, new_load.num_periods)
        unittest.assertEqual(load.bus.index, new_load.bus.index)
        unittest.assertLess(norminf(load.P-new_load.P), eps)
        unittest.assertLess(norminf(load.P_max-new_load.P_max), eps)
        unittest.assertLess(norminf(load.P_min-new_load.P_min), eps)
        unittest.assertLess(norminf(load.Q-new_load.Q), eps)
        unittest.assertLess(norminf(load.target_power_factor-new_load.target_power_factor), eps)
        unittest.assertLess(norminf(load.util_coeff_Q0-new_load.util_coeff_Q0), eps)
        unittest.assertLess(norminf(load.util_coeff_Q1-new_load.util_coeff_Q1), eps)
        unittest.assertLess(norminf(load.util_coeff_Q2-new_load.util_coeff_Q2), eps)

    # Batteries
    unittest.assertEqual(net.num_batteries, new_net.num_batteries)
    for i in range(net.num_batteries):
        bat = net.batteries[i]
        new_bat = new_net.batteries[i]
        unittest.assertTrue(bat is not new_bat)
        unittest.assertEqual(bat.num_periods, new_bat.num_periods)
        unittest.assertEqual(bat.bus.index, new_bat.bus.index)
        unittest.assertLess(norminf(bat.P-new_bat.P), eps)
        unittest.assertLess(norminf(bat.P_max-new_bat.P_max), eps)
        unittest.assertLess(norminf(bat.P_min-new_bat.P_min), eps)
        unittest.assertLess(norminf(bat.eta_c-new_bat.eta_c), eps)
        unittest.assertLess(norminf(bat.eta_d-new_bat.eta_d), eps)
        unittest.assertLess(norminf(bat.E-new_bat.E), eps)
        unittest.assertLess(norminf(bat.E_init-new_bat.E_init), eps)
        unittest.assertLess(norminf(bat.E_final-new_bat.E_final), eps)
        unittest.assertLess(norminf(bat.E_max-new_bat.E_max), eps)

    print 'compared two networks'
