#***************************************************#
#
#
#
#
#***************************************************#

"""
This script is for holding methods commonly used by unittests
(and potentially others)
"""

import pfnet as pf
import numpy as np

#####

def main():
    data_file = 'data/ieee14_v33.raw'
    data_file2 = 'data/ieee25.raw'
    # Using C raw parser
    raw_parser=pf.ParserRAW()
    net1 = raw_parser.parse(data_file)
    net2 = raw_parser.parse(data_file)
    net3 = raw_parser.parse(data_file2)

    print "Testing identical networks"
    try:
        compare_two_networks(net1, net2)
        print "Success"
    except AssertionError:
        # self.fail("myFunc() raised ExceptionType unexpectedly!")
        print "Failed"

    print "Testing different networks"
    try:
        compare_two_networks(net1, net3)
        print "Failed"
    except AssertionError:
        print "Success"

#####

def compare_two_networks(net1, net2):
    """
    Method for checking if two :class:`Network` objects are held in different memory locations but otherwise identical

    Parameters
	----------
	net1, net2 : Network
		The two Network objects to be compared

    Raises an AssertionError exception if all network elements do not match or objects are not distinct
    """

    norminf = lambda x: norm(x,np.inf) if not np.isscalar(x) else np.abs(x)
    eps = 1e-10

    # Network
    assert (net1 is not net2)
    assert (net1.num_periods == net1.num_periods)
    assert (net1.base_power == net2.base_power)
    assert (net1.num_vars == net2.num_vars)

    # Buses
    assert(net1.num_buses == net2.num_buses)
    for i in range(net1.num_buses):
        bus = net1.buses[i]
        new_bus = net2.buses[i]
        assert(bus is not new_bus)
        assert(bus.number == new_bus.number)
        assert(bus.name == new_bus.name)
        assert(norminf(bus.v_mag - new_bus.v_mag) < eps)
        assert(norminf(bus.v_ang - new_bus.v_ang) < eps)
        assert(norminf(bus.v_set - new_bus.v_set) < eps)
        assert(norminf(bus.v_max_reg - new_bus.v_max_reg) < eps)
        assert(norminf(bus.v_min_reg - new_bus.v_min_reg) < eps)
        assert(norminf(bus.v_max_norm - new_bus.v_max_norm) < eps)
        assert(norminf(bus.v_min_norm - new_bus.v_min_norm) < eps)
        assert(norminf(bus.v_max_emer - new_bus.v_max_emer) < eps)
        assert(norminf(bus.v_min_emer - new_bus.v_min_emer) < eps)
        assert(bus.is_slack() == new_bus.is_slack())
        assert(bus.is_regulated_by_gen() == new_bus.is_regulated_by_gen())
        assert(bus.is_regulated_by_tran() == new_bus.is_regulated_by_tran())
        assert(bus.is_regulated_by_shunt() == new_bus.is_regulated_by_shunt())
        assert(norminf(bus.price-new_bus.price) < eps)
        assert(set([o.index for o in bus.generators]) == set([o.index for o in new_bus.generators]))
        assert(set([o.index for o in bus.reg_generators]) == set([o.index for o in new_bus.reg_generators]))
        assert(set([o.index for o in bus.loads]) == set([o.index for o in new_bus.loads]))
        assert(set([o.index for o in bus.shunts]) == set([o.index for o in new_bus.shunts]))
        assert(set([o.index for o in bus.reg_shunts]) == set([o.index for o in new_bus.reg_shunts]))
        assert(set([o.index for o in bus.branches_k]) == set([o.index for o in new_bus.branches_k]))
        assert(set([o.index for o in bus.branches_m]) == set([o.index for o in new_bus.branches_m]))
        assert(set([o.index for o in bus.reg_trans]) == set([o.index for o in new_bus.reg_trans]))
        assert(set([o.index for o in bus.var_generators]) == set([o.index for o in new_bus.var_generators]))
        assert(set([o.index for o in bus.batteries]) == set([o.index for o in new_bus.batteries]))

    # Branches
    assert(net1.num_branches == net2.num_branches)
    for i in range(net1.num_branches):
        branch = net1.branches[i]
        new_branch = net2.branches[i]
        assert(branch is not new_branch)
        assert(branch.num_periods == new_branch.num_periods)
        assert(branch.bus_k.index == new_branch.bus_k.index)
        assert(branch.bus_m.index == new_branch.bus_m.index)
        assert(branch.is_fixed_tran() == new_branch.is_fixed_tran())
        assert(branch.is_line() == new_branch.is_line())
        assert(branch.is_phase_shifter() == new_branch.is_phase_shifter())
        assert(branch.is_tap_changer() == new_branch.is_tap_changer())
        assert(branch.is_tap_changer_v() == new_branch.is_tap_changer_v())
        assert(branch.is_tap_changer_Q() == new_branch.is_tap_changer_Q())
        if branch.is_tap_changer_v():
            assert(branch.reg_bus.index == new_branch.reg_bus.index)
        assert(norminf(branch.g-new_branch.g) < eps*(1+norminf(branch.g)))
        assert(norminf(branch.g_k-new_branch.g_k) < eps*(1+norminf(branch.g_k)))
        assert(norminf(branch.g_m-new_branch.g_m) < eps*(1+norminf(branch.g_m)))
        assert(norminf(branch.b-new_branch.b) < eps*(1+norminf(branch.b)))
        assert(norminf(branch.b_k-new_branch.b_k) < eps*(1+norminf(branch.b_k)))
        assert(norminf(branch.b_m-new_branch.b_m) < eps*(1+norminf(branch.b_m)))
        assert(norminf(branch.ratio-new_branch.ratio) < eps)
        assert(norminf(branch.ratio_max-new_branch.ratio_max) < eps)
        assert(norminf(branch.ratio_min-new_branch.ratio_min) < eps)
        assert(norminf(branch.phase-new_branch.phase) < eps)
        assert(norminf(branch.phase_max-new_branch.phase_max) < eps)
        assert(norminf(branch.phase_min-new_branch.phase_min) < eps)
        assert(norminf(branch.ratingA-new_branch.ratingA) < eps)
        assert(norminf(branch.ratingB-new_branch.ratingB) < eps)
        assert(norminf(branch.ratingC-new_branch.ratingC) < eps)
        assert(branch.is_on_outage() == new_branch.is_on_outage())
        assert(branch.has_pos_ratio_v_sens() == new_branch.has_pos_ratio_v_sens())

    # Generators
    assert(net1.num_generators == net2.num_generators)
    for i in range(net1.num_generators):
        gen = net1.generators[i]
        new_gen = net2.generators[i]
        assert(gen is not new_gen)
        assert(gen.num_periods == new_gen.num_periods)
        assert(gen.bus.index == new_gen.bus.index)
        assert(gen.is_on_outage() == new_gen.is_on_outage())
        assert(gen.is_slack() == new_gen.is_slack())
        assert(gen.is_regulator() == new_gen.is_regulator())
        assert(gen.is_P_adjustable() == new_gen.is_P_adjustable())
        if gen.is_regulator():
            assert(gen.reg_bus.index == new_gen.reg_bus.index)
        assert(norminf(gen.P-new_gen.P) < eps)
        assert(norminf(gen.P_max-new_gen.P_max) < eps)
        assert(norminf(gen.P_min-new_gen.P_min) < eps)
        assert(norminf(gen.dP_max-new_gen.dP_max) < eps)
        assert(norminf(gen.P_prev-new_gen.P_prev) < eps)
        assert(norminf(gen.Q-new_gen.Q) < eps)
        assert(norminf(gen.Q_max-new_gen.Q_max) < eps)
        assert(norminf(gen.Q_min-new_gen.Q_min) < eps)
        assert(norminf(gen.cost_coeff_Q0-new_gen.cost_coeff_Q0) < eps)
        assert(norminf(gen.cost_coeff_Q1-new_gen.cost_coeff_Q1) < eps)
        assert(norminf(gen.cost_coeff_Q2-new_gen.cost_coeff_Q2) < eps)

    # Var generators
    assert(net1.num_var_generators == net2.num_var_generators)
    for i in range(net1.num_var_generators):
        vargen = net1.var_generators[i]
        new_vargen = net2.var_generators[i]
        assert(vargen is not new_vargen)
        assert(vargen.num_periods == new_vargen.num_periods)
        assert(vargen.bus.index == new_vargen.bus.index)
        assert(vargen.name == new_vargen.name)
        assert(norminf(vargen.P-new_vargen.P) < eps)
        assert(norminf(vargen.P_ava-new_vargen.P_ava) < eps)
        assert(norminf(vargen.P_max-new_vargen.P_max) < eps)
        assert(norminf(vargen.P_min-new_vargen.P_min) < eps)
        assert(norminf(vargen.P_std-new_vargen.P_std) < eps)
        assert(norminf(vargen.Q-new_vargen.Q) < eps)
        assert(norminf(vargen.Q_max-new_vargen.Q_max) < eps)
        assert(norminf(vargen.Q_min-new_vargen.Q_min) < eps)

    # Shunts
    assert(net1.num_shunts == net2.num_shunts)
    for i in range(net1.num_shunts):
        shunt = net1.shunts[i]
        new_shunt = net2.shunts[i]
        assert(shunt is not new_shunt)
        assert(shunt.num_periods == new_shunt.num_periods)
        assert(shunt.bus.index == new_shunt.bus.index)
        assert(shunt.is_fixed() == new_shunt.is_fixed())
        assert(shunt.is_switched_v() == new_shunt.is_switched_v())
        if shunt.is_switched_v():
            assert(shunt.reg_bus.index == new_shunt.reg_bus.index)
        assert(norminf(shunt.g-new_shunt.g) < eps*(1+norminf(shunt.g)))
        # assert(norminf(shunt.b-new_shunt.b) < eps*(1+norminf(shunt.b)))
        assert(norminf(shunt.b_max-new_shunt.b_max) < eps*(1+norminf(shunt.b_max)))
        assert(norminf(shunt.b_min-new_shunt.b_min) < eps*(1+norminf(shunt.b_min)))

    # Loads
    assert(net1.num_loads == net2.num_loads)
    for i in range(net1.num_loads):
        load = net1.loads[i]
        new_load = net2.loads[i]
        assert(load is not new_load)
        assert(load.num_periods == new_load.num_periods)
        assert(load.bus.index == new_load.bus.index)
        assert(norminf(load.P-new_load.P) < eps)
        assert(norminf(load.P_max-new_load.P_max) < eps)
        assert(norminf(load.P_min-new_load.P_min) < eps)
        assert(norminf(load.Q-new_load.Q) < eps)
        assert(norminf(load.target_power_factor-new_load.target_power_factor) < eps)
        assert(norminf(load.util_coeff_Q0-new_load.util_coeff_Q0) < eps)
        assert(norminf(load.util_coeff_Q1-new_load.util_coeff_Q1) < eps)
        assert(norminf(load.util_coeff_Q2-new_load.util_coeff_Q2) < eps)

    # Batteries
    assert(net1.num_batteries == net2.num_batteries)
    for i in range(net1.num_batteries):
        bat = net1.batteries[i]
        new_bat = net2.batteries[i]
        assert(bat is not new_bat)
        assert(bat.num_periods == new_bat.num_periods)
        assert(bat.bus.index == new_bat.bus.index)
        assert(norminf(bat.P-new_bat.P) < eps)
        assert(norminf(bat.P_max-new_bat.P_max) < eps)
        assert(norminf(bat.P_min-new_bat.P_min) < eps)
        assert(norminf(bat.eta_c-new_bat.eta_c) < eps)
        assert(norminf(bat.eta_d-new_bat.eta_d) < eps)
        assert(norminf(bat.E-new_bat.E) < eps)
        assert(norminf(bat.E_init-new_bat.E_init) < eps)
        assert(norminf(bat.E_final-new_bat.E_final) < eps)
        assert(norminf(bat.E_max-new_bat.E_max) < eps)

    # Hashes
    for bus in net1.buses:
        assert(bus.index == net1.get_bus_by_number(bus.number).index)
        assert(bus.name == net1.get_bus_by_name(bus.name).name)
    for vargen in net1.var_generators:
        assert(vargen.index == net1.get_var_generator_by_name(vargen.name).index)



if __name__ == '__main__':
    main()
