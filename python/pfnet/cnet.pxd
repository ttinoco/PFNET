#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cvec
cimport cmat
cimport cbus
cimport cbranch
cimport cgen
cimport cload
cimport cshunt
cimport cvargen
cimport cbat

cdef extern from "pfnet/net.h":

    ctypedef struct Net
    ctypedef struct Bus
    ctypedef double REAL
 
    void NET_add_vargens(Net* net, cbus.Bus* bus_list, REAL power_capacity, REAL power_base, REAL power_std, REAL corr_radius, REAL corr_value)
    void NET_add_batteries(Net* net, cbus.Bus* bus_list, REAL power_capacity,  REAL energy_capacity, REAL eta_c, REAL eta_d)        
    void NET_adjust_generators(Net* net)
    cbus.Bus* NET_bus_hash_number_find(Net* net, int number)
    cbus.Bus* NET_bus_hash_name_find(Net* net, char* name)
    void NET_bus_hash_number_add(Net* net, cbus.Bus* bus)
    void NET_bus_hash_name_add(Net* net, cbus.Bus* bus)
    void NET_clear_error(Net* net)
    void NET_clear_flags(Net* net)
    void NET_clear_properties(Net* net)
    void NET_clear_sensitivities(Net* net)
    cmat.Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr)
    void NET_copy_from_net(Net* net, Net* other_net)
    void NET_del(Net* net)
    Net* NET_get_copy(Net* net)
    REAL NET_get_base_power(Net* net)
    char* NET_get_error_string(Net* net)

    cbus.Bus* NET_get_bus(Net* net, int index)
    cbranch.Branch* NET_get_branch(Net* net, int index)
    cgen.Gen* NET_get_gen(Net* net, int index)
    cshunt.Shunt* NET_get_shunt(Net* net, int index)
    cload.Load* NET_get_load(Net* net, int index)
    cvargen.Vargen* NET_get_vargen(Net* net, int index)
    cbat.Bat* NET_get_bat(Net* net, int index)
    cbus.Bus* NET_get_load_buses(Net* net)
    cbus.Bus* NET_get_gen_buses(Net* net)

    cgen.Gen* NET_get_gen_from_name_and_bus_number(Net* net, char* name, int number)
    cbranch.Branch* NET_get_branch_from_name_and_bus_numbers(Net* net, char* name, int number1, int number2)
    cshunt.Shunt* NET_get_shunt_from_name_and_bus_number(Net* net, char* name, int number)
    cload.Load* NET_get_load_from_name_and_bus_number(Net* net, char* name, int number)
    cvargen.Vargen* NET_get_vargen_from_name_and_bus_number(Net* net, char* name, int number)
    cbat.Bat* NET_get_bat_from_name_and_bus_number(Net* net, char* name, int number)

    REAL NET_get_total_load_P(Net* net, int t)
    int NET_get_num_periods(Net* net)
    int NET_get_num_buses(Net* net)
    int NET_get_num_slack_buses(Net* net)
    int NET_get_num_buses_reg_by_gen(Net* net)
    int NET_get_num_buses_reg_by_tran(Net* net)
    int NET_get_num_buses_reg_by_tran_only(Net* net)
    int NET_get_num_buses_reg_by_shunt(Net* net)
    int NET_get_num_buses_reg_by_shunt_only(Net* net)
    int NET_get_num_branches(Net* net)
    int NET_get_num_branches_not_on_outage(Net* net)
    int NET_get_num_fixed_trans(Net* net)
    int NET_get_num_lines(Net* net)
    int NET_get_num_phase_shifters(Net* net)
    int NET_get_num_tap_changers(Net* net)
    int NET_get_num_tap_changers_v(Net* net)
    int NET_get_num_tap_changers_Q(Net* net)
    int NET_get_num_gens(Net* net)
    int NET_get_num_gens_not_on_outage(Net* net)
    int NET_get_num_reg_gens(Net* net)
    int NET_get_num_slack_gens(Net* net)
    int NET_get_num_P_adjust_gens(Net* net)
    int NET_get_num_loads(Net* net)
    int NET_get_num_P_adjust_loads(Net* net)
    int NET_get_num_shunts(Net* net)
    int NET_get_num_fixed_shunts(Net* net)
    int NET_get_num_switched_shunts(Net* net)
    int NET_get_num_vargens(Net* net)
    int NET_get_num_bats(Net* net)
    int NET_get_num_vars(Net* net)
    int NET_get_num_fixed(Net* net)
    int NET_get_num_bounded(Net* net)
    int NET_get_num_sparse(Net* net)
    REAL NET_get_bus_v_max(Net* net, int t)
    REAL NET_get_bus_v_min(Net* net, int t)
    REAL NET_get_bus_v_vio(Net* net, int t)
    REAL NET_get_bus_P_mis(Net* net, int t)
    REAL NET_get_bus_Q_mis(Net* net, int t)
    REAL NET_get_gen_P_cost(Net* net, int t)
    REAL NET_get_gen_v_dev(Net* net, int t)
    REAL NET_get_gen_Q_vio(Net* net, int t)
    REAL NET_get_gen_P_vio(Net* net, int t)
    REAL NET_get_tran_v_vio(Net* net, int t)
    REAL NET_get_tran_r_vio(Net* net, int t)
    REAL NET_get_tran_p_vio(Net* net, int t)
    REAL NET_get_shunt_v_vio(Net* net, int t)
    REAL NET_get_shunt_b_vio(Net* net, int t)
    REAL NET_get_load_P_util(Net* net, int t)
    REAL NET_get_load_P_vio(Net* net, int t)
    int NET_get_num_actions(Net* net, int t)
    REAL NET_get_vargen_corr_radius(Net* net)
    REAL NET_get_vargen_corr_value(Net* net)
    cvec.Vec* NET_get_var_values(Net* net, int code)
    char* NET_get_var_info_string(Net* net, int index)
    cmat.Mat* NET_get_var_projection(Net* net, char obj_type, char prop_mask, char var, int t_start, int t_end)
    char* NET_get_json_string(Net* net)
    bint NET_has_error(Net* net)
    Net* NET_new(int num_periods)
    void NET_set_base_power(Net* net, REAL base_power)
    void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, char val_mask)
    void NET_set_flags_of_component(Net* net, void* obj, char obj_type, char flag_mask, char val_mask)
    void NET_set_var_values(Net* net, cvec.Vec* values)
    void NET_show_components(Net* net)
    char* NET_get_show_components_str(Net* net)
    void NET_show_properties(Net* net, int t)
    char* NET_get_show_properties_str(Net* net, int t)
    void NET_update_properties(Net* net, cvec.Vec* values)
    void NET_propagate_data_in_time(Net* net, int start, int end)
    void NET_update_set_points(Net* net)

    void NET_set_bus_array(Net* net, cbus.Bus* bus_list, int num_buses)
    void NET_set_branch_array(Net* net, cbranch.Branch* branch_list, int num_branches)
    void NET_set_gen_array(Net* net, cgen.Gen* gen_list, int num_generators)
    void NET_set_load_array(Net* net, cload.Load* load_list, int num_loads)
    void NET_set_shunt_array(Net* net, cshunt.Shunt* shunt_list, int num_shunts)
    void NET_set_vargen_array(Net* net, cvargen.Vargen* vargen_list, int num_vargens)
    void NET_set_bat_array(Net* net, cbat.Bat* bat_list, int num_batteries)
