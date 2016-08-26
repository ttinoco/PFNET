#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
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
    ctypedef double REAL
 
    void NET_add_vargens(Net* net, cbus.Bus* bus_list, REAL penetration, REAL uncertainty, REAL corr_radius, REAL corr_value)
    void NET_adjust_generators(Net* net)
    cbus.Bus* NET_bus_hash_number_find(Net* net, int number)
    cbus.Bus* NET_bus_hash_name_find(Net* net, char* name)
    cbus.Vargen* NET_vargen_hash_name_find(Net* net, char* name)
    void NET_clear_error(Net* net)
    void NET_clear_flags(Net* net)
    void NET_clear_properties(Net* net)
    void NET_clear_sensitivities(Net* net)
    cbus.Bus* NET_create_sorted_bus_list(Net* net, int sort_by, int t)
    cmat.Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr)
    void NET_del(Net* net)
    REAL NET_get_base_power(Net* net)
    cbus.Bus* NET_get_bus(Net* net, int index)
    cbranch.Branch* NET_get_branch(Net* net, int index)
    char* NET_get_error_string(Net* net)
    cgen.Gen* NET_get_gen(Net* net, int index)
    cshunt.Shunt* NET_get_shunt(Net* net, int index)
    cload.Load* NET_get_load(Net* net, int index)
    cvargen.Vargen* NET_get_vargen(Net* net, int index)
    cbat.Bat* NET_get_bat(Net* net, int index)
    cbus.Bus* NET_get_load_buses(Net* net)
    cbus.Bus* NET_get_gen_buses(Net* net)
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
    cmat.Mat* NET_get_var_projection(Net* net, char obj_type, char var)
    bint NET_has_error(Net* net)
    void NET_load(Net* net, char* filename, int output_level)
    Net* NET_new(int num_periods)
    void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, char val_mask)
    void NET_set_flags_of_component(Net* net, void* obj, char obj_type, char flag_mask, char val_mask)
    void NET_set_var_values(Net* net, cvec.Vec* values)
    void NET_show_components(Net* net)
    char* NET_get_show_components_str(Net* net)
    void NET_show_properties(Net* net, int t)
    char* NET_get_show_properties_str(Net* net, int t)
    void NET_show_buses(Net* net, int number, int sort_by, int t)
    void NET_update_properties(Net* net, cvec.Vec* values)
    void NET_update_set_points(Net* net)
    
     
          
