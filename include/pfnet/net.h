/** @file net.h
 *  @brief This file lists the constants and routines associated with the Net data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __NET_HEADER__
#define __NET_HEADER__

#include <stdio.h>
#include "types.h"
#include "bus.h"
#include "branch.h"
#include "load.h"
#include "gen.h"
#include "shunt.h"
#include "vargen.h"
#include "bat.h"
#include "bus_dc.h"
#include "branch_dc.h"
#include "conv_vsc.h"
#include "conv_csc.h"
#include "facts.h"
#include "reg_obj.h"
#include "vector.h"
#include "matrix.h"
#include "utils.h"

// Base power
#define NET_BASE_POWER 100 /**< @brief Default system base power (MVA). */

// Buffer
#define NET_BUFFER_SIZE 1024 /**< @brief Default network buffer size for strings */

// Pre contingency status
#define PRE_CONT_UNSET -1         /**< @brief Pre contingency flag of the component is not set */
#define PRE_CONT_ONLINE 1         /**< @brief Component was online pre contingency */
#define PRE_CONT_OFFLINE 0        /**< @brief Component was offline pre contingency */

// Net
typedef struct Net Net;

// Prototypes
void NET_inc_state_tag(Net* net);
unsigned long int NET_get_state_tag(Net* net);

/** @brief Adjusts regulating and slack generators to balance bus power mismatches. */
void NET_adjust_generators(Net* net);

void NET_add_buses(Net* net, Bus** bus_ptr_array, int size);
void NET_del_buses(Net* net, Bus** bus_ptr_array, int size);

void NET_add_branches(Net* net, Branch** br_ptr_array, int size);
void NET_del_branches(Net* net, Branch** br_ptr_array, int size);

void NET_add_gens(Net* net, Gen** gen_ptr_array, int size);
void NET_del_gens(Net* net, Gen** gen_ptr_array, int size);

void NET_add_loads(Net* net, Load** load_ptr_array, int size);
void NET_del_loads(Net* net, Load** load_ptr_array, int size);

void NET_add_shunts(Net* net, Shunt** shunt_ptr_array, int size);
void NET_del_shunts(Net* net, Shunt** shunt_ptr_array, int size);

void NET_add_bats(Net* net, Bat** bat_ptr_array, int size);
void NET_del_bats(Net* net, Bat** bat_ptr_array, int size);

void NET_add_vargens(Net* net, Vargen** vargen_ptr_array, int size);
void NET_del_vargens(Net* net, Vargen** vargen_ptr_array, int size);

void NET_add_vsc_convs(Net* net, ConvVSC** conv_ptr_array, int size);
void NET_del_vsc_convs(Net* net, ConvVSC** conv_ptr_array, int size);

void NET_add_csc_convs(Net* net, ConvCSC** conv_ptr_array, int size);
void NET_del_csc_convs(Net* net, ConvCSC** conv_ptr_array, int size);

void NET_add_dc_buses(Net* net, BusDC** bus_ptr_array, int size);
void NET_del_dc_buses(Net* net, BusDC** bus_ptr_array, int size);

void NET_add_dc_branches(Net* net, BranchDC** br_ptr_array, int size);
void NET_del_dc_branches(Net* net, BranchDC** br_ptr_array, int size);

void NET_add_facts(Net* net, Facts** facts_ptr_array, int size);
void NET_del_facts(Net* net, Facts** facts_ptr_array, int size);

void NET_add_vargens_from_params(Net* net, Bus* bus_list, REAL power_capacity, REAL power_base, REAL power_std, REAL corr_radius, REAL corr_value);
void NET_add_batteries_from_params(Net* net, Bus* bus_list, REAL power_capacity,  REAL energy_capacity, REAL eta_c, REAL eta_d);

void NET_add_red_bus(Net* net, Bus* bus);

void NET_bus_hash_number_add(Net* net, Bus* bus);
Bus* NET_bus_hash_number_find(Net* net, int number);
void NET_bus_hash_name_add(Net* net, Bus* bus);
Bus* NET_bus_hash_name_find(Net* net, char* name);
void NET_dc_bus_hash_number_add(Net* net, BusDC* bus);
BusDC* NET_dc_bus_hash_number_find(Net* net, int number);
void NET_dc_bus_hash_name_add(Net* net, BusDC* bus);
BusDC* NET_dc_bus_hash_name_find(Net* net, char* name);
BOOL NET_check(Net* net, BOOL verbose);
void NET_clear_data(Net* net);
void NET_clear_error(Net* net);
void NET_clear_flags(Net* net);
void NET_clear_properties(Net* net);
void NET_clear_sensitivities(Net* net);
Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr);
void NET_copy_from_net(Net* net, Net* other, int* bus_index_map, int* branch_index_map, int mode);
void NET_del(Net* net);
Net* NET_extract_subnet(Net* net, Bus** bus_ptr_array, int size);
void NET_init(Net* net, int num_periods);
Net* NET_get_copy(Net* net, BOOL merge_buses);
int NET_get_bus_neighbors(Net* net, Bus* bus, int spread, int* neighbors, char* queued);
REAL NET_get_base_power(Net* net);
Branch* NET_get_branch(Net* net, int index);
Bus* NET_get_bus(Net* net, int index);
Bus* NET_get_bus_hash_number(Net* net);
Bus* NET_get_bus_hash_name(Net* net);
BusDC* NET_get_dc_bus_hash_number(Net* net);
BusDC* NET_get_dc_bus_hash_name(Net* net);
char* NET_get_error_string(Net* net);
Gen* NET_get_gen(Net* net, int index);
Load* NET_get_load(Net* net, int index);
Shunt* NET_get_shunt(Net* net, int index);
Vargen* NET_get_vargen(Net* net, int index);
Bat* NET_get_bat(Net* net, int index);
ConvVSC* NET_get_vsc_conv(Net* net, int index);
ConvCSC* NET_get_csc_conv(Net* net, int index);
BusDC* NET_get_dc_bus(Net* net, int index);
BranchDC* NET_get_dc_branch(Net* net, int index);
Facts* NET_get_facts(Net* net, int index);
Bus* NET_get_gen_buses(Net* net);
Bus* NET_get_load_buses(Net* net);

Gen* NET_get_gen_from_name_and_bus_number(Net* net, char* name, int number);
Branch* NET_get_branch_from_name_and_bus_numbers(Net* net, char* name, int number1, int number2);
Shunt* NET_get_shunt_from_name_and_bus_number(Net* net, char* name, int number);
Shunt* NET_get_fixed_shunt_from_name_and_bus_number(Net* net, char* name, int number);
Shunt* NET_get_switched_shunt_from_name_and_bus_number(Net* net, char* name, int number);
Load* NET_get_load_from_name_and_bus_number(Net* net, char* name, int number);
Vargen* NET_get_vargen_from_name_and_bus_number(Net* net, char* name, int number);
Bat* NET_get_bat_from_name_and_bus_number(Net* net, char* name, int number);
ConvCSC* NET_get_csc_conv_from_name_and_ac_bus_number(Net* net, char* name, int number);
ConvCSC* NET_get_csc_conv_from_name_and_dc_bus_number(Net* net, char* name, int number);
ConvCSC* NET_get_csc_conv_from_name_and_dc_bus_name(Net* net, char* name, char* bus_name);
ConvVSC* NET_get_vsc_conv_from_name_and_ac_bus_number(Net* net, char* name, int number);
ConvVSC* NET_get_vsc_conv_from_name_and_dc_bus_number(Net* net, char* name, int number);
ConvVSC* NET_get_vsc_conv_from_name_and_dc_bus_name(Net* net, char* name, char* bus_name);
BranchDC* NET_get_dc_branch_from_name_and_dc_bus_numbers(Net* net, char* name, int number1, int number2);
BranchDC* NET_get_dc_branch_from_name_and_dc_bus_names(Net* net, char* name, char* bus1_name, char* bus2_name);
Facts* NET_get_facts_from_name_and_bus_numbers(Net* net, char* name, int number1, int number2);

int NET_get_num_buses(Net* net, BOOL only_in_service);
int NET_get_num_buses_out_of_service(Net* net);
int NET_get_num_slack_buses(Net* net, BOOL only_in_service);
int NET_get_num_star_buses(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_gen(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_tran(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_tran_only(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_shunt(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_shunt_only(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_vsc_conv(Net* net, BOOL only_in_service);
int NET_get_num_buses_reg_by_facts(Net* net, BOOL only_in_service);
int NET_get_num_red_buses(Net* net);
int NET_get_num_branches(Net* net, BOOL only_in_service);
int NET_get_num_branches_out_of_service(Net* net);
int NET_get_num_fixed_trans(Net* net, BOOL only_in_service);
int NET_get_num_lines(Net* net, BOOL only_in_service);
int NET_get_num_zero_impedance_lines(Net* net, BOOL only_in_service);
int NET_get_num_phase_shifters(Net* net, BOOL only_in_service);
int NET_get_num_tap_changers(Net* net, BOOL only_in_service);
int NET_get_num_tap_changers_v(Net* net, BOOL only_in_service);
int NET_get_num_tap_changers_Q(Net* net, BOOL only_in_service);
int NET_get_num_gens(Net* net, BOOL only_in_service);
int NET_get_num_gens_out_of_service(Net* net);
int NET_get_num_reg_gens(Net* net, BOOL only_in_service);
int NET_get_num_slack_gens(Net* net, BOOL only_in_service);
int NET_get_num_P_adjust_gens(Net* net, BOOL only_in_service);
int NET_get_num_loads(Net* net, BOOL only_in_service);
int NET_get_num_loads_out_of_service(Net* net);
int NET_get_num_P_adjust_loads(Net* net, BOOL only_in_service);
int NET_get_num_vdep_loads(Net* ne, BOOL only_in_servicet);
int NET_get_num_shunts(Net* net, BOOL only_in_service);
int NET_get_num_shunts_out_of_service(Net* net);
int NET_get_num_fixed_shunts(Net* net, BOOL only_in_service);
int NET_get_num_switched_shunts(Net* net, BOOL only_in_service);
int NET_get_num_switched_v_shunts(Net* net, BOOL only_in_service);
int NET_get_num_vargens(Net* net, BOOL only_in_service);
int NET_get_num_vargens_out_of_service(Net* net);
int NET_get_num_bats(Net* net, BOOL only_in_service);
int NET_get_num_bats_out_of_service(Net* net);
int NET_get_num_csc_convs(Net* net, BOOL only_in_service);
int NET_get_num_csc_convs_out_of_service(Net* net);
int NET_get_num_vsc_convs(Net* net, BOOL only_in_service);
int NET_get_num_vsc_convs_out_of_service(Net* net);
int NET_get_num_vsc_convs_in_P_dc_mode(Net* net, BOOL only_in_service);
int NET_get_num_vsc_convs_in_v_dc_mode(Net* net, BOOL only_in_service);
int NET_get_num_vsc_convs_in_v_ac_mode(Net* net, BOOL only_in_service);
int NET_get_num_vsc_convs_in_f_ac_mode(Net* net, BOOL only_in_service);
int NET_get_num_dc_buses(Net* net, BOOL only_in_service);
int NET_get_num_dc_buses_out_of_service(Net* net);
int NET_get_num_dc_branches(Net* net, BOOL only_in_service);
int NET_get_num_dc_branches_out_of_service(Net* net);
int NET_get_num_facts(Net* net, BOOL only_in_service);
int NET_get_num_facts_out_of_service(Net* net);
int NET_get_num_facts_in_normal_series_mode(Net* net, BOOL only_in_service);
int NET_get_num_reg_facts(Net* net, BOOL only_in_service);

int NET_get_num_periods(Net* net);
int NET_get_num_vars(Net* net);
int NET_get_num_fixed(Net* net);
int NET_get_num_bounded(Net* net);
int NET_get_num_sparse(Net* net);
REAL NET_get_total_gen_P(Net* net, int t);
REAL NET_get_total_gen_Q(Net* net, int t);
REAL NET_get_total_load_P(Net* net, int t);
REAL NET_get_total_load_Q(Net* net, int t);
Vec* NET_get_var_values(Net* net, int code);
char* NET_get_var_info_string(Net* net, int index);
Mat* NET_get_var_projection(Net* net, char obj_type, char prop_mask, unsigned char var, int t_start, int t_end);
REAL NET_get_bus_v_max(Net* net, int t);
REAL NET_get_bus_v_min(Net* net, int t);
REAL NET_get_bus_v_vio(Net* net, int t);
REAL NET_get_bus_P_mis(Net* net, int t);
REAL NET_get_bus_Q_mis(Net* net, int t);
REAL NET_get_gen_P_cost(Net* net, int t);
REAL NET_get_gen_v_dev(Net* net, int t);
REAL NET_get_gen_Q_vio(Net* net, int t);
REAL NET_get_gen_P_vio(Net* net, int t);
REAL NET_get_tran_v_vio(Net* net, int t);
REAL NET_get_tran_r_vio(Net* net, int t);
REAL NET_get_tran_p_vio(Net* net, int t);
REAL NET_get_shunt_v_vio(Net* net, int t);
REAL NET_get_shunt_b_vio(Net* net, int t);
REAL NET_get_load_P_util(Net* net, int t);
REAL NET_get_load_P_vio(Net* net, int t);
REAL NET_get_vargen_corr_radius(Net* net);
REAL NET_get_vargen_corr_value(Net* net);
char* NET_get_json_string(Net* net);
BOOL NET_has_error(Net* net);
Net* NET_new(int num_periods);
void NET_make_all_in_service(Net* net);
char* NET_mark_reachable_dc_buses(Net* net, BusDC* dc_bus);
void NET_propagate_data_in_time(Net* net, int start, int end);
int NET_round_discrete_switched_shunts_b(Net* net, int t);
void NET_clip_switched_shunts_b(Net* net, int t);

void NET_set_base_power(Net* net, REAL base_power);
void NET_set_branch_array(Net* net, Branch* branch, int num);
void NET_set_bus_array(Net* net, Bus* bus, int num);
void NET_set_gen_array(Net* net, Gen* gen, int num);
void NET_set_load_array(Net* net, Load* load, int num);
void NET_set_shunt_array(Net* net, Shunt* shunt, int num);
void NET_set_vargen_array(Net* net, Vargen* gen, int num);
void NET_set_bat_array(Net* net, Bat* bat, int num);
void NET_set_vsc_conv_array(Net* net, ConvVSC* conv, int num);
void NET_set_csc_conv_array(Net* net, ConvCSC* conv, int num);
void NET_set_dc_bus_array(Net* net, BusDC* bus, int num);
void NET_set_dc_branch_array(Net* net, BranchDC* branch, int num);
void NET_set_facts_array(Net* net, Facts* facts, int num);
void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, unsigned char val_mask);
void NET_set_flags_of_component(Net* net, void* obj, char obj_type, char flag_mask, unsigned char val_mask);
void NET_set_var_values(Net* net, Vec* values);
void NET_set_vargen_buses(Net* net, Bus* bus_list);
void NET_set_bat_buses(Net* net, Bus* bus_list);
void NET_set_equiv_buses(Net* net);
void NET_show_components(Net* net, int output_level);
char* NET_get_show_components_str(Net* net, int output_level);
void NET_show_properties(Net* net, int t);
void NET_show_equiv_buses(Net* net);
void NET_show_red_buses(Net* net);
char* NET_get_show_properties_str(Net* net, int t);
void NET_update_properties_step(Net* net, Bus* bus, BusDC* busdc, int t, Vec* values);
void NET_update_properties(Net* net, Vec* values);
void NET_update_reg_Q_participations(Net* net, int t);
void NET_update_set_points(Net* net);
void NET_update_hash_tables(Net* net);
void NET_localize_gen_regulation(Net* net, int max_dist);

#endif
