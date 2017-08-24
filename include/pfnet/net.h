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
#include "vector.h"
#include "matrix.h"
#include "utils.h"

// Controls
#define NET_CONTROL_EPS 1e-4      /**< @brief Safeguard for small control ranges (p.u.). */
#define NET_CONTROL_ACTION_PCT 2. /**< @brief Percent threshold for counting control actions (%). */

// Base power
#define NET_BASE_POWER 100 /**< @brief Default system base power (MVA). */

// Buffer
#define NET_BUFFER_SIZE 1024 /**< @brief Default network buffer size for strings */

// Net
typedef struct Net Net;

// Prototypes
/** @brief Adjust generator powers to obtain correct participations without affecting total injections. */
void NET_adjust_generators(Net* net);

void NET_add_vargens(Net* net, Bus* bus_list, REAL power_capacity, REAL power_base, REAL power_std, REAL corr_radius, REAL corr_value);
void NET_add_batteries(Net* net, Bus* bus_list, REAL power_capacity,  REAL energy_capacity, REAL eta_c, REAL eta_d);
void NET_bus_hash_number_add(Net* net, Bus* bus);
Bus* NET_bus_hash_number_find(Net* net, int number);
void NET_bus_hash_name_add(Net* net, Bus* bus);
Bus* NET_bus_hash_name_find(Net* net, char* name);
void NET_vargen_hash_name_add(Net* net, Vargen* vg);
Vargen* NET_vargen_hash_name_find(Net* net, char* name);
BOOL NET_check(Net* net, BOOL verbose);
void NET_clear_data(Net* net);
void NET_clear_error(Net* net);
void NET_clear_flags(Net* net);
void NET_clear_outages(Net* net);
void NET_clear_properties(Net* net);
void NET_clear_sensitivities(Net* net);
Bus* NET_create_sorted_bus_list(Net* net, int sort_by, int t);
Mat* NET_create_vargen_P_sigma(Net* net, int spread, REAL corr);
void NET_propagate_data_in_time(Net* net, int start, int end);
int NET_get_bus_neighbors(Net* net, Bus* bus, int spread, int* neighbors, char* queued);
void NET_del(Net* net);
void NET_init(Net* net, int num_periods);
REAL NET_get_base_power(Net* net);
Branch* NET_get_branch(Net* net, int index);
Bus* NET_get_bus(Net* net, int index);
Bus* NET_get_bus_hash_number(Net* net);
Bus* NET_get_bus_hash_name(Net* net);
Vargen* NET_get_vargen_hash_name(Net* net);
char* NET_get_error_string(Net* net);
Gen* NET_get_gen(Net* net, int index);
Load* NET_get_load(Net* net, int index);
Shunt* NET_get_shunt(Net* net, int index);
Vargen* NET_get_vargen(Net* net, int index);
Bat* NET_get_bat(Net* net, int index);
Bus* NET_get_gen_buses(Net* net);
Bus* NET_get_load_buses(Net* net);
int NET_get_num_periods(Net* net);
int NET_get_num_buses(Net* net);
int NET_get_num_slack_buses(Net* net);
int NET_get_num_buses_reg_by_gen(Net* net);
int NET_get_num_buses_reg_by_tran(Net* net);
int NET_get_num_buses_reg_by_tran_only(Net* net);
int NET_get_num_buses_reg_by_shunt(Net* net);
int NET_get_num_buses_reg_by_shunt_only(Net* net);
int NET_get_num_branches(Net* net);
int NET_get_num_branches_not_on_outage(Net* net);
int NET_get_num_fixed_trans(Net* net);
int NET_get_num_lines(Net* net);
int NET_get_num_phase_shifters(Net* net);
int NET_get_num_tap_changers(Net* net);
int NET_get_num_tap_changers_v(Net* net);
int NET_get_num_tap_changers_Q(Net* net);
int NET_get_num_gens(Net* net);
int NET_get_num_gens_not_on_outage(Net* net);
int NET_get_num_reg_gens(Net* net);
int NET_get_num_slack_gens(Net* net);
int NET_get_num_P_adjust_gens(Net* net);
int NET_get_num_loads(Net* net);
int NET_get_num_P_adjust_loads(Net* net);
int NET_get_num_shunts(Net* net);
int NET_get_num_fixed_shunts(Net* net);
int NET_get_num_switched_shunts(Net* net);
int NET_get_num_vargens(Net* net);
int NET_get_num_bats(Net* net);
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
int NET_get_num_actions(Net* net, int t);
REAL NET_get_vargen_corr_radius(Net* net);
REAL NET_get_vargen_corr_value(Net* net);
char* NET_get_json_string(Net* net);
BOOL NET_has_error(Net* net);
Net* NET_new(int num_periods);
void NET_set_base_power(Net* net, REAL base_power);
void NET_set_branch_array(Net* net, Branch* branch, int num);
void NET_set_bus_array(Net* net, Bus* bus, int num);
void NET_set_gen_array(Net* net, Gen* gen, int num);
void NET_set_load_array(Net* net, Load* load, int num);
void NET_set_shunt_array(Net* net, Shunt* shunt, int num);
void NET_set_vargen_array(Net* net, Vargen* gen, int num);
void NET_set_bat_array(Net* net, Bat* bat, int num);
void NET_set_flags(Net* net, char obj_type, char flag_mask, char prop_mask, unsigned char val_mask);
void NET_set_flags_of_component(Net* net, void* obj, char obj_type, char flag_mask, unsigned char val_mask);
void NET_set_var_values(Net* net, Vec* values);
void NET_set_vargen_buses(Net* net, Bus* bus_list);
void NET_set_bat_buses(Net* net, Bus* bus_list);
void NET_show_components(Net* net);
char* NET_get_show_components_str(Net* net);
void NET_show_properties(Net* net, int t);
char* NET_get_show_properties_str(Net* net, int t);
void NET_show_buses(Net* net, int number, int sort_by, int t);
void NET_update_properties_step(Net* net, Branch* br, int t, Vec* values);
void NET_update_properties(Net* net, Vec* values);
void NET_update_set_points(Net* net);

#endif
