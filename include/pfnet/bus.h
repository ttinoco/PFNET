/** @file bus.h
 *  @brief This file lists the constants and routines associated with the Bus data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BUS_HEADER__
#define __BUS_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "node.h"
#include "vector.h"
#include "uthash.h"

// Limits
#define BUS_DEFAULT_V_MAX 1.1 /**< @brief Default maximum voltage magnitude. */
#define BUS_DEFAULT_V_MIN 0.9 /**< @brief Default minimum voltage magnitude. */

// Infinity
#define BUS_INF_V_MAG 1e8 /**< @brief Infinite voltage magnitude (p.u.) */
#define BUS_INF_V_ANG 1e2 /**< @brief Infinite voltage angle (radians) */

// Variables
/** \defgroup bus_vars Bus Variable Masks
 *  @{
 */
#define BUS_VAR_VMAG 0x01 /**< @brief Variable: voltage magnitude. */
#define BUS_VAR_VANG 0x02 /**< @brief Variable: volatge angle. */
/** @} */

// Properties
/** \defgroup bus_props Bus Property Masks
 *  @{
 */
#define BUS_PROP_ANY 0x00            /**< @brief Property: any */
#define BUS_PROP_SLACK 0x01          /**< @brief Property: slack bus */
#define BUS_PROP_VSET_REG 0x02       /**< @brief Property: voltage set point regulated */
#define BUS_PROP_REG_BY_GEN 0x04     /**< @brief Property: bus regulated by generator */
#define BUS_PROP_REG_BY_TRAN 0x08    /**< @brief Property: bus regulated by transformer */
#define BUS_PROP_REG_BY_SHUNT 0x10   /**< @brief Property: bus regualted by shunt device */
#define BUS_PROP_NOT_REG_BY_GEN 0x20 /**< @brief Property: bus not regulated by generator */
#define BUS_PROP_NOT_SLACK 0x40      /**< @brief Property: non-slack bus */
/** @} */

// Sensitivities
/** \defgroup bus_sens Bus Sensitivities
 *  @{
 */
#define BUS_SENS_LARGEST 0
#define BUS_SENS_P_BALANCE 1      /**< @brief Sensitivity: active power balance */
#define BUS_SENS_Q_BALANCE 2      /**< @brief Sensitivity: reactive power balance */
#define BUS_SENS_V_MAG_U_BOUND 3  /**< @brief Sensitivity: voltage magnitude upper bound */
#define BUS_SENS_V_MAG_L_BOUND 4  /**< @brief Sensitivity: voltage magnitude lower bound */
#define BUS_SENS_V_ANG_U_BOUND 5  /**< @brief Sensitivity: voltage angle upper bound */
#define BUS_SENS_V_ANG_L_BOUND 6  /**< @brief Sensitivity: voltage angle lower bound */
#define BUS_SENS_V_SET_REG 7      /**< @brief Sensitivity: voltage set point regulation */
#define BUS_SENS_V_REG_BY_TRAN 8  /**< @brief Sensitivity: voltage magnitude regulation by transformer */
#define BUS_SENS_V_REG_BY_SHUNT 9 /**< @brief Sensitivity: voltage magnitude regulation by shunt device */
/** @} */

// Mismatches
/** \defgroup bus_mis Bus Power Mismatches
 *  @{
 */
#define BUS_MIS_LARGEST 10
#define BUS_MIS_ACTIVE 11   /**< @brief Mismatch: active power */
#define BUS_MIS_REACTIVE 12 /**< @brief Mismatch: reactive power */
/** @} */

// Constants
/** \defgroup bus_const Bus Constants
 *  @{
 */
#define BUS_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define BUS_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define BUS_NUM_JSON_FIELDS 45   /**< @brief Constant: max number of json fields */
/** @} */

// Bus
typedef struct Bus Bus;

// Other
typedef struct Gen Gen;
typedef struct Load Load;
typedef struct Branch Branch;
typedef struct Shunt Shunt;
typedef struct Vargen Vargen;
typedef struct Bat Bat;
typedef struct ConvVSC ConvVSC;
typedef struct ConvCSC ConvCSC;
typedef struct Facts Facts;
typedef struct Vec Vec;
typedef struct Net Net;

// Prototypes
/** @brief Adds generator to list of generators connected to bus. */
void BUS_add_gen(Bus* bus, Gen* gen);
void BUS_del_gen(Bus* bus, Gen* gen);

/** @brief Adds load to list of loads connected to bus. */
void BUS_add_load(Bus* bus, Load* load);
void BUS_del_load(Bus* bus, Load* load);

/** @brief Adds generator to list of generators regulating bus voltage. */
void BUS_add_reg_gen(Bus* bus, Gen* reg_gen);
void BUS_del_reg_gen(Bus* bus, Gen* reg_gen);

/** @brief Adds transformer to list of transformers regulating bus voltage. */
void BUS_add_reg_tran(Bus* bus, Branch* reg_tran);
void BUS_del_reg_tran(Bus* bus, Branch* reg_tran);

/** @brief Adds switched shunt to list of shunts regulating bus voltage. */
void BUS_add_reg_shunt(Bus* bus, Shunt* reg_shunt);
void BUS_del_reg_shunt(Bus* bus, Shunt* reg_shunt);

/** @brief Adds shunt to list of shunts connected to bus. */
void BUS_add_shunt(Bus* bus, Shunt* shunt);
void BUS_del_shunt(Bus* bus, Shunt* shunt);

/** @brief Adds variable generator to list of variable generators connected to bus. */
void BUS_add_vargen(Bus* bus, Vargen* gen);
void BUS_del_vargen(Bus* bus, Vargen* gen);

/** @brief Adds battery to list of batteries connected to bus. */
void BUS_add_bat(Bus* bus, Bat* bat);
void BUS_del_bat(Bus* bus, Bat* bat);

/** @brief Adds branch to list of branches connected at "k" bus. */
void BUS_add_branch_k(Bus* bus, Branch* branch);

/** @brief Deletes branch from list of branches connected at "k" bus. */
void BUS_del_branch_k(Bus* bus, Branch* branch);

/** @brief Adds branch to list of branches connected at "m" bus. */
void BUS_add_branch_m(Bus* bus, Branch* branch);

/** @brief Deletes branch from list of branches connected at "m" bus. */
void BUS_del_branch_m(Bus* bus, Branch* branch);

/** @brief Adds VSC converter to list of VSC converters connected to bus. */
void BUS_add_vsc_conv(Bus* bus, ConvVSC* conv);
void BUS_del_vsc_conv(Bus* bus, ConvVSC* conv);

/** @brief Adds VSC converter to list of VSC converters regulating bus. */
void BUS_add_reg_vsc_conv(Bus* bus, ConvVSC* conv);
void BUS_del_reg_vsc_conv(Bus* bus, ConvVSC* conv);

/** @brief Adds CSC converter to list of CSC converters connected to bus. */
void BUS_add_csc_conv(Bus* bus, ConvCSC* conv);
void BUS_del_csc_conv(Bus* bus, ConvCSC* conv);

void BUS_add_facts_k(Bus* bus, Facts* facts);
void BUS_del_facts_k(Bus* bus, Facts* facts);
void BUS_add_facts_m(Bus* bus, Facts* facts);
void BUS_del_facts_m(Bus* bus, Facts* facts);
void BUS_add_reg_facts(Bus* bus, Facts* facts);
void BUS_del_reg_facts(Bus* bus, Facts* facts);

void BUS_del_all_connections(Bus* bus);

void BUS_array_del(Bus* bus_array, int size);
BOOL BUS_array_check(Bus* bus_array, int size, BOOL verbose);
void* BUS_array_get(void* bus_array, int index);
Bus* BUS_array_new(int size, int num_periods);
void BUS_array_show(Bus* bus_array, int size, int t);
void BUS_array_get_max_mismatches(Bus* bus_array, int size, REAL* P, REAL* Q, int t);
BOOL BUS_check(Bus* bus, BOOL verbose);
void BUS_clear_flags(Bus* bus, char flag_type);
void BUS_clear_sensitivities(Bus* bus);
void BUS_clear_mismatches(Bus* bus);
void BUS_copy_from_bus(Bus* bus, Bus* other, int mode, BOOL propagate);

char BUS_get_flags_vars(Bus* bus);
char BUS_get_flags_fixed(Bus* bus);
char BUS_get_flags_bounded(Bus* bus);
char BUS_get_flags_sparse(Bus* bus);

int BUS_get_area(Bus* bus);
int BUS_get_zone(Bus* bus);
char BUS_get_obj_type(void* bus);
int BUS_get_degree(Bus* bus);
REAL BUS_get_price(Bus* bus, int t);
REAL* BUS_get_price_array(Bus* bus);
int BUS_get_index(Bus* bus);
int BUS_get_index_t(Bus* bus, int t);
int BUS_get_index_v_mag(Bus* bus, int t);
int* BUS_get_index_v_mag_array(Bus* bus);
int BUS_get_index_v_ang(Bus* bus, int t);
int* BUS_get_index_v_ang_array(Bus* bus);
int BUS_get_index_P(Bus* bus, int t);
int BUS_get_index_Q(Bus* bus, int t);
Bus* BUS_get_next(Bus* bus);
Node* BUS_get_equiv(Bus* bus);
int BUS_get_number(Bus* bus);
char* BUS_get_name(Bus* bus);
int BUS_get_alt_number(Bus* bus);
char* BUS_get_alt_name(Bus* bus);
int BUS_get_num_periods(Bus* bus);
int BUS_get_num_gens(Bus* bus);
int BUS_get_num_loads(Bus* bus);
int BUS_get_num_shunts(Bus* bus);
int BUS_get_num_vargens(Bus* bus);
int BUS_get_num_bats(Bus* bus);
int BUS_get_num_csc_convs(Bus* bus);
int BUS_get_num_vsc_convs(Bus* bus);
int BUS_get_num_reg_gens(Bus* bus);
int BUS_get_num_reg_trans(Bus* bus);
int BUS_get_num_reg_shunts(Bus* bus);
int BUS_get_num_reg_vsc_convs(Bus* bus);
int BUS_get_num_reg_facts(Bus* bus);
Gen* BUS_get_gen(Bus* bus);
Load* BUS_get_load(Bus* bus);
Gen* BUS_get_reg_gen(Bus* bus);
Branch* BUS_get_reg_tran(Bus* bus);
Shunt* BUS_get_reg_shunt(Bus* bus);
Shunt* BUS_get_shunt(Bus* bus);
Branch* BUS_get_branch_k(Bus* bus);
Branch* BUS_get_branch_m(Bus* bus);
Vargen* BUS_get_vargen(Bus* bus);
Bat* BUS_get_bat(Bus* bus);
ConvCSC* BUS_get_csc_conv(Bus* bus);
ConvVSC* BUS_get_vsc_conv(Bus* bus);
ConvVSC* BUS_get_reg_vsc_conv(Bus* bus);
Facts* BUS_get_facts_k(Bus* bus);
Facts* BUS_get_facts_m(Bus* bus);
Facts* BUS_get_reg_facts(Bus* bus);
REAL BUS_get_P_mis(Bus* bus, int t);
REAL BUS_get_Q_mis(Bus* bus, int t);
REAL BUS_get_total_gen_P(Bus* bus, int t);
REAL BUS_get_total_gen_Q(Bus* bus, int t);
REAL BUS_get_total_gen_Q_max(Bus* bus);
REAL BUS_get_total_gen_Q_min(Bus* bus);
REAL BUS_get_total_reg_gen_Q(Bus* bus, int t);
REAL BUS_get_total_reg_gen_Q_max(Bus* bus);
REAL BUS_get_total_reg_gen_Q_min(Bus* bus);
REAL BUS_get_total_load_P(Bus* bus, int t);
REAL BUS_get_total_load_Q(Bus* bus, int t);
REAL BUS_get_total_shunt_g(Bus* bus);
REAL BUS_get_total_shunt_b(Bus* bus, int t);
REAL BUS_get_v_base(Bus* bus);
REAL BUS_get_v_mag(Bus* bus, int t);
REAL* BUS_get_v_mag_array(Bus* bus);
REAL BUS_get_v_ang(Bus* bus, int t);
REAL* BUS_get_v_ang_array(Bus* bus);
REAL BUS_get_v_set(Bus* bus, int t);
REAL* BUS_get_v_set_array(Bus* bus);
REAL BUS_get_v_max_reg(Bus* bus);
REAL BUS_get_v_min_reg(Bus* bus);
REAL BUS_get_v_max_norm(Bus* bus);
REAL BUS_get_v_min_norm(Bus* bus);
REAL BUS_get_v_max_emer(Bus* bus);
REAL BUS_get_v_min_emer(Bus* bus);
void BUS_get_var_values(Bus* bus, Vec* values, int code);
char* BUS_get_var_info_string(Bus* bus, int index);
int BUS_get_num_vars(void* bus, unsigned char var, int t_start, int t_end);
Vec* BUS_get_var_indices(void* bus, unsigned char var, int t_start, int t_end);

REAL BUS_get_sens_P_balance(Bus* bus, int t);
REAL* BUS_get_sens_P_balance_array(Bus* bus);
REAL BUS_get_sens_Q_balance(Bus* bus, int t);
REAL* BUS_get_sens_Q_balance_array(Bus* bus);
REAL BUS_get_sens_v_mag_u_bound(Bus* bus, int t);
REAL* BUS_get_sens_v_mag_u_bound_array(Bus* bus);
REAL BUS_get_sens_v_mag_l_bound(Bus* bus, int t);
REAL* BUS_get_sens_v_mag_l_bound_array(Bus* bus);
REAL BUS_get_sens_v_ang_u_bound(Bus* bus, int t);
REAL* BUS_get_sens_v_ang_u_bound_array(Bus* bus);
REAL BUS_get_sens_v_ang_l_bound(Bus* bus, int t);
REAL* BUS_get_sens_v_ang_l_bound_array(Bus* bus);
REAL BUS_get_sens_v_set_reg(Bus* bus, int t);
REAL* BUS_get_sens_v_set_reg_array(Bus* bus);
REAL BUS_get_sens_v_reg_by_tran(Bus* bus, int t);
REAL* BUS_get_sens_v_reg_by_tran_array(Bus* bus);
REAL BUS_get_sens_v_reg_by_shunt(Bus* bus, int t);
REAL* BUS_get_sens_v_reg_by_shunt_array(Bus* bus);

REAL BUS_get_largest_sens(Bus* bus, int t);
int BUS_get_largest_sens_type(Bus* bus, int t);
REAL BUS_get_largest_mis(Bus* bus, int t);
int BUS_get_largest_mis_type(Bus* bus, int t);
REAL BUS_get_quantity(Bus* bus, int qtype, int t);
char* BUS_get_json_string(Bus* bus, char* output);
BOOL BUS_has_flags(void* bus, char flag_type, unsigned char mask);
BOOL BUS_has_properties(void* bus, char prop);
Bus* BUS_hash_number_add(Bus* bus_hash, Bus* bus);
void BUS_hash_number_del(Bus* bus_hash);
Bus* BUS_hash_number_find(Bus* bus_hash, int number);
int BUS_hash_number_len(Bus* bus_hash);
Bus* BUS_hash_name_add(Bus* bus_hash, Bus* bus);
void BUS_hash_name_del(Bus* bus_hash);
Bus* BUS_hash_name_find(Bus* bus_hash, char* name);
int BUS_hash_name_len(Bus* bus_hash);
void BUS_init(Bus* bus, int num_periods);
void BUS_inject_P(Bus* bus, REAL P, int t);
void BUS_inject_Q(Bus* bus, REAL Q, int t);
BOOL BUS_is_equal(Bus* bus, Bus* other);
BOOL BUS_is_regulated_by_gen(Bus* bus);
BOOL BUS_is_regulated_by_tran(Bus* bus);
BOOL BUS_is_regulated_by_shunt(Bus* bus);
BOOL BUS_is_regulated_by_vsc_conv(Bus* bus);
BOOL BUS_is_regulated_by_facts(Bus* bus);
BOOL BUS_is_v_set_regulated(Bus* bus);
BOOL BUS_is_slack(Bus* bus);
BOOL BUS_is_star(Bus* bus);
BOOL BUS_is_redundant(Bus* bus);

void BUS_equiv_add_to_net(Bus* bus, Net* net);
void BUS_equiv_add(Bus* bus, Bus* other_bus);
BOOL BUS_equiv_has(Bus* bus, Bus* other_bus);
BOOL BUS_equiv_has_slack(Bus* bus);
void BUS_equiv_make(Bus* bus1, Bus* bus2);
void BUS_equiv_del(Bus* bus);
int BUS_equiv_len(Bus* bus);
void BUS_equiv_show(Bus* bus);

Bus* BUS_list_add(Bus* bus_list, Bus* bus);
Bus* BUS_list_add_sorting(Bus* bus_list, Bus* bus, int sort_by, int t);
int BUS_list_len(Bus* bus_list);
void BUS_list_del(Bus* bus_list);

Bus* BUS_new(int num_periods);
void BUS_propagate_data_in_time(Bus* bus, int start, int end);
void BUS_set_area(Bus* bus, int area);
void BUS_set_zone(Bus* bus, int zone);
void BUS_set_network(Bus* bus, void* net);
void BUS_set_next(Bus* bus, Bus* next_bus);
void BUS_set_number(Bus* bus, int number);
void BUS_set_name(Bus* bus, char* name);
void BUS_set_alt_number(Bus* bus, int number);
void BUS_set_alt_name(Bus* bus, char* name);
void BUS_set_price(Bus* bus, REAL price, int t);
void BUS_set_v_base(Bus* bus, REAL v_base);
void BUS_set_v_mag(Bus* bus, REAL v_mag, int t);
void BUS_set_v_ang(Bus* bus, REAL v_ang, int t);
void BUS_set_v_set(Bus* bus, REAL v_set, int t);
void BUS_set_v_max_reg(Bus* bus, REAL v_max_reg);
void BUS_set_v_min_reg(Bus* bus, REAL v_min_reg);
void BUS_set_v_max_norm(Bus* bus, REAL v_max_norm);
void BUS_set_v_min_norm(Bus* bus, REAL v_min_norm);
void BUS_set_v_max_emer(Bus* bus, REAL v_max_emer);
void BUS_set_v_min_emer(Bus* bus, REAL v_min_emer);
void BUS_set_slack_flag(Bus* bus, BOOL slack);
void BUS_set_star_flag(Bus* bus, BOOL star);
void BUS_set_index(Bus* bus, int index);
int BUS_set_flags(void* bus, char flag_type, unsigned char mask, int index);
void BUS_set_var_values(Bus* bus, Vec* values);
void BUS_set_sens_P_balance(Bus* bus, REAL value, int t);
void BUS_set_sens_Q_balance(Bus* bus, REAL value, int t);
void BUS_set_sens_v_mag_u_bound(Bus* bus, REAL value, int t);
void BUS_set_sens_v_mag_l_bound(Bus* bus, REAL value, int t);
void BUS_set_sens_v_ang_u_bound(Bus* bus, REAL value, int t);
void BUS_set_sens_v_ang_l_bound(Bus* bus, REAL value, int t);
void BUS_set_sens_v_set_reg(Bus* bus, REAL value, int t);
void BUS_set_sens_v_reg_by_tran(Bus* bus, REAL value, int t);
void BUS_set_sens_v_reg_by_shunt(Bus* bus, REAL value, int t);
void BUS_set_P_mis(Bus* bus, REAL mis, int t);
void BUS_set_Q_mis(Bus* bus, REAL mis, int t);
void BUS_show(Bus* bus, int t);

int BUS_get_dPdw_index(Bus* bus, int t);
int BUS_get_dQdw_index(Bus* bus, int t);
int BUS_get_dPdv_index(Bus* bus, int t);
int BUS_get_dQdv_index(Bus* bus, int t);
int BUS_get_dwdw_index(Bus* bus, int t);
int BUS_get_dwdv_index(Bus* bus, int t);
int BUS_get_dvdv_index(Bus* bus, int t);

void BUS_set_dPdw_index(Bus* bus, int idx, int t);
void BUS_set_dQdw_index(Bus* bus, int idx, int t);
void BUS_set_dPdv_index(Bus* bus, int idx, int t);
void BUS_set_dQdv_index(Bus* bus, int idx, int t);
void BUS_set_dwdw_index(Bus* bus, int idx, int t);
void BUS_set_dwdv_index(Bus* bus, int idx, int t);
void BUS_set_dvdv_index(Bus* bus, int idx, int t);

int BUS_get_oindex(Bus* bus);
void BUS_set_oindex(Bus* bus, int oindex);

#endif
