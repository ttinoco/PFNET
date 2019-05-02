/** @file bus_dc.h
 *  @brief This file lists the constants and routines associated with the BusDC data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BUSDC_HEADER__
#define __BUSDC_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"
#include "uthash.h"

// Infinity
#define BUSDC_INF_V 1e8 /**< @brief Infinite voltage */

// Variables
/** \defgroup bus_dc_vars BusDC Variable Masks 
 *  @{
 */
#define BUSDC_VAR_V 0x01     /**< @brief Variable: DC bus voltage */
/** @} */

// Properties
/** \defgroup bus_dc_props BusDC Property Masks
 *  @{
 */
#define BUSDC_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup bus_dc_const BusDC Constants
 *  @{
 */
#define BUSDC_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define BUSDC_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define BUSDC_NUM_JSON_FIELDS 12   /**< @brief Constant: max number of json fields */
/** @} */

// BusDC
typedef struct BusDC BusDC;

// Others
typedef struct ConvCSC ConvCSC;
typedef struct ConvVSC ConvVSC;
typedef struct BranchDC BranchDC;

void BUSDC_add_branch_k(BusDC* bus, BranchDC* branch);
void BUSDC_add_branch_m(BusDC* bus, BranchDC* branch);
void BUSDC_add_vsc_conv(BusDC* bus, ConvVSC* conv);
void BUSDC_add_csc_conv(BusDC* bus, ConvCSC* conv);
void BUSDC_array_del(BusDC* bus_array, int size);
void* BUSDC_array_get(void* bus_array, int index);
void BUSDC_array_get_max_mismatches(BusDC* bus_array, int size, REAL* P, int t);
BusDC* BUSDC_array_new(int size, int num_periods);

void BUSDC_clear_flags(BusDC* bus, char flag_mask);
void BUSDC_clear_mismatches(BusDC* bus);
void BUSDC_clear_sensitivities(BusDC* bus);

void BUSDC_copy_from_dc_bus(BusDC* bus, BusDC* other);

void BUSDC_del_all_connections(BusDC* bus);
void BUSDC_del_branch_k(BusDC* bus, BranchDC* branch);
void BUSDC_del_branch_m(BusDC* bus, BranchDC* branch);
void BUSDC_del_vsc_conv(BusDC* bus, ConvVSC* conv);
void BUSDC_del_csc_conv(BusDC* bus, ConvCSC* conv);

BranchDC* BUSDC_get_branch_k(BusDC* bus);
BranchDC* BUSDC_get_branch_m(BusDC* bus);
ConvVSC* BUSDC_get_vsc_conv(BusDC* bus);
ConvCSC* BUSDC_get_csc_conv(BusDC* bus);
char BUSDC_get_flags_vars(BusDC* bus);
char BUSDC_get_flags_fixed(BusDC* bus);
char BUSDC_get_flags_bounded(BusDC* bus);
char BUSDC_get_flags_sparse(BusDC* bus);
int BUSDC_get_index(BusDC* bus);
int BUSDC_get_index_t(BusDC* bus, int t);
int BUSDC_get_index_v(BusDC* bus, int t);
char* BUSDC_get_json_string(BusDC* bus, char* output);
int BUSDC_get_num_csc_convs(BusDC* bus);
int BUSDC_get_num_vsc_convs(BusDC* bus);
int BUSDC_get_num_periods(BusDC* bus);
int BUSDC_get_num_vars(void* bus, unsigned char var, int t_start, int t_end);
int BUSDC_get_number(BusDC* bus);
char* BUSDC_get_name(BusDC* bus);
char BUSDC_get_obj_type(void* bus);
REAL BUSDC_get_P_mis(BusDC* bus, int t);
REAL BUSDC_get_v(BusDC* bus, int t);
REAL BUSDC_get_v_base(BusDC* bus);
Vec* BUSDC_get_var_indices(void* bus, unsigned char var, int t_start, int t_end);
char* BUSDC_get_var_info_string(BusDC* bus, int index);
void BUSDC_get_var_values(BusDC* bus, Vec* values, int code);

BOOL BUSDC_has_flags(void* vbus, char flag_type, unsigned char mask);
BOOL BUSDC_has_properties(void* bus, char prop);
BusDC* BUSDC_hash_number_add(BusDC* bus_hash, BusDC* bus);
void BUSDC_hash_number_del(BusDC* bus_hash);
BusDC* BUSDC_hash_number_find(BusDC* bus_hash, int number);
int BUSDC_hash_number_len(BusDC* bus_hash);
BusDC* BUSDC_hash_name_add(BusDC* bus_hash, BusDC* bus);
void BUSDC_hash_name_del(BusDC* bus_hash);
BusDC* BUSDC_hash_name_find(BusDC* bus_hash, char* name);
int BUSDC_hash_name_len(BusDC* bus_hash);

void BUSDC_init(BusDC* bus, int num_periods);
BOOL BUSDC_is_in_service(BusDC* bus);
BOOL BUSDC_is_equal(BusDC* bus, BusDC* other);

BusDC* BUSDC_new(int num_periods);
void BUSDC_propagate_data_in_time(BusDC* bus, int start, int end);

int BUSDC_set_flags(void* bus, char flag_type, unsigned char mask, int index);
void BUSDC_set_in_service(BusDC* bus, BOOL in_service);
void BUSDC_set_index(BusDC* bus, int index);
void BUSDC_set_number(BusDC* bus, int number);
void BUSDC_set_name(BusDC* bus, char* name);
void BUSDC_set_network(BusDC* bus, void* net);
void BUSDC_set_v(BusDC* bus, REAL v, int t);
void BUSDC_set_v_base(BusDC* bus, REAL v_base);
void BUSDC_set_var_values(BusDC* bus, Vec* values);

#endif
