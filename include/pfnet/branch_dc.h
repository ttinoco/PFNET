/** @file branch_dc.h
 *  @brief This file lists the constants and routines associated with the BranchDC data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __BRANCHDC_HEADER__
#define __BRANCHDC_HEADER__

#include "types.h"
#include "list.h"
#include "vector.h"

// Properties
/** \defgroup branch_dc_props BranchDC Property Masks
 *  @{
 */
#define BRANCHDC_PROP_ANY 0x00           /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup branch_dc_const BranchDC Constants
 *  @{
 */
#define BRANCHDC_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define BRANCHDC_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define BRANCHDC_NUM_JSON_FIELDS 10   /**< @brief Constant: max number of json fields */
/** @} */

// BranchDC
typedef struct BranchDC BranchDC;

// Other
typedef struct BusDC BusDC;

// Prototypes
void BRANCHDC_array_del(BranchDC* br_array, int size);
void* BRANCHDC_array_get(void* br, int index);
BranchDC* BRANCHDC_array_new(int size, int num_periods);

void BRANCHDC_clear_sensitivities(BranchDC* br);
void BRANCHDC_clear_flags(BranchDC* br, char flag_type);

void BRANCHDC_copy_from_dc_branch(BranchDC* br, BranchDC* other);

BusDC* BRANCHDC_get_bus_k(BranchDC* br);
BusDC* BRANCHDC_get_bus_m(BranchDC* br);
char BRANCHDC_get_flags_vars(BranchDC* br);
char BRANCHDC_get_flags_fixed(BranchDC* br);
char BRANCHDC_get_flags_bounded(BranchDC* br);
char BRANCHDC_get_flags_sparse(BranchDC* br);
REAL BRANCHDC_get_i_km(BranchDC* br, Vec* var_values, int t);
REAL BRANCHDC_get_i_mk(BranchDC* br, Vec* var_values, int t);
int BRANCHDC_get_index(BranchDC* br);
char* BRANCHDC_get_json_string(BranchDC* br, char* output);
char* BRANCHDC_get_name(BranchDC* br);
BranchDC* BRANCHDC_get_next_k(BranchDC* br);
BranchDC* BRANCHDC_get_next_m(BranchDC* br);
int BRANCHDC_get_num_periods(BranchDC* br);
int BRANCHDC_get_num_vars(void* br, unsigned char var, int t_start, int t_end);
char BRANCHDC_get_obj_type(void* br);
REAL BRANCHDC_get_P_km(BranchDC* br, Vec* var_values, int t);
REAL BRANCHDC_get_P_mk(BranchDC* br, Vec* var_values, int t);
REAL BRANCHDC_get_r(BranchDC* br);
Vec* BRANCHDC_get_var_indices(void* br, unsigned char var, int t_start, int t_end);
char* BRANCHDC_get_var_info_string(BranchDC* br, int index);
void BRANCHDC_get_var_values(BranchDC* br, Vec* values, int code);

BOOL BRANCHDC_has_flags(void* br, char flag_type, unsigned char mask);
BOOL BRANCHDC_has_properties(void* br, char prop);

void BRANCHDC_init(BranchDC* br, int num_periods);
BOOL BRANCHDC_is_in_service(void* br);
BOOL BRANCHDC_is_equal(BranchDC* br, BranchDC* other);

BranchDC* BRANCHDC_list_k_add(BranchDC* k_br_list, BranchDC* br);
BranchDC* BRANCHDC_list_k_del(BranchDC* k_br_list, BranchDC* br);
int BRANCHDC_list_k_len(BranchDC* k_br_list);
BranchDC* BRANCHDC_list_m_add(BranchDC* m_br_list, BranchDC* br);
BranchDC* BRANCHDC_list_m_del(BranchDC* m_br_list, BranchDC* br);
int BRANCHDC_list_m_len(BranchDC* m_br_list);

BranchDC* BRANCHDC_new(int num_periods);
void BRANCHDC_propagate_data_in_time(BranchDC* br, int start, int end);

void BRANCHDC_set_network(BranchDC* br, void* net);
void BRANCHDC_set_in_service(BranchDC* br, BOOL in_service);
void BRANCHDC_set_bus_k(BranchDC* br, BusDC* bus_k);
void BRANCHDC_set_bus_m(BranchDC* br, BusDC* bus_m);
int BRANCHDC_set_flags(void* vbr, char flag_type, unsigned char mask, int index);
void BRANCHDC_set_index(BranchDC* br, int index);
void BRANCHDC_set_name(BranchDC* br, char* name);
void BRANCHDC_set_r(BranchDC* br, REAL r);
void BRANCHDC_set_var_values(BranchDC* br, Vec* values);

#endif
