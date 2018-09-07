/** @file conv_csc.h
 *  @brief This file lists the constants and routines associated with the ConvCSC data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONVCSC_HEADER__
#define __CONVCSC_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// CSC converter types
#define CONVCSC_TYPE_REC 0     /**< @brief Type: Rectifier */
#define CONVCSC_TYPE_INV 1     /**< @brief Type: Inverter */
#define CONVCSC_TYPE_UNKNOWN 2 /**< @brief Type: Unknown */

// CSC converter modes
#define CONVCSC_MODE_DC_NC 0 /**< @brief Mode: No control */
#define CONVCSC_MODE_DC_CP 1 /**< @brief Mode: Constant power */
#define CONVCSC_MODE_DC_CC 2 /**< @brief Mode: Constant current */
#define CONVCSC_MODE_DC_CV 3 /**< @brief Mode: Constant voltage */

// Infinity
#define CONVCSC_INF_P 1e8     /**< @brief Infinite active power */
#define CONVCSC_INF_Q 1e8     /**< @brief Infinite reactive power */
#define CONVCSC_INF_PDC 1e8   /**< @brief Infinite DC power */
#define CONVCSC_INF_RATIO 1e8 /**< @brief Infinite tap ratio */
#define CONVCSC_INF_ANGLE 1e8 /**< @brief Infinite phase shift */

// Variables
/** \defgroup conv_csc_vars ConvCSC Variable Masks 
 *  @{
 */
#define CONVCSC_VAR_P 0x01     /**< @brief Variable: converter active power injection into AC bus */
#define CONVCSC_VAR_Q 0x02     /**< @brief Variable: converter reactive power injection into AC bus */
#define CONVCSC_VAR_PDC 0x04   /**< @brief Variable: converter power injection into DC bus*/
#define CONVCSC_VAR_RATIO 0x08 /**< @brief Variable: converter transformer taps ratio */
#define CONVCSC_VAR_ANGLE 0x10 /**< @brief Variable: converter ignition or extinction angle */
/** @} */

// Properties
/** \defgroup conv_csc_props ConvCSC Property Masks
 *  @{
 */
#define CONVCSC_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup conv_csc_const ConvCSC Constants
 *  @{
 */
#define CONVCSC_BUFFER_SIZE 100    /**< @brief Constant: buffer size for strings */
#define CONVCSC_NUM_JSON_FIELDS 25 /**< @brief Constant: max number of json fields */
/** @} */

// ConvCSC
typedef struct ConvCSC ConvCSC;

// Others
typedef struct Bus Bus;
typedef struct BusDC BusDC;

void CONVCSC_array_del(ConvCSC* conv_array, int size);
void* CONVCSC_array_get(void* conv_array, int index);
ConvCSC* CONVCSC_array_new(int size, int num_periods);

void CONVCSC_clear_flags(ConvCSC* conv, char flag_mask);
void CONVCSC_clear_sensitivities(ConvCSC* conv);

void CONVCSC_copy_from_conv(ConvCSC* conv, ConvCSC* other);

Bus* CONVCSC_get_ac_bus(ConvCSC* conv);
BusDC* CONVCSC_get_dc_bus(ConvCSC* conv);
REAL CONVCSC_get_angle(ConvCSC* conv, int t);
REAL CONVCSC_get_angle_max(ConvCSC* conv);
REAL CONVCSC_get_angle_min(ConvCSC* conv);
char CONVCSC_get_flags_vars(ConvCSC* bus);
char CONVCSC_get_flags_fixed(ConvCSC* bus);
char CONVCSC_get_flags_bounded(ConvCSC* bus);
char CONVCSC_get_flags_sparse(ConvCSC* bus);
REAL CONVCSC_get_i_dc_set(ConvCSC* conv, int t);
int CONVCSC_get_index(ConvCSC* conv);
int CONVCSC_get_index_P(ConvCSC* conv, int t);
int CONVCSC_get_index_Q(ConvCSC* conv, int t);
int CONVCSC_get_index_P_dc(ConvCSC* conv, int t);
int CONVCSC_get_index_i_dc(ConvCSC* conv, int t);
int CONVCSC_get_index_ratio(ConvCSC* conv, int t);
int CONVCSC_get_index_angle(ConvCSC* conv, int t);
char* CONVCSC_get_json_string(ConvCSC* conv, char* output);
char CONVCSC_get_mode_dc(ConvCSC* conv);
char* CONVCSC_get_name(ConvCSC* conv);
ConvCSC* CONVCSC_get_next_ac(ConvCSC* conv);
ConvCSC* CONVCSC_get_next_dc(ConvCSC* conv);
int CONVCSC_get_num_bridges(ConvCSC* conv);
int CONVCSC_get_num_periods(ConvCSC* conv);
int CONVCSC_get_num_vars(void* conv, unsigned char var, int t_start, int t_end);
char CONVCSC_get_obj_type(void* conv);
REAL CONVCSC_get_P(ConvCSC* conv, int t);
REAL CONVCSC_get_P_dc(ConvCSC* conv, int t);
REAL CONVCSC_get_i_dc(ConvCSC* conv, int t);
REAL CONVCSC_get_P_dc_set(ConvCSC* conv, int t);
REAL CONVCSC_get_Q(ConvCSC* conv, int t);
REAL CONVCSC_get_ratio(ConvCSC* conv, int t);
REAL CONVCSC_get_ratio_max(ConvCSC* conv);
REAL CONVCSC_get_ratio_min(ConvCSC* conv);
REAL CONVCSC_get_v_base_p(ConvCSC* conv);
REAL CONVCSC_get_v_base_s(ConvCSC* conv);
REAL CONVCSC_get_v_dc_set(ConvCSC* conv, int t);
Vec* CONVCSC_get_var_indices(void* conv, unsigned char var, int t_start, int t_end);
char* CONVCSC_get_var_info_string(ConvCSC* conv, int index);
void CONVCSC_get_var_values(ConvCSC* conv, Vec* values, int code);
REAL CONVCSC_get_x(ConvCSC* conv);
REAL CONVCSC_get_r(ConvCSC* conv);
REAL CONVCSC_get_x_cap(ConvCSC* conv);

BOOL CONVCSC_has_flags(void* vconv, char flag_type, unsigned char mask);
BOOL CONVCSC_has_properties(void* conv, char prop);

void CONVCSC_init(ConvCSC* conv, int num_periods);
BOOL CONVCSC_is_equal(ConvCSC* conv, ConvCSC* other);
BOOL CONVCSC_is_inverter(ConvCSC* conv);
BOOL CONVCSC_is_rectifier(ConvCSC* conv);
BOOL CONVCSC_is_in_P_dc_mode(ConvCSC* conv);
BOOL CONVCSC_is_in_i_dc_mode(ConvCSC* conv);
BOOL CONVCSC_is_in_v_dc_mode(ConvCSC* conv);

ConvCSC* CONVCSC_list_ac_add(ConvCSC* conv_list, ConvCSC* conv);
ConvCSC* CONVCSC_list_ac_del(ConvCSC* conv_list, ConvCSC* conv);
int CONVCSC_list_ac_len(ConvCSC* conv_list);
ConvCSC* CONVCSC_list_dc_add(ConvCSC* conv_list, ConvCSC* conv);
ConvCSC* CONVCSC_list_dc_del(ConvCSC* conv_list, ConvCSC* conv);
int CONVCSC_list_dc_len(ConvCSC* conv_list);

ConvCSC* CONVCSC_new(int num_periods);
void CONVCSC_propagate_data_in_time(ConvCSC* conv, int start, int end);

void CONVCSC_set_ac_bus(ConvCSC* conv, Bus* bus);
void CONVCSC_set_dc_bus(ConvCSC* conv, BusDC* bus);
void CONVCSC_set_angle(ConvCSC* conv, REAL angle, int t);
void CONVCSC_set_angle_max(ConvCSC* conv, REAL angle_max);
void CONVCSC_set_angle_min(ConvCSC* conv, REAL angle_min);
int CONVCSC_set_flags(void* conv, char flag_type, unsigned char mask, int index);
void CONVCSC_set_i_dc_set(ConvCSC* conv, REAL i, int t);
void CONVCSC_set_index(ConvCSC* conv, int index);
void CONVCSC_set_mode_dc(ConvCSC* conv, char mode);
void CONVCSC_set_name(ConvCSC* conv, char* name);
void CONVCSC_set_num_bridges(ConvCSC* conv, int num);
void CONVCSC_set_P(ConvCSC* conv, REAL P, int t);
void CONVCSC_set_P_dc(ConvCSC* conv, REAL P, int t);
void CONVCSC_set_P_dc_set(ConvCSC* conv, REAL P, int t);
void CONVCSC_set_Q(ConvCSC* conv, REAL Q, int t);
void CONVCSC_set_ratio(ConvCSC* conv, REAL ratio, int t);
void CONVCSC_set_ratio_max(ConvCSC* conv, REAL ratio_max);
void CONVCSC_set_ratio_min(ConvCSC* conv, REAL ratio_min);
void CONVCSC_set_type(ConvCSC* conv, char type);
void CONVCSC_set_v_base_p(ConvCSC* conv, REAL v_base_p);
void CONVCSC_set_v_base_s(ConvCSC* conv, REAL v_base_s);
void CONVCSC_set_v_dc_set(ConvCSC* conv, REAL v, int t);
void CONVCSC_set_var_values(ConvCSC* conv, Vec* values);
void CONVCSC_set_x_cap(ConvCSC* conv, REAL x_cap);
void CONVCSC_set_x(ConvCSC* conv, REAL x);
void CONVCSC_set_r(ConvCSC* conv, REAL r);

#endif
