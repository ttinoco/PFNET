/** @file conv_vsc.h
 *  @brief This file lists the constants and routines associated with the ConvVSC data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONVVSC_HEADER__
#define __CONVVSC_HEADER__

#include "stdio.h"
#include "types.h"
#include "list.h"
#include "vector.h"

// VSC converter AC modes
#define CONVVSC_MODE_AC_NC 0 /**< @brief Mode: No control */
#define CONVVSC_MODE_AC_CV 1 /**< @brief Mode: Constant AC voltage */
#define CONVVSC_MODE_AC_CF 2 /**< @brief Mode: Constant AC power factor */

// VSC converter DC modes
#define CONVVSC_MODE_DC_NC 0 /**< @brief Mode: No control */
#define CONVVSC_MODE_DC_CV 1 /**< @brief Mode: Constant DC voltage */
#define CONVVSC_MODE_DC_CP 2 /**< @brief Mode: Constant DC power */

// Infinity
#define CONVVSC_INF_P 1e8     /**< @brief Infinite active power */
#define CONVVSC_INF_Q 1e8     /**< @brief Infinite reactive power */
#define CONVVSC_INF_PDC 1e8   /**< @brief Infinite DC power */

// Variables
/** \defgroup conv_vsc_vars ConvVSC Variable Masks 
 *  @{
 */
#define CONVVSC_VAR_P 0x01     /**< @brief Variable: converter active power injection into AC bus */
#define CONVVSC_VAR_Q 0x02     /**< @brief Variable: converter reactive power injection into AC bus */
#define CONVVSC_VAR_PDC 0x04   /**< @brief Variable: converter power injection into DC bus*/
/** @} */

// Properties
/** \defgroup conv_vsc_props ConvVSC Property Masks
 *  @{
 */
#define CONVVSC_PROP_ANY 0x00     /**< @brief Property: any */
/** @} */

// Constants
/** \defgroup conv_vsc_const ConvVSC Constants
 *  @{
 */
#define CONVVSC_BUFFER_SIZE 100      /**< @brief Constant: buffer size for strings */
#define CONVVSC_JSON_BUFFER_SIZE 200 /**< @brief Constant: buffer size for json strings */
#define CONVVSC_NUM_JSON_FIELDS 25   /**< @brief Constant: max number of json fields */
#define CONVVSC_MIN_TARGET_PF 1e-4   /**< @brief Minimum target power factor */
/** @} */

// ConvVSC
typedef struct ConvVSC ConvVSC;

// Others
typedef struct Bus Bus;
typedef struct BusDC BusDC;

void CONVVSC_array_del(ConvVSC* conv_array, int size);
void* CONVVSC_array_get(void* conv_array, int index);
ConvVSC* CONVVSC_array_new(int size, int num_periods);

void CONVVSC_clear_flags(ConvVSC* conv, char flag_mask);
void CONVVSC_clear_sensitivities(ConvVSC* conv);

void CONVVSC_copy_from_conv(ConvVSC* conv, ConvVSC* other);

Bus* CONVVSC_get_ac_bus(ConvVSC* conv);
BusDC* CONVVSC_get_dc_bus(ConvVSC* conv);
Bus* CONVVSC_get_reg_bus(ConvVSC* conv);
char CONVVSC_get_flags_vars(ConvVSC* bus);
char CONVVSC_get_flags_fixed(ConvVSC* bus);
char CONVVSC_get_flags_bounded(ConvVSC* bus);
char CONVVSC_get_flags_sparse(ConvVSC* bus);

int CONVVSC_get_index(ConvVSC* conv);
int CONVVSC_get_index_P(ConvVSC* conv, int t);
int CONVVSC_get_index_Q(ConvVSC* conv, int t);
int CONVVSC_get_index_P_dc(ConvVSC* conv, int t);
int CONVVSC_get_index_i_dc(ConvVSC* conv, int t);

int* CONVVSC_get_index_P_array(ConvVSC* conv);
int* CONVVSC_get_index_Q_array(ConvVSC* conv);
int* CONVVSC_get_index_P_dc_array(ConvVSC* conv);
int* CONVVSC_get_index_i_dc_array(ConvVSC* conv);
  
char* CONVVSC_get_json_string(ConvVSC* conv, char* output);
char CONVVSC_get_mode_ac(ConvVSC* conv);
char CONVVSC_get_mode_dc(ConvVSC* conv);
char* CONVVSC_get_name(ConvVSC* conv);
ConvVSC* CONVVSC_get_next_ac(ConvVSC* conv);
ConvVSC* CONVVSC_get_next_dc(ConvVSC* conv);
ConvVSC* CONVVSC_get_reg_next(ConvVSC* conv);
int CONVVSC_get_num_periods(ConvVSC* conv);
int CONVVSC_get_num_vars(void* conv, unsigned char var, int t_start, int t_end);
char CONVVSC_get_obj_type(void* conv);
REAL CONVVSC_get_P(ConvVSC* conv, int t);
REAL CONVVSC_get_P_dc(ConvVSC* conv, int t);
REAL CONVVSC_get_i_dc(ConvVSC* conv, int t);
REAL CONVVSC_get_Q(ConvVSC* conv, int t);
REAL CONVVSC_get_P_max(ConvVSC* conv);
REAL CONVVSC_get_P_min(ConvVSC* conv);
REAL CONVVSC_get_Q_max(ConvVSC* conv);
REAL CONVVSC_get_Q_min(ConvVSC* conv);
REAL CONVVSC_get_Q_par(ConvVSC* conv);
REAL CONVVSC_get_loss_coeff_A(ConvVSC* conv);
REAL CONVVSC_get_loss_coeff_B(ConvVSC* conv);
REAL CONVVSC_get_target_power_factor(ConvVSC* conv);
REAL CONVVSC_get_v_dc_set(ConvVSC* conv, int t);
REAL CONVVSC_get_P_dc_set(ConvVSC* conv, int t);

REAL* CONVVSC_get_P_array(ConvVSC* conv);
REAL* CONVVSC_get_P_dc_array(ConvVSC* conv);
REAL* CONVVSC_get_Q_array(ConvVSC* conv);
REAL* CONVVSC_get_v_dc_set_array(ConvVSC* conv);
REAL* CONVVSC_get_P_dc_set_array(ConvVSC* conv);

Vec* CONVVSC_get_var_indices(void* conv, unsigned char var, int t_start, int t_end);
char* CONVVSC_get_var_info_string(ConvVSC* conv, int index);
void CONVVSC_get_var_values(ConvVSC* conv, Vec* values, int code);

BOOL CONVVSC_has_flags(void* vconv, char flag_type, unsigned char mask);
BOOL CONVVSC_has_properties(void* conv, char prop);

void CONVVSC_init(ConvVSC* conv, int num_periods);
BOOL CONVVSC_is_in_service(ConvVSC* conv);
BOOL CONVVSC_is_equal(ConvVSC* conv, ConvVSC* other);
BOOL CONVVSC_is_in_f_ac_mode(ConvVSC* conv);
BOOL CONVVSC_is_in_v_ac_mode(ConvVSC* conv);
BOOL CONVVSC_is_in_v_dc_mode(ConvVSC* conv);
BOOL CONVVSC_is_in_P_dc_mode(ConvVSC* conv);

ConvVSC* CONVVSC_list_ac_add(ConvVSC* conv_list, ConvVSC* conv);
ConvVSC* CONVVSC_list_ac_del(ConvVSC* conv_list, ConvVSC* conv);
int CONVVSC_list_ac_len(ConvVSC* conv_list);
ConvVSC* CONVVSC_list_dc_add(ConvVSC* conv_list, ConvVSC* conv);
ConvVSC* CONVVSC_list_dc_del(ConvVSC* conv_list, ConvVSC* conv);
int CONVVSC_list_dc_len(ConvVSC* conv_list);
ConvVSC* CONVVSC_list_reg_add(ConvVSC* conv_list, ConvVSC* conv);
ConvVSC* CONVVSC_list_reg_del(ConvVSC* conv_list, ConvVSC* conv);
int CONVVSC_list_reg_len(ConvVSC* conv_list);

ConvVSC* CONVVSC_new(int num_periods);
void CONVVSC_propagate_data_in_time(ConvVSC* conv, int start, int end);

void CONVVSC_set_in_service(ConvVSC* conv, BOOL in_service);
void CONVVSC_set_ac_bus(ConvVSC* conv, Bus* bus);
void CONVVSC_set_dc_bus(ConvVSC* conv, BusDC* bus);
void CONVVSC_set_reg_bus(ConvVSC* conv, Bus* reg_bus);
int CONVVSC_set_flags(void* conv, char flag_type, unsigned char mask, int index);
void CONVVSC_set_index(ConvVSC* conv, int index);
void CONVVSC_set_mode_ac(ConvVSC* conv, char mode);
void CONVVSC_set_mode_dc(ConvVSC* conv, char mode);
void CONVVSC_set_name(ConvVSC* conv, char* name);
void CONVVSC_set_P(ConvVSC* conv, REAL P, int t);
void CONVVSC_set_P_dc(ConvVSC* conv, REAL P, int t);
void CONVVSC_set_Q(ConvVSC* conv, REAL Q, int t);
void CONVVSC_set_P_max(ConvVSC* conv, REAL P_max);
void CONVVSC_set_P_min(ConvVSC* conv, REAL P_min);
void CONVVSC_set_Q_max(ConvVSC* conv, REAL Q_max);
void CONVVSC_set_Q_min(ConvVSC* conv, REAL Q_min);
void CONVVSC_set_Q_par(ConvVSC* conv, REAL Q_par);
void CONVVSC_set_loss_coeff_A(ConvVSC* conv, REAL A);
void CONVVSC_set_loss_coeff_B(ConvVSC* conv, REAL B);

void CONVVSC_set_var_values(ConvVSC* conv, Vec* values);

void CONVVSC_set_v_dc_set(ConvVSC* conv, REAL v, int t);
void CONVVSC_set_P_dc_set(ConvVSC* conv, REAL P, int t);
void CONVVSC_set_target_power_factor(ConvVSC* conv, REAL pf);
  
#endif
