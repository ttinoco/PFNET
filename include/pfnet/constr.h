/** @file constr.h
 *  @brief This file lists the constants and routines associated with the Constr data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_HEADER__
#define __CONSTR_HEADER__

#include "net.h"
#include "types.h"
#include "list.h"
#include "vector.h"
#include "matrix.h"

// Buffer
#define CONSTR_BUFFER_SIZE 1024 /**< @brief Default constraint buffer size for strings */

// Constraint types
/** \defgroup constr_types Constraint Types
 *  @{
 */
#define CONSTR_TYPE_UNKNOWN -1     /**< @brief Constraint type: Unknown. */
#define CONSTR_TYPE_PF 0           /**< @brief Constraint type: Power flow equations. */
#define CONSTR_TYPE_DCPF 1         /**< @brief Constraint type: DC power flow equations. */
#define CONSTR_TYPE_LINPF 2        /**< @brief Constraint type: Linearized power flow equations. */
#define CONSTR_TYPE_FIX 3          /**< @brief Constraint type: Variable fixing. */
#define CONSTR_TYPE_BOUND 4        /**< @brief Constraint type: Variable bounds as nonlinear equality constraints. */ 
#define CONSTR_TYPE_PAR_GEN_P 5    /**< @brief Constraint type: Generator participation (active power). */
#define CONSTR_TYPE_PAR_GEN_Q 6    /**< @brief Constraint type: Generator participation (reactive power). */
#define CONSTR_TYPE_REG_GEN 7      /**< @brief Constraint type: Voltage regualtion by generators. */
#define CONSTR_TYPE_REG_TRAN 8     /**< @brief Constraint type: Voltage regulation by transformers. */
#define CONSTR_TYPE_REG_SHUNT 9    /**< @brief Constraint type: Voltage regulation by shunt devices. */
#define CONSTR_TYPE_DC_FLOW_LIM 10 /**< @brief Constraint type: DC branch flow limits. */
#define CONSTR_TYPE_LBOUND 11      /**< @brief Constraint type: Variable bounds as linear inequality constraints. */
#define CONSTR_TYPE_GEN_RAMP 12    /**< @brief Constraint type: Generator active power ramping constraints. */
/** @} */

// Constraint types strings
#define CONSTR_TYPE_UNKNOWN_STR "UNKNOWN"
#define CONSTR_TYPE_PF_STR "PF"
#define CONSTR_TYPE_DCPF_STR "DCPF"
#define CONSTR_TYPE_LINPF_STR "LINPF"
#define CONSTR_TYPE_FIX_STR "FIX"
#define CONSTR_TYPE_BOUND_STR "BOUND"
#define CONSTR_TYPE_PAR_GEN_P_STR "PAR_GEN_P"
#define CONSTR_TYPE_PAR_GEN_Q_STR "PAR_GEN_Q"
#define CONSTR_TYPE_REG_GEN_STR "REG_GEN"
#define CONSTR_TYPE_REG_TRAN_STR "REG_TRAN"
#define CONSTR_TYPE_REG_SHUNT_STR "REG_SHUNT"
#define CONSTR_TYPE_DC_FLOW_LIM_STR "DC_FLOW_LIM"
#define CONSTR_TYPE_LBOUND_STR "LBOUND"
#define CONSTR_TYPE_GEN_RAMP_STR "GEN_RAMP"

// Constraint
typedef struct Constr Constr;

// Function prototypes
void CONSTR_clear_Hcounter(Constr* c);
void CONSTR_clear_bus_counted(Constr* c);
void CONSTR_combine_H(Constr* c, Vec* coeff, BOOL ensure_psd);
void CONSTR_del(Constr* constr);
void CONSTR_del_matvec(Constr* constr);
int CONSTR_get_type(Constr* c);
char* CONSTR_get_type_str(Constr* c);
Vec* CONSTR_get_b(Constr* c);
Mat* CONSTR_get_A(Constr* c);
Vec* CONSTR_get_l(Constr* c);
Vec* CONSTR_get_u(Constr* c);
Mat* CONSTR_get_G(Constr* c);
Vec* CONSTR_get_f(Constr* c);
Mat* CONSTR_get_J(Constr* c);
Mat* CONSTR_get_H_array(Constr* c);
int CONSTR_get_H_array_size(Constr* c);
Mat* CONSTR_get_H_single(Constr* c, int i);
Mat* CONSTR_get_H_combined(Constr* c);
int CONSTR_get_Acounter(Constr* c);
int* CONSTR_get_Acounter_ptr(Constr* c);
int CONSTR_get_Gcounter(Constr* c);
int* CONSTR_get_Gcounter_ptr(Constr* c);
int CONSTR_get_Jcounter(Constr* c);
int* CONSTR_get_Jcounter_ptr(Constr* c);
int* CONSTR_get_Hcounter(Constr* c);
int CONSTR_get_Hcounter_size(Constr* c);
int CONSTR_get_Aconstr_index(Constr* c);
int* CONSTR_get_Aconstr_index_ptr(Constr* c);
int CONSTR_get_Gconstr_index(Constr* c);
int* CONSTR_get_Gconstr_index_ptr(Constr* c);
int CONSTR_get_Jconstr_index(Constr* c);
int* CONSTR_get_Jconstr_index_ptr(Constr* c);
char* CONSTR_get_bus_counted(Constr *c);
int CONSTR_get_bus_counted_size(Constr* c);
void* CONSTR_get_data(Constr* c);
Constr* CONSTR_get_next(Constr* c);
Constr* CONSTR_list_add(Constr* clist, Constr* nc);
int CONSTR_list_len(Constr* clist);
void CONSTR_list_del(Constr* clist);
void CONSTR_list_combine_H(Constr* clist, Vec* coeff, BOOL ensure_psd);
void CONSTR_list_count_step(Constr* clist, Branch* br, int t);
void CONSTR_list_allocate(Constr* clist);
void CONSTR_list_clear(Constr* clist);
void CONSTR_list_analyze_step(Constr* clist, Branch* br, int t);
void CONSTR_list_eval_step(Constr* clist, Branch* br, int t, Vec* var_values);
void CONSTR_list_store_sens_step(Constr* clist, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
Constr* CONSTR_new(int type, Net* net);
void CONSTR_set_b(Constr* c, Vec* b);
void CONSTR_set_A(Constr* c, Mat* A);
void CONSTR_set_l(Constr* c, Vec* l);
void CONSTR_set_u(Constr* c, Vec* u);
void CONSTR_set_G(Constr* c, Mat* G);
void CONSTR_set_f(Constr* c, Vec* f);
void CONSTR_set_J(Constr* c, Mat* J);
void CONSTR_set_H_array(Constr* c, Mat* H_array, int size);
void CONSTR_set_H_combined(Constr* c, Mat* H_combined);
void CONSTR_set_Acounter(Constr* c, int counter);
void CONSTR_set_Gcounter(Constr* c, int counter);
void CONSTR_set_Jcounter(Constr* c, int counter);
void CONSTR_set_Hcounter(Constr* c, int* counter, int size);
void CONSTR_set_Aconstr_index(Constr* c, int index);
void CONSTR_set_Gconstr_index(Constr* c, int index);
void CONSTR_set_Jconstr_index(Constr* c, int index);
void CONSTR_set_bus_counted(Constr* c, char* counted, int size);
void CONSTR_set_data(Constr* c, void* data);
void CONSTR_count(Constr* c);
void CONSTR_count_step(Constr* c, Branch* br, int t);
void CONSTR_allocate(Constr* c);
void CONSTR_clear(Constr* c);
void CONSTR_analyze(Constr* c);
void CONSTR_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_eval(Constr* c, Vec* var_values);
void CONSTR_eval_step(Constr* c, Branch* br, int t, Vec* var_values);
void CONSTR_store_sens(Constr* c, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
BOOL CONSTR_is_safe_to_count(Constr* c);
BOOL CONSTR_is_safe_to_analyze(Constr* c);
BOOL CONSTR_is_safe_to_eval(Constr* c, Vec* values);
BOOL CONSTR_has_error(Constr* c);
void CONSTR_set_error(Constr* c, char* string);
void CONSTR_clear_error(Constr* c);
char* CONSTR_get_error_string(Constr* c);
void CONSTR_update_network(Constr* c);
Net* CONSTR_get_network(Constr* c);

#endif
