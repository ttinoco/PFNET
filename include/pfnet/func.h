/** @file func.h
 *  @brief This file lists the constants and routines associated with the Func data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_HEADER__
#define __FUNC_HEADER__

#include "net.h"
#include "types.h"
#include "list.h"
#include "vector.h"
#include "matrix.h"

// Buffer
#define FUNC_BUFFER_SIZE 1024 /**< @brief Default function buffer size for strings */

// Function types
/** \defgroup func_types Function Types
 *  @{
 */
#define FUNC_TYPE_UNKNOWN -1    /**< @brief Function type: unknown */
#define FUNC_TYPE_REG_VMAG 0    /**< @brief Function type: voltage magnitude regulation. */
#define FUNC_TYPE_REG_VANG 1    /**< @brief Function type: vooltage angle regulation. */
#define FUNC_TYPE_REG_PQ 2      /**< @brief Function type: generator power regulation. */
#define FUNC_TYPE_REG_RATIO 3   /**< @brief Function type: transformer tap ratio regularization. */
#define FUNC_TYPE_REG_PHASE 4   /**< @brief Function type: transformer phase shift regularization. */
#define FUNC_TYPE_REG_SUSC 5    /**< @brief Function type: shunt susceptance regularization. */
#define FUNC_TYPE_GEN_COST 6    /**< @brief Function type: power generation cost. */
#define FUNC_TYPE_SP_CONTROLS 7 /**< @brief Function type: sparse controls. */
#define FUNC_TYPE_SLIM_VMAG 8   /**< @brief Function type: soft limits of bus voltage magnitudes. */
/** @} */

// Function
typedef struct Func Func;

// Function prototypes
void FUNC_clear_bus_counted(Func* f);
void FUNC_del(Func* f);
int FUNC_get_type(Func* f);
REAL FUNC_get_weight(Func* f);
REAL FUNC_get_phi(Func* f);
REAL* FUNC_get_phi_ptr(Func* f);
Vec* FUNC_get_gphi(Func* f);
Mat* FUNC_get_Hphi(Func* f);
int FUNC_get_Hcounter(Func* f);
int* FUNC_get_Hcounter_ptr(Func* f);
char* FUNC_get_bus_counted(Func* f);
int FUNC_get_bus_counted_size(Func* f);
Func* FUNC_get_next(Func* f);
int FUNC_get_branch_counter(Func* f);
void FUNC_inc_branch_counter(Func* f);
Func* FUNC_list_add(Func* flist, Func* nf);
int FUNC_list_len(Func* flist);
void FUNC_list_del(Func* flist);
void FUNC_list_count_branch(Func* f,Branch* br);
void FUNC_list_allocate(Func* f);
void FUNC_list_clear(Func* f);
void FUNC_list_analyze_branch(Func* f, Branch* br);
void FUNC_list_eval_branch(Func* f, Branch* br, Vec* var_values);
Func* FUNC_new(int type, REAL weight, Net* net);
void FUNC_set_phi(Func* f, REAL phi);
void FUNC_set_gphi(Func* f, Vec* gphi);
void FUNC_set_Hphi(Func* f, Mat* Hphi);
void FUNC_set_Hcounter(Func* f, int counter);
void FUNC_set_bus_counted(Func* f, char* counted, int size);
void FUNC_set_branch_counter(Func* f, int counter);
void FUNC_count(Func* f);
void FUNC_count_branch(Func* f, Branch* b);
void FUNC_allocate(Func* f);
void FUNC_clear(Func* f);
void FUNC_analyze(Func* f);
void FUNC_analyze_branch(Func* f, Branch* b);
void FUNC_eval(Func* f, Vec* var_values);
void FUNC_eval_branch(Func* f, Branch *b, Vec* var_values);
BOOL FUNC_is_safe_to_count(Func* f);
BOOL FUNC_is_safe_to_analyze(Func* f);
BOOL FUNC_is_safe_to_eval(Func* f, Vec* values);
BOOL FUNC_has_error(Func* f);
void FUNC_clear_error(Func * f);
char* FUNC_get_error_string(Func* f);
void FUNC_update_network(Func* f);
Net* FUNC_get_network(Func* f);

#endif
