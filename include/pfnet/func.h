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

// Function
typedef struct Func Func;

// Function prototypes
void FUNC_clear_bus_counted(Func* f);
void FUNC_del(Func* f);
void FUNC_del_matvec(Func* f);
char* FUNC_get_name(Func* f);
REAL FUNC_get_weight(Func* f);
REAL FUNC_get_phi(Func* f);
REAL* FUNC_get_phi_ptr(Func* f);
Vec* FUNC_get_gphi(Func* f);
Mat* FUNC_get_Hphi(Func* f);
int FUNC_get_Hphi_nnz(Func* f);
int* FUNC_get_Hphi_nnz_ptr(Func* f);
char* FUNC_get_bus_counted(Func* f);
int FUNC_get_bus_counted_size(Func* f);
Func* FUNC_get_next(Func* f);
void* FUNC_get_data(Func* f);
void FUNC_list_clear_error(Func * flist);
BOOL FUNC_list_has_error(Func* flist);
char* FUNC_list_get_error_string(Func* flist);
Func* FUNC_list_add(Func* flist, Func* nf);
int FUNC_list_len(Func* flist);
void FUNC_list_del(Func* flist);
void FUNC_list_count_step(Func* f, Branch* br, int t);
void FUNC_list_allocate(Func* f);
void FUNC_list_clear(Func* f);
void FUNC_list_analyze_step(Func* f, Branch* br, int t);
void FUNC_list_eval_step(Func* f, Branch* br, int t, Vec* var_values);
void FUNC_list_finalize_structure_of_Hessian(Func* flist);
void FUNC_finalize_structure_of_Hessian(Func* f);
Func* FUNC_new(REAL weight, Net* net);
void FUNC_set_parameter(Func* f, char* key, void* value);
void FUNC_set_name(Func* f, char* name);
void FUNC_set_phi(Func* f, REAL phi);
void FUNC_set_gphi(Func* f, Vec* gphi);
void FUNC_set_Hphi(Func* f, Mat* Hphi);
void FUNC_set_Hphi_nnz(Func* f, int nnz);
void FUNC_set_bus_counted(Func* f, char* counted, int size);
void FUNC_set_data(Func* f, void* data);
void FUNC_init(Func* f);
void FUNC_count(Func* f);
void FUNC_count_step(Func* f, Branch* br, int t);
void FUNC_allocate(Func* f);
void FUNC_clear(Func* f);
void FUNC_analyze(Func* f);
void FUNC_analyze_step(Func* f, Branch* br, int t);
void FUNC_eval(Func* f, Vec* var_values);
void FUNC_eval_step(Func* f, Branch* br, int t, Vec* var_values);
BOOL FUNC_is_safe_to_count(Func* f);
BOOL FUNC_is_safe_to_analyze(Func* f);
BOOL FUNC_is_safe_to_eval(Func* f, Vec* values);
void FUNC_clear_error(Func * f);
BOOL FUNC_has_error(Func* f);
char* FUNC_get_error_string(Func* f);
void FUNC_update_network(Func* f);
Net* FUNC_get_network(Func* f);
void FUNC_set_func_init(Func* f, void (*func)(Func* f));
void FUNC_set_func_count_step(Func* f, void (*func)(Func* f, Branch* br, int t));
void FUNC_set_func_allocate(Func* f, void (*func)(Func* f));
void FUNC_set_func_clear(Func* f, void (*func)(Func* f));
void FUNC_set_func_analyze_step(Func* f, void (*func)(Func* f, Branch* br, int t));
void FUNC_set_func_eval_step(Func* f, void (*func)(Func* f, Branch* br, int t, Vec* v));
void FUNC_set_func_free(Func* f, void (*func)(Func* f));
void FUNC_set_func_set_parameter(Func* f, void (*func)(Func* f, char* key, void* value));

#endif
