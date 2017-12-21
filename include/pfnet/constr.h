/** @file constr.h
 *  @brief This file lists the constants and routines associated with the Constr data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
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
#define CONSTR_BUFFER_SIZE 1024     /**< @brief Default constraint buffer size for general strings */
#define CONSTR_INFO_BUFFER_SIZE 100 /**< @brife Default buffer size for row info strings */

// Constraint
typedef struct Constr Constr;

// Function prototypes
void CONSTR_allocate_H_array(Constr* c, int size);
void CONSTR_allocate_H_combined(Constr* c);
void CONSTR_finalize_structure_of_Hessians(Constr* c);
void CONSTR_clear_H_nnz(Constr* c);
void CONSTR_clear_bus_counted(Constr* c);
void CONSTR_combine_H(Constr* c, Vec* coeff, BOOL ensure_psd);
void CONSTR_del(Constr* constr);
void CONSTR_del_matvec(Constr* constr);
char* CONSTR_get_name(Constr* c);
Vec* CONSTR_get_b(Constr* c);
Mat* CONSTR_get_A(Constr* c);
Vec* CONSTR_get_l(Constr* c);
Vec* CONSTR_get_u(Constr* c);
Vec* CONSTR_get_l_extra_vars(Constr* c);
Vec* CONSTR_get_u_extra_vars(Constr* c);
Vec* CONSTR_get_init_extra_vars(Constr* c);
Mat* CONSTR_get_G(Constr* c);
Vec* CONSTR_get_f(Constr* c);
Mat* CONSTR_get_J(Constr* c);
Mat* CONSTR_get_H_array(Constr* c);
int CONSTR_get_H_array_size(Constr* c);
Mat* CONSTR_get_H_single(Constr* c, int i);
Mat* CONSTR_get_H_combined(Constr* c);
int CONSTR_get_A_nnz(Constr* c);
int* CONSTR_get_A_nnz_ptr(Constr* c);
int CONSTR_get_G_nnz(Constr* c);
int* CONSTR_get_G_nnz_ptr(Constr* c);
int CONSTR_get_J_nnz(Constr* c);
int* CONSTR_get_J_nnz_ptr(Constr* c);
int* CONSTR_get_H_nnz(Constr* c);
int CONSTR_get_H_nnz_size(Constr* c);
int CONSTR_get_A_row(Constr* c);
int* CONSTR_get_A_row_ptr(Constr* c);
int CONSTR_get_G_row(Constr* c);
int* CONSTR_get_G_row_ptr(Constr* c);
int CONSTR_get_J_row(Constr* c);
int* CONSTR_get_J_row_ptr(Constr* c);
char* CONSTR_get_bus_counted(Constr *c);
int CONSTR_get_bus_counted_size(Constr* c);
void* CONSTR_get_data(Constr* c);
Constr* CONSTR_get_next(Constr* c);
Mat* CONSTR_get_var_projection(Constr* c);
Mat* CONSTR_get_extra_var_projection(Constr* c);
char* CONSTR_get_A_row_info_string(Constr* c, int index);
char* CONSTR_get_J_row_info_string(Constr* c, int index);
char* CONSTR_get_G_row_info_string(Constr* c, int index);
void CONSTR_list_finalize_structure_of_Hessians(Constr* clist);
void CONSTR_list_clear_error(Constr* clist);
BOOL CONSTR_list_has_error(Constr* clist);
char* CONSTR_list_get_error_string(Constr* clist);
Constr* CONSTR_list_add(Constr* clist, Constr* nc);
int CONSTR_list_len(Constr* clist);
void CONSTR_list_del(Constr* clist);
void CONSTR_list_combine_H(Constr* clist, Vec* coeff, BOOL ensure_psd);
void CONSTR_list_count_step(Constr* clist, Branch* br, int t);
void CONSTR_list_allocate(Constr* clist);
void CONSTR_list_clear(Constr* clist);
void CONSTR_list_analyze_step(Constr* clist, Branch* br, int t);
void CONSTR_list_eval_step(Constr* clist, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_list_store_sens_step(Constr* clist, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
Constr* CONSTR_new(Net* net);
void CONSTR_set_parameter(Constr* c, char* key, void* value);
void CONSTR_set_name(Constr* c, char* name);
void CONSTR_set_b(Constr* c, Vec* b);
void CONSTR_set_A(Constr* c, Mat* A);
void CONSTR_set_l(Constr* c, Vec* l);
void CONSTR_set_u(Constr* c, Vec* u);
void CONSTR_set_l_extra_vars(Constr* c, Vec* l);
void CONSTR_set_u_extra_vars(Constr* c, Vec* u);
void CONSTR_set_init_extra_vars(Constr* c, Vec* init);
void CONSTR_set_G(Constr* c, Mat* G);
void CONSTR_set_f(Constr* c, Vec* f);
void CONSTR_set_J(Constr* c, Mat* J);
void CONSTR_set_H_array(Constr* c, Mat* H_array, int size);
void CONSTR_set_H_combined(Constr* c, Mat* H_combined);
void CONSTR_set_H_single(Constr* c, int i, Mat* m);
void CONSTR_set_A_nnz(Constr* c, int nnz);
void CONSTR_set_G_nnz(Constr* c, int nnz);
void CONSTR_set_J_nnz(Constr* c, int nnz);
void CONSTR_set_H_nnz(Constr* c, int* nnz, int size);
void CONSTR_set_A_row(Constr* c, int index);
void CONSTR_set_G_row(Constr* c, int index);
void CONSTR_set_J_row(Constr* c, int index);
void CONSTR_set_bus_counted(Constr* c, char* counted, int size);
void CONSTR_set_data(Constr* c, void* data);
void CONSTR_set_A_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time);
void CONSTR_set_J_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time);
void CONSTR_set_G_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time);
void CONSTR_init(Constr* c);
void CONSTR_count(Constr* c);
void CONSTR_count_step(Constr* c, Branch* br, int t);
void CONSTR_allocate(Constr* c);
void CONSTR_clear(Constr* c);
void CONSTR_analyze(Constr* c);
void CONSTR_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_eval(Constr* c, Vec* v, Vec* ve);
void CONSTR_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_store_sens(Constr* c, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
BOOL CONSTR_is_safe_to_count(Constr* c);
BOOL CONSTR_is_safe_to_analyze(Constr* c);
BOOL CONSTR_is_safe_to_eval(Constr* c, Vec* v, Vec* ve);
void CONSTR_set_error(Constr* c, char* string);
void CONSTR_clear_error(Constr* c);
BOOL CONSTR_has_error(Constr* c);
char* CONSTR_get_error_string(Constr* c);
void CONSTR_update_network(Constr* c);
Net* CONSTR_get_network(Constr* c);
int CONSTR_get_num_extra_vars(Constr* c);
void CONSTR_set_num_extra_vars(Constr* c, int num);
void CONSTR_set_func_init(Constr* c, void (*func)(Constr* c));
void CONSTR_set_func_count_step(Constr* c, void (*func)(Constr* c, Branch* br, int t));
void CONSTR_set_func_allocate(Constr* c, void (*func)(Constr* c));
void CONSTR_set_func_clear(Constr* c, void (*func)(Constr* c));
void CONSTR_set_func_analyze_step(Constr* c, void (*func)(Constr* c, Branch* br, int t));
void CONSTR_set_func_eval_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* v, Vec* ve));
void CONSTR_set_func_store_sens_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl));
void CONSTR_set_func_free(Constr* c, void (*func)(Constr* c));
void CONSTR_set_func_set_parameter(Constr* c, void (*func)(Constr* c, char* key, void* value));

#endif
