/** @file problem.h
 *  @brief This file lists the constants and routines associated with the Prob data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PROBLEM_HEADER__
#define __PROBLEM_HEADER__

#include "net.h"
#include "constr.h"
#include "func.h"
#include "heur.h"

// Buffer
#define PROB_BUFFER_SIZE 1024 /**< @brief Default problem buffer size for strings */

// Inf
#define PROB_EXTRA_VAR_INF 1e8 /**< @brief Large constant for lower and upper bounds */

// Problem
typedef struct Prob Prob;

// Function prototypes
void PROB_add_constr(Prob* p, int type);
void PROB_add_func(Prob* p, int type, REAL weight);
void PROB_add_heur(Prob* p, int type);
void PROB_analyze(Prob* p);
void PROB_apply_heuristics(Prob* p, Vec* point);
void PROB_eval(Prob* p, Vec* point);
void PROB_store_sens(Prob* p, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void PROB_del(Prob* p);
void PROB_del_matvec(Prob* p);
void PROB_clear(Prob* p);
void PROB_clear_error(Prob* p);
void PROB_combine_H(Prob* p, Vec* coeff, BOOL ensure_psd);
Constr* PROB_find_constr(Prob* p, int constr_type);
Constr* PROB_get_constr(Prob* p);
char* PROB_get_error_string(Prob* p);
Func* PROB_get_func(Prob* p);
Heur* PROB_get_heur(Prob* p);
Vec* PROB_get_init_point(Prob* p);
Vec* PROB_get_upper_limits(Prob* p);
Vec* PROB_get_lower_limits(Prob* p);
Net* PROB_get_network(Prob* p);
REAL PROB_get_phi(Prob* p);
Vec* PROB_get_gphi(Prob* p);
Mat* PROB_get_Hphi(Prob* p);
Vec* PROB_get_b(Prob* p);
Mat* PROB_get_A(Prob* p);
Vec* PROB_get_l(Prob* p);
Vec* PROB_get_u(Prob* p);
Mat* PROB_get_G(Prob* p);
Vec* PROB_get_f(Prob* p);
Mat* PROB_get_J(Prob* p);
Mat* PROB_get_H_combined(Prob* p);
BOOL PROB_has_error(Prob* p);
void PROB_init(Prob* p);
Prob* PROB_new(void);
void PROB_set_network(Prob* p, Net* net);
void PROB_show(Prob* p);
char* PROB_get_show_str(Prob* p);
void PROB_update_lin(Prob* p);
void PROB_update_nonlin_struc(Prob* p);
void PROB_update_nonlin_data(Prob* p, Vec* point);
int PROB_get_num_primal_variables(Prob* p);
int PROB_get_num_linear_equality_constraints(Prob* p);
int PROB_get_num_nonlinear_equality_constraints(Prob* p);
int PROB_get_num_extra_vars(Prob* p);

#endif
