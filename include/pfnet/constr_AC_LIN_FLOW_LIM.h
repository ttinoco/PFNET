/** @file constr_AC_LIN_FLOW_LIM.h
 *  @brief This file lists the constants and routines associated with the constraint of type AC_LIN_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __CONSTR_AC_LIN_FLOW_LIM_HEADER__
#define __CONSTR_AC_LIN_FLOW_LIM_HEADER__

#include <math.h>
#include "constr.h"
#include "pfnet_config.h"

#define CONSTR_AC_LIN_FLOW_LIM_INF 1e8

// Line flow interface
typedef struct LINE_FLOW_Results LINE_FLOW_Results;
typedef struct LINE_FLOW_Params LINE_FLOW_Params;
typedef enum { error_line_data, error_params, error_patch, infeasible, 
	       non_binding, success, zero_limit, error_other } LINE_FLOW_Flag;

LINE_FLOW_Results* LINE_FLOW_construct(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
				       double g, double b, double B_sh, double Kt_real, double Kt_shift, 
				       double I_max_user, int flow_side,
				       LINE_FLOW_Params* params);

void LINE_FLOW_free_results(LINE_FLOW_Results* results);
double* LINE_FLOW_get_A_matrix(LINE_FLOW_Results* results);
double* LINE_FLOW_get_b_vector(LINE_FLOW_Results* results);
LINE_FLOW_Flag LINE_FLOW_get_flag(LINE_FLOW_Results* results);
int LINE_FLOW_get_number_constraints(LINE_FLOW_Results* results);
double LINE_FLOW_get_error(LINE_FLOW_Results* results);
char* LINE_FLOW_get_message(LINE_FLOW_Results* results);

// Data
typedef struct Constr_AC_LIN_FLOW_LIM_Data Constr_AC_LIN_FLOW_LIM_Data;

// Function prototypes
Constr* CONSTR_AC_LIN_FLOW_LIM_new(Net* net);
void CONSTR_AC_LIN_FLOW_LIM_init(Constr* c);
void CONSTR_AC_LIN_FLOW_LIM_count_step(Constr* c, Branch* br, int t);
void CONSTR_AC_LIN_FLOW_LIM_allocate(Constr* c);
void CONSTR_AC_LIN_FLOW_LIM_clear(Constr* c);
void CONSTR_AC_LIN_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t);
void CONSTR_AC_LIN_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve);
void CONSTR_AC_LIN_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);
void CONSTR_AC_LIN_FLOW_LIM_free(Constr* c);

#endif
