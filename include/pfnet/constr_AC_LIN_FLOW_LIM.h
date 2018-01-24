/** @file constr_AC_LIN_FLOW_LIM.h
 *  @brief This file lists the constants and routines associated with the constraint of type AC_LIN_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
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
typedef struct LF_Results LF_Results;
typedef struct LF_Options LF_Options;

typedef struct LF_Branch { //contains branch parameters
  double V_i_min;
  double V_i_max;
  double V_j_min;
  double V_j_max;
  double g;
  double b;
  double b_sh;
  double t_ratio;
  double t_shift;
  double I_max;
} LF_Branch;

typedef enum LF_ResultFlag {
  non_binding,
  infeasible,
  success,
  error_branch_data, 
  error_options,
  zero_limit,
  error_other
} LF_ResultFlag;

void LF_set_branch_parameters(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
			      double g, double b, double b_sh, double t_ratio, double t_shift,
			      double I_max, LF_Branch* branch);
LF_Results* LF_construct(LF_Branch* branch, int flow_side, LF_Options* options);
void LF_free_results(LF_Results* results);
double* LF_get_A_matrix(LF_Results* results);
double* LF_get_b_vector(LF_Results* results);
LF_ResultFlag LF_get_flag(LF_Results* results);
int LF_get_number_constraints(LF_Results* results);
double LF_get_error(LF_Results* results);
char* LF_get_message(LF_Results* results);

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
