/** @file func_SLIM_VMAG.h
 *  @brief This file lists the constants and routines associated with the function of type SLIM_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_SLIM_VMAG_HEADER__
#define __FUNC_SLIM_VMAG_HEADER__

#include <math.h>
#include "func.h"

#define FUNC_SLIM_VMAG_PARAM 1e-4

// Function prototypes
void FUNC_SLIM_VMAG_init(Func* f);
void FUNC_SLIM_VMAG_count_branch(Func* f, Branch *branch);
void FUNC_SLIM_VMAG_allocate(Func* f);
void FUNC_SLIM_VMAG_clear(Func* f);
void FUNC_SLIM_VMAG_analyze_branch(Func* f, Branch *branch);
void FUNC_SLIM_VMAG_eval_branch(Func* f, Branch* branch, Vec* var_values);
void FUNC_SLIM_VMAG_free(Func* f);

#endif
