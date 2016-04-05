/** @file func_LOAD_UTIL.h
 *  @brief This file lists the constants and routines associated with the function of type LOAD_UTIL.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __FUNC_LOAD_UTIL_HEADER__
#define __FUNC_LOAD_UTIL_HEADER__

#include <math.h>
#include "func.h"

// Function prototypes
void FUNC_LOAD_UTIL_init(Func* f);
void FUNC_LOAD_UTIL_count_branch(Func* f, Branch *branch);
void FUNC_LOAD_UTIL_allocate(Func* f);
void FUNC_LOAD_UTIL_clear(Func* f);
void FUNC_LOAD_UTIL_analyze_branch(Func* f, Branch *branch);
void FUNC_LOAD_UTIL_eval_branch(Func* f, Branch* branch, Vec* var_values);
void FUNC_LOAD_UTIL_free(Func* f);

#endif
