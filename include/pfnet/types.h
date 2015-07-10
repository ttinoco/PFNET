/** @file types.h
 *  @brief This files defines scalar types.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __TYPES_HEADER__
#define __TYPES_HEADER__

#include <stdlib.h>
#include <stddef.h>
#include "flag_types.h"
#include "obj_types.h"

// Real
typedef double REAL;

// Bool
typedef char BOOL;
enum { FALSE = 0, TRUE = 1};

#endif
