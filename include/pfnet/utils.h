/** @file utils.h
 *  @brief This file lists utility functions.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __UTILS_HEADER__
#define __UTILS_HEADER__

#include <stdio.h>
#include <string.h>
#include <ctype.h>

char* trim(char* s);
char* strtoupper(char s[]);
char* strtolower(char s[]);

#endif
