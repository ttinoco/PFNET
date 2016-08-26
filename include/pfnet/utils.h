/** @file utils.h
 *  @brief This file some utiliy functions.
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

static char* trim(char *s) {
  /* Trims string inplace. */

  char *ptr;
  if (!s)
    return NULL;   // handle NULL string
  if (!*s)
    return s;      // handle empty string
  for (ptr = s + strlen(s) - 1; (ptr >= s) && isspace(*ptr); --ptr);
  ptr[1] = '\0';
  return s;
}

#endif
