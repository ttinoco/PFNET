/** @file utils.c
 *  @brief This file defines utility functions.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/utils.h>

char* trim(char* s) {
  /* Trims string inplace. */
  
  char* ptr;
  if (!s)
    return NULL;   // handle NULL string
  if (!*s)
    return s;      // handle empty string
  for (ptr = s + strlen(s) - 1; (ptr >= s) && isspace(*ptr); --ptr);
  ptr[1] = '\0';
  return s;
}

char* strtoupper(char s[]) {
  /* Replaces lowercase string chars with uppercase values. */
  int i = 0;
  while(s[i]) {
    s[i] = toupper(s[i]);
    i = i + 1;
  }
  return s;
}
