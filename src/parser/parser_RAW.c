/** @file parser_RAW.c
 *  @brief This file defines the (dummy) RAW_Parser data structure and its associated (empty) methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_RAW.h>

#if !HAVE_RAW_PARSER

Parser* RAW_PARSER_new(void) {
  return NULL;
}

void RAW_PARSER_init(Parser* p, BOOL params) {
  // pass
}

Net* RAW_PARSER_parse(Parser* p, char* filename, int num_periods) {
  return NULL;
}

void RAW_PARSER_set(Parser* p, char* key, REAL value) {
  // pass
}

void RAW_PARSER_show(Parser* p) {
  // pass
}

void RAW_PARSER_write(Parser* p, Net* net, char* f) {
  // nothing
}

void RAW_PARSER_free(Parser* p, BOOL del_params) {
  // pass
}

#endif
