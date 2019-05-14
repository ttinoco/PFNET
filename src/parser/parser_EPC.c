/** @file parser_EPC.c
 *  @brief This file defines the (dummy) EPC_Parser data structure and its associated (empty) methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_EPC.h>

#if !HAVE_EPC_PARSER

Parser* EPC_PARSER_new(void) {
  return NULL;
}

void EPC_PARSER_init(Parser* p, BOOL params) {
  // pass
}

Net* EPC_PARSER_parse(Parser* p, char* filename, int num_periods) {
  return NULL;
}

void EPC_PARSER_set(Parser* p, char* key, REAL value) {
  // pass
}

void EPC_PARSER_show(Parser* p) {
  // pass
}

void EPC_PARSER_write(Parser* p, Net* net, char* f) {
  // nothing
}

void EPC_PARSER_free(Parser* p, BOOL del_params) {
  // pass
}

#endif
