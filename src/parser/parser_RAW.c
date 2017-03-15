/** @file parser_RAW.c
 *  @brief This file defines the (dummy) RAW_Parser data structure and its associated (empty) methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_RAW.h>

#if !HAVE_RAW_PARSER

struct RAW_Parser;

RAW_Parser* RAW_PARSER_new(void) {
  return NULL;
}

void RAW_PARSER_read(RAW_Parser* parser, char* filename) {
  // pass
}

void RAW_PARSER_show(RAW_Parser* parser) {
  // pass
}

void RAW_PARSER_load(RAW_Parser* parser, Net* net) {
  // pass
}

void RAW_PARSER_del(RAW_Parser* parser) {
  // pass
}

void RAW_PARSER_set(RAW_Parser* parser, char* key, REAL value) {
  // pass
}

BOOL RAW_PARSER_has_error(RAW_Parser* parser) {
  return TRUE;
}

char* RAW_PARSER_get_error_string(RAW_Parser* parser) {
  return "this RAW parser does not do anything";
}

#endif
