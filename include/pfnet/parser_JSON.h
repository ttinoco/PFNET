/** @file parser_JSON.h
 *  @brief This file list the constants and routines associated with the JSON_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_JSON_HEADER__
#define __PARSER_JSON_HEADER__

#include <stdio.h>
#include <string.h>
#include "parser.h"

// Structs
typedef struct JSON_Parser JSON_Parser;

// Interface
Parser* JSON_PARSER_new(void);
void JSON_PARSER_init(Parser* p);
Net* JSON_PARSER_parse(Parser* p, char* f, int num_periods);
void JSON_PARSER_set(Parser* p, char* key, REAL value);
void JSON_PARSER_show(Parser* p);
void JSON_PARSER_write(Parser* p, Net* net, char* f);
void JSON_PARSER_free(Parser* p);

#endif
