/** @file parser_RAW.h
 *  @brief This file list the constants and routines associated with the RAW_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_RAW_HEADER__
#define __PARSER_RAW_HEADER__

#include "parser.h"
#include "pfnet_config.h"

#ifndef HAVE_RAW_PARSER
#define HAVE_RAW_PARSER 0
#endif

// Interface
Parser* RAW_PARSER_new(void);
void RAW_PARSER_init(Parser* p);
Net* RAW_PARSER_parse(Parser* p, char* f, int num_periods);
void RAW_PARSER_set(Parser* p, char* key, REAL value);
void RAW_PARSER_show(Parser* p);
void RAW_PARSER_write(Parser* p, Net* net, char* f);
void RAW_PARSER_free(Parser* p);

#endif
