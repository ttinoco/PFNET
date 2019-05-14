/** @file parser_EPC.h
 *  @brief This file list the constants and routines associated with the EPC_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_EPC_HEADER__
#define __PARSER_EPC_HEADER__

#include "parser.h"
#include "pfnet_config.h"

#ifndef HAVE_EPC_PARSER
#define HAVE_EPC_PARSER 0
#endif

// Interface
Parser* EPC_PARSER_new(void);
void EPC_PARSER_init(Parser* p, BOOL init_params);
Net* EPC_PARSER_parse(Parser* p, char* f, int num_periods);
void EPC_PARSER_set(Parser* p, char* key, REAL value);
void EPC_PARSER_show(Parser* p);
void EPC_PARSER_write(Parser* p, Net* net, char* f);
void EPC_PARSER_free(Parser* p, BOOL del_parser);

#endif
