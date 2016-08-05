/** @file parser_RAW.h
 *  @brief This file list the constants and (dummy) routines associated with the RAW_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_RAW_HEADER__
#define __PARSER_RAW_HEADER__

#include "net.h"

// Struct
typedef struct RAW_Parser RAW_Parser;

// Prototypes
RAW_Parser* RAW_PARSER_new(void);
void RAW_PARSER_read(RAW_Parser* parser, char* filename);
void RAW_PARSER_show(RAW_Parser* parser);
void RAW_PARSER_load(RAW_Parser* parser, Net* net);
void RAW_PARSER_del(RAW_Parser* parser);
void RAW_PARSER_set(RAW_Parser* parser, char* key, REAL value);
BOOL RAW_PARSER_has_error(RAW_Parser* parser);
char* RAW_PARSER_get_error_string(RAW_Parser* parser);

#endif
