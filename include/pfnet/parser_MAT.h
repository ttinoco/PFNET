/** @file parser_MAT.h
 *  @brief This file list the constants and routines associated with the MAT_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_MAT_HEADER__
#define __PARSER_MAT_HEADER__

#include <stdio.h>
#include <string.h>
#include "parser.h"
#include "parser_CSV.h"

// Buffer
#define MAT_PARSER_BUFFER_SIZE 1024
#define MAT_BUS_NAME_BUFFER_SIZE 25

// State
#define MAT_PARSER_STATE_TOKEN -1
#define MAT_PARSER_STATE_TITLE 0
#define MAT_PARSER_STATE_BUS 1
#define MAT_PARSER_STATE_GEN 2
#define MAT_PARSER_STATE_BRANCH 3
#define MAT_PARSER_STATE_COST 4
#define MAT_PARSER_STATE_UTIL 5

// Defaults
#define MAT_PARSER_BASE_POWER 100

// Tokens
#define MAT_BUS_TOKEN  "BUS"
#define MAT_END_TOKEN  "END"
#define MAT_GEN_TOKEN  "GEN"
#define MAT_BRANCH_TOKEN "BRANCH"
#define MAT_COST_TOKEN "COST"
#define MAT_UTIL_TOKEN "UTIL"

// Bus types
#define MAT_BUS_TYPE_PQ 1
#define MAT_BUS_TYPE_PV 2
#define MAT_BUS_TYPE_SL 3
#define MAT_BUS_TYPE_IS 4

// Structs
typedef struct MAT_Bus MAT_Bus;
typedef struct MAT_Gen MAT_Gen;
typedef struct MAT_Branch MAT_Branch;
typedef struct MAT_Parser MAT_Parser;
typedef struct MAT_Cost MAT_Cost;
typedef struct MAT_Util MAT_Util;

// Interface
Parser* MAT_PARSER_new(void);
void MAT_PARSER_init(Parser* p, BOOL init_params);
Net* MAT_PARSER_parse(Parser* p, char* f, int num_periods);
void MAT_PARSER_set(Parser* p, char* key, REAL value);
void MAT_PARSER_show(Parser* p);
void MAT_PARSER_write(Parser* p, Net* net, char* f);
void MAT_PARSER_free(Parser* p, BOOL del_parser);

// MAT-specific
void MAT_PARSER_load(MAT_Parser* p, Net* net);
void MAT_PARSER_clear_token(MAT_Parser* p);
BOOL MAT_PARSER_has_error(MAT_Parser* p);
char* MAT_PARSER_get_error_string(MAT_Parser* p);
void MAT_PARSER_callback_field(char* s, void* data);
void MAT_PARSER_callback_row(void* data);
void MAT_PARSER_parse_token_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_token_row(MAT_Parser* p);
void MAT_PARSER_parse_title_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_title_row(MAT_Parser* p);
void MAT_PARSER_parse_bus_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_bus_row(MAT_Parser* p);
void MAT_PARSER_parse_gen_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_gen_row(MAT_Parser* p);
void MAT_PARSER_parse_branch_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_branch_row(MAT_Parser* p);
void MAT_PARSER_parse_cost_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_cost_row(MAT_Parser* p);
void MAT_PARSER_parse_util_field(char* s, MAT_Parser* p);
void MAT_PARSER_parse_util_row(MAT_Parser* p);

#endif
