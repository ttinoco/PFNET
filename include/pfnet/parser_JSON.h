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
#include <pfnet/parser.h>
#include <pfnet/json.h>

// Structs
typedef struct JSON_Parser JSON_Parser;

// Interface
Parser* JSON_PARSER_new(void);
void JSON_PARSER_init(Parser* p, BOOL init_params);
Net* JSON_PARSER_parse(Parser* p, char* filename, int num_periods);
void JSON_PARSER_set(Parser* p, char* key, REAL value);
void JSON_PARSER_show(Parser* p);
void JSON_PARSER_write(Parser* p, Net* net, char* filename);
void JSON_PARSER_free(Parser* p, BOOL del_parser);

// Others
void JSON_PARSER_process_json_bus_array(Parser* p, Net* net, json_value* json_bus_array);
void JSON_PARSER_process_json_branch_array(Parser* p, Net* net, json_value* json_branch_array);
void JSON_PARSER_process_json_gen_array(Parser* p, Net* net, json_value* json_gen_array);
void JSON_PARSER_process_json_vargen_array(Parser* p, Net* net, json_value* json_vargen_array);
void JSON_PARSER_process_json_shunt_array(Parser* p, Net* net, json_value* json_shunt_array);
void JSON_PARSER_process_json_load_array(Parser* p, Net* net, json_value* json_load_array);
void JSON_PARSER_process_json_bat_array(Parser* p, Net* net, json_value* json_bat_array);

#endif
