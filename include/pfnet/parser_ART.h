/** @file parser_ART.h
 *  @brief This file list the constants and routines associated with the ART_Parser data structure.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_ART_HEADER__
#define __PARSER_ART_HEADER__

#include <stdio.h>
#include <string.h>
#include "net.h"
#include "types.h"
#include "parser_CSV.h"

// Buffer
#define ART_PARSER_BUFFER_SIZE 1024

// State
#define ART_PARSER_STATE_INIT 0
#define ART_PARSER_STATE_BUS 1

// Defaults
#define ART_PARSER_BASE_POWER 100

// Tokens
#define ART_BUS_TOKEN  "BUS"

// Structs
typedef struct ART_Bus ART_Bus;
typedef struct ART_Parser ART_Parser;

// Prototypes
ART_Parser* ART_PARSER_new(void);
void ART_PARSER_clear_token(ART_Parser* parser);
void ART_PARSER_read(ART_Parser* parser, char* filename);
void ART_PARSER_show(ART_Parser* parser);
void ART_PARSER_load(ART_Parser* parser, Net* net);
void ART_PARSER_del(ART_Parser* parser);
BOOL ART_PARSER_has_error(ART_Parser* parser);
char* ART_PARSER_get_error_string(ART_Parser* parser);
void ART_PARSER_callback_field(char* s, void* data);
void ART_PARSER_callback_row(void* data);
void ART_PARSER_parse_bus_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_bus_row(ART_Parser* parser);

#endif
