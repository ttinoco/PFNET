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
#define ART_PARSER_STATE_LINE 2
#define ART_PARSER_STATE_TRANSFO 3
#define ART_PARSER_STATE_LTCV 4
#define ART_PARSER_STATE_TRFO 5
#define ART_PARSER_STATE_PSHIFTP 6

// Defaults
#define ART_PARSER_BASE_POWER 100

// Tokens
#define ART_BUS_TOKEN  "BUS"
#define ART_LINE_TOKEN  "LINE"
#define ART_TRANSFO_TOKEN  "TRANSFO"
#define ART_LTCV_TOKEN  "LTC-V"
#define ART_TRFO_TOKEN  "TRFO"
#define ART_PSHIFTP_TOKEN  "PSHIFT-P"

// Structs
typedef struct ART_Bus ART_Bus;
typedef struct ART_Line ART_Line;
typedef struct ART_Transfo ART_Transfo;
typedef struct ART_Ltcv ART_Ltcv;
typedef struct ART_Trfo ART_Trfo;
typedef struct ART_Pshiftp ART_Pshiftp;
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
void ART_PARSER_callback_record(void* data);
void ART_PARSER_parse_bus_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_bus_record(ART_Parser* parser);
void ART_PARSER_parse_line_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_line_record(ART_Parser* parser);
void ART_PARSER_parse_transfo_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_transfo_record(ART_Parser* parser);
void ART_PARSER_parse_ltcv_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_ltcv_record(ART_Parser* parser);
void ART_PARSER_parse_trfo_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_trfo_record(ART_Parser* parser);
void ART_PARSER_parse_pshiftp_field(char* s, ART_Parser* parser);
void ART_PARSER_parse_pshiftp_record(ART_Parser* parser);

#endif
