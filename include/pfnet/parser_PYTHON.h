/** @file parser_PYTHON.h
 *  @brief This file lists the routines associated with the PYTHON parser interface.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_PYTHON_HEADER__
#define __PARSER_PYTHON_HEADER__

#include "net.h"

// Buffer
#define PYTHON_PARSER_BUFFER_SIZE 1024

// Structure
typedef struct PYTHON_Parser PYTHON_Parser;

// Prototypes
PYTHON_Parser* PYTHON_PARSER_new(char* classname);
void PYTHON_PARSER_read(PYTHON_Parser* parser, char* filename);
void PYTHON_PARSER_show(PYTHON_Parser* parser);
void PYTHON_PARSER_load(PYTHON_Parser* parser, Net* net);
void PYTHON_PARSER_del(PYTHON_Parser* parser);
BOOL PYTHON_PARSER_has_error(PYTHON_Parser* parser);
char* PYTHON_PARSER_get_error_string(PYTHON_Parser* parser);

#endif
