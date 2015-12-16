/** @file parser_CSV.h
 *  @brief This file lists the constants and routines associated with the CSV_Parser data structures.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __PARSER_CSV_HEADER__
#define __PARSER_CSV_HEADER__

#include <stdio.h>
#include <string.h>
#include "net.h"
#include "types.h"

#define CSV_PARSER_BUFFER_SIZE 1024

// Structs
typedef struct CSV_Parser CSV_Parser;

// Prototypes
void CSV_PARSER_clear_field(CSV_Parser* p);
CSV_Parser* CSV_PARSER_new(void);
size_t CSV_PARSER_parse(CSV_Parser* p, 
			char* buffer,
			size_t len,
			BOOL last,
			char delimeter,
			char end_of_record,
			void (*cfield)(char*,void*),
			void (*crecord)(void*),
			void* data);
void CSV_PARSER_del(CSV_Parser* p);
			   
#endif
