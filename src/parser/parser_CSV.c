/** @file parser_CSV.c
 *  @brief This file defines the CSV_Parser data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_CSV.h>

struct CSV_Parser {
  char field[CSV_PARSER_BUFFER_SIZE];
  int field_index;
  BOOL in_string_single;
  BOOL in_string_double;
};

void CSV_PARSER_clear_field(CSV_Parser* p) {
  int i;
  if (p) {
    for (i = 0; i < CSV_PARSER_BUFFER_SIZE; i++)
      p->field[i] = 0;
  }
}

CSV_Parser* CSV_PARSER_new(void) {
  CSV_Parser* csv = (CSV_Parser*)malloc(sizeof(CSV_Parser));
  CSV_PARSER_clear_field(csv);
  csv->field_index = 0;
  csv->in_string_single = FALSE;
  csv->in_string_double = FALSE;
  return csv;
}

size_t CSV_PARSER_parse(CSV_Parser* p, 
			char* buffer,
			size_t len,
			BOOL last,
			char delimeter,
			char end_of_record,
			void (*cfield)(char*,void*),
			void (*crecord)(void*),
			void* data) {
  
  // Local variables
  size_t buffer_index;

  // Parse
  buffer_index = 0;
  while (buffer_index < len) {

    // single quote
    if (buffer[buffer_index] == '\'' && !p->in_string_double) {
      if (!p->in_string_single)
	p->in_string_single = TRUE;
      else
	p->in_string_single = FALSE;
    }

    // double quote
    else if (buffer[buffer_index] == '"' && !p->in_string_single) {
      if (!p->in_string_double)
	p->in_string_double = TRUE;
      else
	p->in_string_double = FALSE; 
    }
	
    // end of field
    else if (buffer[buffer_index] == delimeter && !p->in_string_single && !p->in_string_double) {
      p->field[p->field_index] = 0;
      cfield(p->field,data);
      p->field_index = 0;
      
      // skip remaining white if white is delimeter
      if (delimeter == ' ') {
	while ((buffer[buffer_index] == ' ' || 
		buffer[buffer_index] == '\t') &&
	       buffer_index < len)
	  buffer_index++;
	buffer_index--;
      }
    }
    
    // end of record
    else if (buffer[buffer_index] == end_of_record) {
      p->in_string_single = FALSE;
      p->in_string_double = FALSE;
      if (p->field_index > 0) {
	p->field[p->field_index] = 0;
	cfield(p->field,data);
	p->field_index = 0;
      }
      crecord(data);
    }

    // end of line
    else if (buffer[buffer_index] == '\n') {
      if (p->field_index > 0) {
	p->field[p->field_index] = 0;
	cfield(p->field,data);
	p->field_index = 0;
      }
    }
    
    // in field
    else {
      p->field[p->field_index] = buffer[buffer_index];
      p->field_index++;
    }

    buffer_index++;
  }

  // Remaining
  if (last && p->field_index > 0) {
    p->field[p->field_index] = 0;
    cfield(p->field,data);
    p->field_index = 0;
    crecord(data);
  }
  
  return buffer_index;
}

void CSV_PARSER_del(CSV_Parser* p) {
  if (p)
    free(p);
}

