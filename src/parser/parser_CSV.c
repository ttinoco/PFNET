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
  return csv;
}

size_t CSV_PARSER_parse(CSV_Parser* p, 
			char* buffer,
			size_t len,
			BOOL last,
			char del,
			void (*cfield)(char*,void*),
			void (*crow)(void*),
			void* data) {
  
  // Local variables
  size_t buffer_index;

  buffer_index = 0;
  while (buffer_index < len) {
    
    // end of field
    if (buffer[buffer_index] == del) {
      p->field[p->field_index] = 0;
      cfield(p->field,data);
      p->field_index = 0;
      
      // skip remaining white if white is del
      if (del == ' ') {
	while ((buffer[buffer_index] == ' ' || 
		buffer[buffer_index] == '\t') &&
	       buffer_index < len)
	  buffer_index++;
	buffer_index--;
      }
    }
    
    // end of line
    else if (buffer[buffer_index] == '\n') {
      if (p->field_index > 0) {
	p->field[p->field_index] = 0;
	cfield(p->field,data);
	p->field_index = 0;
      }
      crow(data);
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
    crow(data);
  }
  
  return buffer_index;
}

void CSV_PARSER_del(CSV_Parser* p) {
  if (p)
    free(p);
}

