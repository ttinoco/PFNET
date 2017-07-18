/** @file parser_CSV.c
 *  @brief This file defines the CSV_Parser data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_CSV.h>

struct CSV_Parser {
  char field[CSV_PARSER_BUFFER_SIZE];
  int field_index;
  BOOL in_string_single;
  BOOL in_string_double;
  BOOL in_comment;
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
  csv->in_comment = FALSE;
  return csv;
}

size_t CSV_PARSER_parse(CSV_Parser* p,
			char* buffer,
			size_t len,
			BOOL last,
			char delimeter,
			char end_of_record,
			char comment,
			void (*cfield)(char*,void*),
			void (*crecord)(void*),
			void* data) {

  // Local variables
  size_t buffer_index;
  size_t lookahead_index;

  // Parse
  buffer_index = 0;
  while (buffer_index < len) {

    // Single quote
    if (buffer[buffer_index] == '\'' && !p->in_string_double && !p->in_comment) {
      if (!p->in_string_single)
        p->in_string_single = TRUE;
      else
        p->in_string_single = FALSE;
        // account for quote within a quote
        lookahead_index = buffer_index+1;
        while (lookahead_index < len) {
          if (buffer[lookahead_index] == delimeter)
            break;
          else if (buffer[lookahead_index] == '\'') {
            p->in_string_single = TRUE;
            break;
          }
          else
            lookahead_index++;
        }
    }

    // Double quote
    else if (buffer[buffer_index] == '"' && !p->in_string_single && !p->in_comment) {
      if (!p->in_string_double)
        p->in_string_double = TRUE;
      else
        p->in_string_double = FALSE;
        // account for quote within a quote
        lookahead_index = buffer_index+1;
        while (lookahead_index < len) {
          if (buffer[lookahead_index] == delimeter)
            break;
          else if (buffer[lookahead_index] == '"') {
            p->in_string_double = TRUE;
            break;
          }
          else
            lookahead_index++;
        }
    }

    // Comment
    else if (buffer[buffer_index] == comment && !p->in_string_single && !p->in_string_double) {
      p->in_comment = TRUE;
    }

    // End of field
    else if (buffer[buffer_index] == delimeter && !p->in_string_single && !p->in_string_double && !p->in_comment) {
      if (p->field_index > 0) {
        p->field[p->field_index] = 0;
        cfield(p->field,data);
        p->field_index = 0;
      }

      // Skip remaining white if white is delimeter
      if (delimeter == ' ') {
        while ((buffer[buffer_index] == ' ' ||
                buffer[buffer_index] == '\t') &&
                buffer_index < len)
          buffer_index++;
        buffer_index--;
      }
    }

    // End of record
    else if (buffer[buffer_index] == end_of_record) {
      p->in_string_single = FALSE;
      p->in_string_double = FALSE;
      p->in_comment = FALSE;
      if (p->field_index > 0) {
        p->field[p->field_index] = 0;
        cfield(p->field,data);
        p->field_index = 0;
      }
      crecord(data);
    }

    // End of line or return carriage (treat as end of field).
    // Return carriage also supports dos/windows end of line '\r\n'
    else if ((buffer[buffer_index] == '\n') ||
	     (buffer[buffer_index] == '\r')) {
      p->in_comment = FALSE;
      if (p->field_index > 0) {
        p->field[p->field_index] = 0;
        cfield(p->field,data);
        p->field_index = 0;
      }
    }

    // In field
    else if (!p->in_comment) {
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
