/** @file parser_JSON.c
 *  @brief This file defines the JSON_Parser data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/json.h>
#include <pfnet/parser_JSON.h>

Parser* JSON_PARSER_new(void) {
  return NULL;
}

void JSON_PARSER_init(Parser* p) {
  // pass
}

Net* JSON_PARSER_parse(Parser* p, char* filename, int num_periods) {
  return NULL;
}

void JSON_PARSER_set(Parser* p, char* key, REAL value) {
  // pass
}

void JSON_PARSER_show(Parser* p) {
  // pass
}

void JSON_PARSER_write(Parser* p, Net* net, char* filename) {

  // Local variables
  FILE* file;
  char* json_string;

  // Open file
  file = fopen(filename,"w");
  if (file == NULL) {
    PARSER_set_error(p,"unable to open file");
    return;
  }

  // Write
  json_string = NET_get_json_string(net);
  fprintf(file,"%s",json_string);

  // Clean up
  free(json_string);
  fclose(file);
}

void JSON_PARSER_free(Parser* p) {
  // pass
}
