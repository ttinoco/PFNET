/** @file parser_ART.c
 *  @brief This file defines the ART_Parser and related data structures and their associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_ART.h>

struct ART_Bus {
  
  struct ART_Bus* next;
  UT_hash_handle hh;
};

struct ART_Parser {

  // Error
  BOOL error_flag;
  char error_string[ART_PARSER_BUFFER_SIZE];

  // State
  int state;
  int field;
  int record;
  char token[ART_PARSER_BUFFER_SIZE];

  // Base
  REAL base_power;
  
  // Buses
  ART_Bus* bus;
  ART_Bus* bus_list;
  ART_Bus* bus_hash;
};

void ART_PARSER_clear_token(ART_Parser* parser) {
  int i;
  for (i = 0; i < ART_PARSER_BUFFER_SIZE; i++)
    parser->token[i] = 0;
}

ART_Parser* ART_PARSER_new(void) {
  
  // Allocate
  ART_Parser* parser = (ART_Parser*)malloc(sizeof(ART_Parser));

  // Error
  parser->error_flag = FALSE;
  strcpy(parser->error_string,"");

  // State
  parser->state = ART_PARSER_STATE_INIT;
  parser->field = 0;
  parser->record = 0;
  ART_PARSER_clear_token(parser);

  // Base
  parser->base_power = ART_PARSER_BASE_POWER;

  // Buses
  parser->bus = NULL;
  parser->bus_list = NULL;
  parser->bus_hash = NULL;

  // Return
  return parser;
}

void ART_PARSER_read(ART_Parser* parser, char* filename) {
  
  // Local variables
  FILE* file;
  char* line;
  CSV_Parser* csv = CSV_PARSER_new();
  char buffer[ART_PARSER_BUFFER_SIZE];
  size_t bytes_read;

  // No parser
  if (!parser)
    return;
  
  // Open file
  file = fopen(filename,"rb");
  if (!file) {
    sprintf(parser->error_string,"unable to open file %s",filename);
    parser->error_flag = TRUE;
    CSV_PARSER_del(csv);
    return;
  }
  
  // Parse
  while ((bytes_read=fread(buffer,1,ART_PARSER_BUFFER_SIZE,file)) > 0) {
    if (CSV_PARSER_parse(csv,buffer,bytes_read,feof(file),
			 ART_PARSER_callback_field,
			 ART_PARSER_callback_row,
			 parser) != bytes_read) {
      strcpy(parser->error_string,"error parsing buffer");
      parser->error_flag = TRUE;
      break;
    }
  }
 
  // Free and close
  CSV_PARSER_del(csv);
  fclose(file);
}

void ART_PARSER_show(ART_Parser* parser) {

  if (!parser)
    return;

  // Show
  printf("\nParsed Data\n");
  printf("base power : %.2f\n",parser->base_power);
}

void ART_PARSER_load(ART_Parser* parser, Net* net) {


}

void ART_PARSER_del(ART_Parser* parser) {

  if (!parser)
    return;

  // Parser
  free(parser);  
}

BOOL ART_PARSER_has_error(ART_Parser* parser) {
  if (!parser)
    return FALSE;
  else
    return parser->error_flag;
}

char* ART_PARSER_get_error_string(ART_Parser* parser) {
  if (!parser)
    return NULL;
  else
    return parser->error_string;
}

void ART_PARSER_callback_field(char* s, void* data) {

}

void ART_PARSER_callback_row(void *data) {

}


