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

  char name[10];
  REAL vnom;    // kv
  REAL pload;   // mw
  REAL qload;   // mvar
  REAL bshunt;  // mvar
  REAL qshunt;  // mvar 
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

  // Base
  REAL base_power;
  
  // Buses
  ART_Bus* bus;
  ART_Bus* bus_list;
  ART_Bus* bus_hash;
};

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
    if (CSV_PARSER_parse(csv,
			 buffer,
			 bytes_read,
			 feof(file),
			 ' ',
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

  // Local variables
  int len_bus_list;

  if (!parser)
    return;

  // List lengths
  LIST_len(ART_Bus,parser->bus_list,next,len_bus_list);

  // Show
  printf("\nParsed Data\n");
  printf("base power : %.2f\n",parser->base_power);
  printf("bus list   : %d\n",len_bus_list);

  // Debugging
  ART_Bus* bus;
  for (bus = parser->bus_list; bus != NULL; bus = bus->next) {
    printf("Bus %s %.5f %.5f %.5f %.5f %.5f\n",
	   bus->name,
	   bus->vnom,
	   bus->pload,
	   bus->qload,
	   bus->bshunt,
	   bus->qshunt);
  }
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

  // Local variables
  ART_Parser* parser = (ART_Parser*)data;

  // Field callback
  switch (parser->state) {
    
  case ART_PARSER_STATE_INIT:
    if (parser->field == 0) {
      
      // Bus
      if (strstr(s,ART_BUS_TOKEN) != NULL) {
	printf("*** BUS STATE ***\n");
	parser->state = ART_PARSER_STATE_BUS;
      }
	
    }
    break;
    
  case ART_PARSER_STATE_BUS:
    ART_PARSER_parse_bus_field((char*)s,parser);
    break;
  }
    
  // Update field
  parser->field++;
}

void ART_PARSER_callback_row(void *data) {

  // Local variables
  ART_Parser* parser = (ART_Parser*)data;

  // Row callback
  switch (parser->state) {
  case ART_PARSER_STATE_INIT:
    parser->field = 0;
    parser->record = 0;
    break;
  case ART_PARSER_STATE_BUS:
    ART_PARSER_parse_bus_row(parser);
    break;
  }  
}

void ART_PARSER_parse_bus_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New bus
  if (parser->field == 1) {
    parser->bus = (ART_Bus*)malloc(sizeof(ART_Bus));
    parser->bus->next = NULL;
  }

  // Fields
  if (parser->bus) {
    switch (parser->field) {
    case 1:
      strcpy(parser->bus->name,s);
      break;
    case 2:
      parser->bus->vnom = atof(s);
      break;
    case 3:
      parser->bus->pload = atof(s);
      break;
    case 4:
      parser->bus->qload = atof(s);
      break;
    case 5:
      parser->bus->bshunt = atof(s);
      break;
    case 6:
      parser->bus->qshunt = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_bus_row(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->bus && parser->record >= 0) {
    LIST_add(parser->bus_list,parser->bus,next);
    HASH_ADD_STR(parser->bus_hash,name,parser->bus);
  }
  parser->bus = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}
