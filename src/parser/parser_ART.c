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

struct ART_Line {

  char name[22];
  char from_bus[10];
  char to_bus[10];
  REAL r;             // ohms
  REAL x;             // ohms
  REAL wc_half;       // micro siemens
  REAL snom;          // mva
  REAL br;            
  struct ART_Line* next;
};

struct ART_Transfo {

  char name[22];
  char from_bus[10];
  char to_bus[10];
  REAL r;             // % on the Vb1,SNOM base
  REAL x;             // % on the Vb1,SNOM base
  REAL b1;            // % on the Vb1,SNOM base
  REAL b2;            // % on the Vb1,SNOM base
  REAL n;             // % on the Vb1,Vb2 base
  REAL phi;           // degrees
  REAL snom;          // mva
  REAL br;            
  struct ART_Transfo* next;
};

struct ART_Ltcv {

  char name[22];
  char con_bus[10];
  REAL nfirst;      // %           
  REAL nlast;       // %    
  int nbpos; 
  REAL tolv;        // per unit
  REAL vdes;        // per unit
 
  struct ART_Ltcv* next;
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

  // Lines
  ART_Line* line;
  ART_Line* line_list;

  // Transformers
  ART_Transfo* transfo;
  ART_Transfo* transfo_list;

  // LTC-Vs
  ART_Ltcv* ltcv;
  ART_Ltcv* ltcv_list;
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

  // Lines
  parser->line = NULL;
  parser->line_list = NULL;

  // Transformers
  parser->transfo = NULL;
  parser->transfo_list = NULL;

  // LTC-Vs
  parser->ltcv = NULL;
  parser->ltcv_list = NULL;

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
			 ';',
			 ART_PARSER_callback_field,
			 ART_PARSER_callback_record,
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
  int len_line_list;
  int len_transfo_list;
  int len_ltcv_list;

  if (!parser)
    return;

  // List lengths
  LIST_len(ART_Bus,parser->bus_list,next,len_bus_list);
  LIST_len(ART_Line,parser->line_list,next,len_line_list);
  LIST_len(ART_Transfo,parser->transfo_list,next,len_transfo_list);
  LIST_len(ART_Ltcv,parser->ltcv_list,next,len_ltcv_list);
  
  // Show
  printf("\nParsed Data\n");
  printf("base power   : %.2f\n",parser->base_power);
  printf("bus list     : %d\n",len_bus_list);
  printf("line list    : %d\n",len_line_list);
  printf("transfo list : %d\n",len_transfo_list);
  printf("ltc-v list   : %d\n",len_ltcv_list);

  // Debugging BUS
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

  // Debugging LINE
  ART_Line* line;
  for (line = parser->line_list; line != NULL; line = line->next) {
    printf("Line %s %s %s %.5f %.5f %.5f %.5f %.5f\n",
	   line->name,
	   line->from_bus,
	   line->to_bus,
	   line->r,
	   line->x,
	   line->wc_half,
	   line->snom,
	   line->br);
  }
  
  // Debugging TRANSFO
  ART_Transfo* transfo;
  for (transfo = parser->transfo_list; transfo != NULL; transfo = transfo->next) {
    printf("Transfo %s %s %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
	   transfo->name,
	   transfo->from_bus,
	   transfo->to_bus,
	   transfo->r,
	   transfo->x,
	   transfo->b1,
	   transfo->b2,
	   transfo->n,
	   transfo->phi,
	   transfo->snom,
	   transfo->br);
  }

  // Debugging LTC-V
  ART_Ltcv* ltcv;
  for (ltcv = parser->ltcv_list; ltcv != NULL; ltcv = ltcv->next) {
    printf("Ltcv %s %s %.5f %.5f %d %.5f %.5f\n",
	   ltcv->name,
	   ltcv->con_bus,
	   ltcv->nfirst,
	   ltcv->nlast,
	   ltcv->nbpos,
	   ltcv->tolv,
	   ltcv->vdes);
  }
}

void ART_PARSER_load(ART_Parser* parser, Net* net) {


}

void ART_PARSER_del(ART_Parser* parser) {

  if (!parser)
    return;

  // Buses
  while (parser->bus_hash)
    HASH_DEL(parser->bus_hash,parser->bus_hash);
  LIST_map(ART_Bus,parser->bus_list,bus,next,{free(bus);});

  // Line
  LIST_map(ART_Line,parser->line_list,line,next,{free(line);});

  // Transformers
  LIST_map(ART_Transfo,parser->transfo_list,transfo,next,{free(transfo);});

  // LTC-Vs
  LIST_map(ART_Ltcv,parser->ltcv_list,ltcv,next,{free(ltcv);});

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
      
      // Line
      else if (strstr(s,ART_LINE_TOKEN) != NULL) {
	printf("*** LINE STATE ***\n");
	parser->state = ART_PARSER_STATE_LINE;
      }

      // Transformer
      else if (strstr(s,ART_TRANSFO_TOKEN) != NULL) {
	printf("*** TRANSFO STATE ***\n");
	parser->state = ART_PARSER_STATE_TRANSFO;
      }

      // LTC-V
      else if (strstr(s,ART_LTCV_TOKEN) != NULL) {
	printf("*** LTCV STATE ***\n");
	parser->state = ART_PARSER_STATE_LTCV;
      }
	
    }
    break;
    
  case ART_PARSER_STATE_BUS:
    ART_PARSER_parse_bus_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_LINE:
    ART_PARSER_parse_line_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_TRANSFO:
    ART_PARSER_parse_transfo_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_LTCV:
    ART_PARSER_parse_ltcv_field((char*)s,parser);
    break;
  }
    
  // Update field
  parser->field++;
}

void ART_PARSER_callback_record(void *data) {

  // Local variables
  ART_Parser* parser = (ART_Parser*)data;

  // Record callback
  switch (parser->state) {
  case ART_PARSER_STATE_INIT:
    parser->field = 0;
    parser->record = 0;
    break;
  case ART_PARSER_STATE_BUS:
    ART_PARSER_parse_bus_record(parser);
    break;
  case ART_PARSER_STATE_LINE:
    ART_PARSER_parse_line_record(parser);
    break;
  case ART_PARSER_STATE_TRANSFO:
    ART_PARSER_parse_transfo_record(parser);
    break;
  case ART_PARSER_STATE_LTCV:
    ART_PARSER_parse_ltcv_record(parser);
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

void ART_PARSER_parse_bus_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->bus) {
    LIST_add(parser->bus_list,parser->bus,next);
    HASH_ADD_STR(parser->bus_hash,name,parser->bus);
  }
  parser->bus = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_line_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New line
  if (parser->field == 1) {
    parser->line = (ART_Line*)malloc(sizeof(ART_Line));
    parser->line->next = NULL;
  }

  // Fields
  if (parser->line) {
    switch (parser->field) {
    case 1:
      strcpy(parser->line->name,s);
      break;
    case 2:
      strcpy(parser->line->from_bus,s);
      break;
    case 3:
      strcpy(parser->line->to_bus,s);
      break;
    case 4:
      parser->line->r = atof(s);
      break;
    case 5:
      parser->line->x = atof(s);
      break;
    case 6:
      parser->line->wc_half = atof(s);
      break;
    case 7:
      parser->line->snom = atof(s);
      break;
    case 8:
      parser->line->br = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_line_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->line) {
    LIST_add(parser->line_list,parser->line,next);
  }
  parser->line = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_transfo_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New transfo
  if (parser->field == 1) {
    parser->transfo = (ART_Transfo*)malloc(sizeof(ART_Transfo));
    parser->transfo->next = NULL;
  }

  // Fields
  if (parser->transfo) {
    switch (parser->field) {
    case 1:
      strcpy(parser->transfo->name,s);
      break;
    case 2:
      strcpy(parser->transfo->from_bus,s);
      break;
    case 3:
      strcpy(parser->transfo->to_bus,s);
      break;
    case 4:
      parser->transfo->r = atof(s);
      break;
    case 5:
      parser->transfo->x = atof(s);
      break;
    case 6:
      parser->transfo->b1 = atof(s);
      break;
    case 7:
      parser->transfo->b2 = atof(s);
      break;
    case 8:
      parser->transfo->n = atof(s);
      break;
    case 9:
      parser->transfo->phi = atof(s);
      break;
    case 10:
      parser->transfo->snom = atof(s);
      break;
    case 11:
      parser->transfo->br = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_transfo_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->transfo) {
    LIST_add(parser->transfo_list,parser->transfo,next);
  }
  parser->transfo = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_ltcv_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New ltcv
  if (parser->field == 1) {
    parser->ltcv = (ART_Ltcv*)malloc(sizeof(ART_Ltcv));
    parser->ltcv->next = NULL;    
  }

  // Fields
  if (parser->ltcv) {
    switch (parser->field) {
    case 1:
      strcpy(parser->ltcv->name,s);
      break;
    case 2:
      strcpy(parser->ltcv->con_bus,s);
      break;
    case 3:
      parser->ltcv->nfirst = atof(s);
      break;
    case 4:
      parser->ltcv->nlast = atof(s);
      break;
    case 5:
      parser->ltcv->nbpos = atoi(s);
      break;
    case 6:
      parser->ltcv->tolv = atof(s);
      break;
    case 7:
      parser->ltcv->vdes = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_ltcv_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->transfo) {
    LIST_add(parser->transfo_list,parser->transfo,next);
  }
  parser->transfo = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}
