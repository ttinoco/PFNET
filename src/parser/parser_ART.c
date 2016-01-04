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
  int index;
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

struct ART_Trfo {
  char name[22];
  char from_bus[10];
  char to_bus[10];
  char con_bus[10];
  REAL r;             // % on the Vb1,SNOM base
  REAL x;             // % on the Vb1,SNOM base
  REAL b;             // % on the Vb1,SNOM base
  REAL n;             // % on the Vb1,Vb2 base
  REAL snom;          // mva
  REAL nfirst;        // %           
  REAL nlast;         // %    
  int nbpos; 
  REAL tolv;          // per unit
  REAL vdes;          // per unit
  REAL br;            
  struct ART_Trfo* next;
};

struct ART_Pshiftp {
  char contrfo[22];
  char monbranch[10];
  REAL phafirst;      // degress
  REAL phalast;       // degrees
  int nbpos;  
  int sign;           // 1 or -1 
  REAL pdes;          // mw
  REAL tolp;          // mw
  struct ART_Pshiftp* next;
};

struct ART_Gener {
  char name[22];
  char con_bus[10];
  char mon_bus[10];
  REAL p;           // mw 
  REAL q;           // mvar
  REAL vimp;        // per unit
  REAL snom;        // MVA
  REAL qmin;        // mvar
  REAL qmax;        // mvar
  REAL br;  
  struct ART_Gener* next;
};

struct ART_Slack {
  char at_bus[10];
  struct ART_Slack* next;
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
  int bus_counter;

  // Lines
  ART_Line* line;
  ART_Line* line_list;

  // Transformers
  ART_Transfo* transfo;
  ART_Transfo* transfo_list;

  // LTC-Vs
  ART_Ltcv* ltcv;
  ART_Ltcv* ltcv_list;

  // TRFOs
  ART_Trfo* trfo;
  ART_Trfo* trfo_list;

  // Phase shifters
  ART_Pshiftp* pshiftp;
  ART_Pshiftp* pshiftp_list;

  // Generators
  ART_Gener* gener;
  ART_Gener* gener_list;

  // Slacks
  ART_Slack* slack;
  ART_Slack* slack_list;
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
  parser->bus_counter = 0;

  // Lines
  parser->line = NULL;
  parser->line_list = NULL;

  // Transformers
  parser->transfo = NULL;
  parser->transfo_list = NULL;

  // LTC-Vs
  parser->ltcv = NULL;
  parser->ltcv_list = NULL;

  // TRFOs
  parser->trfo = NULL;
  parser->trfo_list = NULL;

  // Phase shifters
  parser->pshiftp = NULL;
  parser->pshiftp_list = NULL;

  // Generators
  parser->gener = NULL;
  parser->gener_list = NULL;

  // Slacks
  parser->slack = NULL;
  parser->slack_list = NULL;

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
			 '#',
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
  int len_trfo_list;
  int len_pshiftp_list;
  int len_gener_list;
  int len_slack_list;

  if (!parser)
    return;

  // List lengths
  LIST_len(ART_Bus,parser->bus_list,next,len_bus_list);
  LIST_len(ART_Line,parser->line_list,next,len_line_list);
  LIST_len(ART_Transfo,parser->transfo_list,next,len_transfo_list);
  LIST_len(ART_Ltcv,parser->ltcv_list,next,len_ltcv_list);
  LIST_len(ART_Trfo,parser->trfo_list,next,len_trfo_list);
  LIST_len(ART_Pshiftp,parser->pshiftp_list,next,len_pshiftp_list);
  LIST_len(ART_Gener,parser->gener_list,next,len_gener_list);
  LIST_len(ART_Slack,parser->slack_list,next,len_slack_list);
  
  // Show
  printf("\nParsed Data\n");
  printf("base power   : %.2f\n",parser->base_power);
  printf("bus list     : %d\n",len_bus_list);
  printf("line list    : %d\n",len_line_list);
  printf("transfo list : %d\n",len_transfo_list);
  printf("ltc-v list   : %d\n",len_ltcv_list);
  printf("trfo list    : %d\n",len_trfo_list);
  printf("pshiftp list : %d\n",len_pshiftp_list);
  printf("gener list   : %d\n",len_gener_list);
  printf("slack list   : %d\n",len_slack_list);

  // Debugging BUS
  ART_Bus* bus;
  for (bus = parser->bus_list; bus != NULL; bus = bus->next) {
    printf("Bus %s %d %.5f %.5f %.5f %.5f %.5f\n",
	   bus->name,
	   bus->index,
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

  // Debugging TRFO
  ART_Trfo* trfo;
  for (trfo = parser->trfo_list; trfo != NULL; trfo = trfo->next) {
    printf("Trfo %s %s %s %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %d %.5f %.5f %.5f\n",
	   trfo->name,
	   trfo->from_bus,
	   trfo->to_bus,
	   trfo->con_bus,
	   trfo->r,
	   trfo->x,
	   trfo->b,
	   trfo->n,
	   trfo->snom,
	   trfo->nfirst,
	   trfo->nlast,
	   trfo->nbpos,
	   trfo->tolv,
	   trfo->vdes,
	   trfo->br);
  }

  // Debugging PSHIFTP
  ART_Pshiftp* pshiftp;
  for (pshiftp = parser->pshiftp_list; pshiftp != NULL; pshiftp = pshiftp->next) {
    printf("Pshiftp %s %s %.5f %.5f %d %d %.5f %.5f\n",
	   pshiftp->contrfo,
	   pshiftp->monbranch,
	   pshiftp->phafirst,
	   pshiftp->phalast,
	   pshiftp->nbpos,
	   pshiftp->sign,
	   pshiftp->pdes,
	   pshiftp->tolp);
  }

  // Debugging GENER
  ART_Gener* gener;
  for (gener = parser->gener_list; gener != NULL; gener = gener->next) {
    printf("Gener %s %s %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
	   gener->name,
	   gener->con_bus,
	   gener->mon_bus,
	   gener->p,
	   gener->q,
	   gener->vimp,
	   gener->snom,
	   gener->qmin,
	   gener->qmax,
	   gener->br);
  }
  
  // Debugging SLACK
  ART_Slack* slack;
  for (slack = parser->slack_list; slack != NULL; slack = slack->next) {
    printf("Slack %s\n",
	   slack->at_bus);
  }
}

void ART_PARSER_load(ART_Parser* parser, Net* net) {
  
  // Local variables
  int index;
  int num_buses;
  int num_loads;
  ART_Bus* art_bus;
  ART_Slack* art_slack;
  Bus* bus;
  Load* load;

  // Check
  if (!parser || !net)
    return;
  
  // Base power
  NET_set_base_power(net,parser->base_power);

  // Buses
  num_buses = 0;
  LIST_len(ART_Bus,parser->bus_list,next,num_buses);
  NET_set_bus_array(net,BUS_array_new(num_buses),num_buses);
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    bus = NET_get_bus(net,art_bus->index);
    BUS_set_number(bus,art_bus->index+1);
    NET_bus_hash_add(net,bus);
  }

  // Slacks
  for (art_slack = parser->slack_list; art_slack != NULL; art_slack = art_slack->next) {
    art_bus = NULL;
    HASH_FIND_STR(parser->bus_hash,art_slack->at_bus,art_bus);
    if (art_bus) {
      bus = NET_get_bus(net,art_bus->index);
      BUS_set_slack(bus,TRUE);
    }
    else {
      sprintf(parser->error_string,"unable to find slack bus %s",art_slack->at_bus);
      parser->error_flag = TRUE;
    }
  }

  // Loads
  index = 0;
  num_loads = 0;
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->pload != 0. || art_bus->qload != 0.) {
      num_loads++;
    }
  }
  NET_set_load_array(net,LOAD_array_new(num_loads),num_loads);
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->pload != 0. || art_bus->qload != 0.) {
      bus = NET_get_bus(net,art_bus->index);
      load = NET_get_load(net,index);
      BUS_add_load(bus,load);                             // connect load to bus
      LOAD_set_bus(load,bus);                             // connect bus to load
      LOAD_set_P(load,art_bus->pload/parser->base_power); // per unit 
      LOAD_set_Q(load,art_bus->qload/parser->base_power); // per unit
      index++;
    }
  }
  
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

  // TRFOs
  LIST_map(ART_Trfo,parser->trfo_list,trfo,next,{free(trfo);});

  // Phase shifters
  LIST_map(ART_Pshiftp,parser->pshiftp_list,pshiftp,next,{free(pshiftp);});

  // Generators
  LIST_map(ART_Gener,parser->gener_list,gener,next,{free(gener);});

  // Slacks
  LIST_map(ART_Slack,parser->slack_list,slack,next,{free(slack);});

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

      // TRFOs
      else if (strstr(s,ART_TRFO_TOKEN) != NULL) {
	printf("*** TRFO STATE ***\n");
	parser->state = ART_PARSER_STATE_TRFO;
      }

      // PSHIFTPs
      else if (strstr(s,ART_PSHIFTP_TOKEN) != NULL) {
	printf("*** PSHIFTP STATE ***\n");
	parser->state = ART_PARSER_STATE_PSHIFTP;
      }

      // Generators
      else if (strstr(s,ART_GENER_TOKEN) != NULL) {
	printf("*** GENER STATE ***\n");
	parser->state = ART_PARSER_STATE_GENER;
      }

      // Slacks
      else if (strstr(s,ART_SLACK_TOKEN) != NULL) {
	printf("*** SLACK STATE ***\n");
	parser->state = ART_PARSER_STATE_SLACK;
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
  case ART_PARSER_STATE_TRFO:
    ART_PARSER_parse_trfo_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_PSHIFTP:
    ART_PARSER_parse_pshiftp_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_GENER:
    ART_PARSER_parse_gener_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_SLACK:
    ART_PARSER_parse_slack_field((char*)s,parser);
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
  case ART_PARSER_STATE_TRFO:
    ART_PARSER_parse_trfo_record(parser);
    break;
  case ART_PARSER_STATE_PSHIFTP:
    ART_PARSER_parse_pshiftp_record(parser);
    break;
  case ART_PARSER_STATE_GENER:
    ART_PARSER_parse_gener_record(parser);
    break;
  case ART_PARSER_STATE_SLACK:
    ART_PARSER_parse_slack_record(parser);
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
    parser->bus->index = parser->bus_counter;
    parser->bus_counter++;
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

  if (parser->ltcv) {
    LIST_add(parser->ltcv_list,parser->ltcv,next);
  }
  parser->ltcv = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_trfo_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New trfo
  if (parser->field == 1) {
    parser->trfo = (ART_Trfo*)malloc(sizeof(ART_Trfo));
    parser->trfo->next = NULL;
  }

  // Fields
  if (parser->trfo) {
    switch (parser->field) {
    case 1:
      strcpy(parser->trfo->name,s);
      break;
    case 2:
      strcpy(parser->trfo->from_bus,s);
      break;
    case 3:
      strcpy(parser->trfo->to_bus,s);
      break;
    case 4:
      strcpy(parser->trfo->con_bus,s);
      break;
    case 5:
      parser->trfo->r = atof(s);
      break;
    case 6:
      parser->trfo->x = atof(s);
      break;
    case 7:
      parser->trfo->b = atof(s);
      break;
    case 8:
      parser->trfo->n = atof(s);
      break;
    case 9:
      parser->trfo->snom = atof(s);
      break;
    case 10:
      parser->trfo->nfirst = atof(s);
      break;
    case 11:
      parser->trfo->nlast = atof(s);
      break;
    case 12:
      parser->trfo->nbpos = atoi(s);
      break;
    case 13:
      parser->trfo->tolv = atof(s);
      break;
    case 14:
      parser->trfo->vdes = atof(s);
      break;
    case 15:
      parser->trfo->br = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_trfo_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->trfo) {
    LIST_add(parser->trfo_list,parser->trfo,next);
  }
  parser->trfo = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_pshiftp_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New pshiftp
  if (parser->field == 1) {
    parser->pshiftp = (ART_Pshiftp*)malloc(sizeof(ART_Pshiftp));
    parser->pshiftp->next = NULL;    
  }
  
  // Fields
  if (parser->pshiftp) {
    switch (parser->field) {
    case 1:
      strcpy(parser->pshiftp->contrfo,s);
      break;
    case 2:
      strcpy(parser->pshiftp->monbranch,s);
      break;
    case 3:
      parser->pshiftp->phafirst = atof(s);
      break;
    case 4:
      parser->pshiftp->phalast = atof(s);
      break;
    case 5:
      parser->pshiftp->nbpos = atoi(s);
      break;
    case 6:
      parser->pshiftp->sign = atoi(s);
      break;
    case 7:
      parser->pshiftp->pdes = atof(s);
      break;
    case 8:
      parser->pshiftp->tolp = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_pshiftp_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->pshiftp) {
    LIST_add(parser->pshiftp_list,parser->pshiftp,next);
  }
  parser->pshiftp = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_gener_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New gener
  if (parser->field == 1) {
    parser->gener = (ART_Gener*)malloc(sizeof(ART_Gener));
    parser->gener->next = NULL;
  }

  // Fields
  if (parser->gener) {
    switch (parser->field) {
    case 1:
      strcpy(parser->gener->name,s);
      break;
    case 2:
      strcpy(parser->gener->con_bus,s);
      break;
    case 3:
      strcpy(parser->gener->mon_bus,s);
      break;
    case 4:
      parser->gener->p = atof(s);
      break;
    case 5:
      parser->gener->q = atof(s);
      break;
    case 6:
      parser->gener->vimp = atof(s);
      break;
    case 7:
      parser->gener->snom = atof(s);
      break;
    case 8:
      parser->gener->qmin = atof(s);
      break;
    case 9:
      parser->gener->qmax = atof(s);
      break;
    case 10:
      parser->gener->br = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_gener_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->gener) {
    LIST_add(parser->gener_list,parser->gener,next);
  }
  parser->gener = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_slack_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New slack
  if (parser->field == 1) {
    parser->slack = (ART_Slack*)malloc(sizeof(ART_Slack));
    parser->slack->next = NULL;
  }

  // Fields
  if (parser->slack) {
    switch (parser->field) {
    case 1:
      strcpy(parser->slack->at_bus,s);
      break;
    }
  }
}

void ART_PARSER_parse_slack_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->slack) {
    LIST_add(parser->slack_list,parser->slack,next);
  }
  parser->slack = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}
