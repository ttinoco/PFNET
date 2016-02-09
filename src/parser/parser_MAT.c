/** @file parser_MAT.c
 *  @brief This file defines the MAT_Parser and related data structures and their associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_MAT.h>

struct MAT_Bus {
  int number;
  int type;               
  REAL Pd;                // MW
  REAL Qd;                // MVAr
  REAL Gs;                // MW
  REAL Bs;                // MVAr
  int area;               
  REAL Vm;                // per unit
  REAL Va;                // degrees
  REAL basekv;            // kV
  int zone;               
  REAL maxVm;             // per unit
  REAL minVm;             // per unit
  struct MAT_Bus* next;
  UT_hash_handle hh;
};

struct MAT_Gen {
  int bus_number;
  REAL Pg;                // MW
  REAL Qg;                // MVAr
  REAL Qmax;              // MVAr
  REAL Qmin;              // MVAr
  REAL Vg;                // per unit
  REAL mBase;             // MVA
  int status;            
  REAL Pmax;              // MW
  REAL Pmin;              // MW
  struct MAT_Gen* next;

};

struct MAT_Branch {
  int bus_from_number;
  int bus_to_number;
  REAL r;                   // per unit
  REAL x;                   // per unit
  REAL b;                   // per unit
  REAL rateA;               // MVA
  REAL rateB;               // MVA
  REAL rateC;               // MVA
  REAL ratio;               
  REAL angle;               // degrees
  int status;                
  struct MAT_Branch* next;
};

struct MAT_Cost {
  int gen_index;
  REAL Q2;       // $/(hr MW^2)
  REAL Q1;       // $/(hr MW)
  REAL Q0;       // $/(hr)
  struct MAT_Cost* next;
};

struct MAT_Parser {

  // Error
  BOOL error_flag;
  char error_string[MAT_PARSER_BUFFER_SIZE];

  // State
  int state;
  int field;
  int record;
  char token[MAT_PARSER_BUFFER_SIZE];

  // Base
  REAL base_power;
  
  // Buses
  MAT_Bus* bus;
  MAT_Bus* bus_list;
  MAT_Bus* bus_hash;

  // Generators
  MAT_Gen* gen;
  MAT_Gen* gen_list;

  // Branches
  MAT_Branch* branch;
  MAT_Branch* branch_list;

  // Cost
  MAT_Cost* cost;
  MAT_Cost* cost_list;
};

void MAT_PARSER_clear_token(MAT_Parser* parser) {
  int i;
  for (i = 0; i < MAT_PARSER_BUFFER_SIZE; i++)
    parser->token[i] = 0;
}

MAT_Parser* MAT_PARSER_new(void) {
  
  // Allocate
  MAT_Parser* parser = (MAT_Parser*)malloc(sizeof(MAT_Parser));

  // Error
  parser->error_flag = FALSE;
  strcpy(parser->error_string,"");

  // State
  parser->state = MAT_PARSER_STATE_TITLE;
  parser->field = 0;
  parser->record = 0;
  MAT_PARSER_clear_token(parser);

  // Base
  parser->base_power = MAT_PARSER_BASE_POWER;

  // Buses
  parser->bus = NULL;
  parser->bus_list = NULL;
  parser->bus_hash = NULL;
  
  // Generators
  parser->gen = NULL;
  parser->gen_list = NULL;

  // Branches
  parser->branch = NULL;
  parser->branch_list = NULL;

  // Cost
  parser->cost = NULL;
  parser->cost_list = NULL;

  // Return
  return parser;
}

void MAT_PARSER_read(MAT_Parser* parser, char* filename) {
  
  // Local variables
  FILE* file;
  char* line;
  CSV_Parser* csv = CSV_PARSER_new();
  char buffer[MAT_PARSER_BUFFER_SIZE];
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
  while ((bytes_read=fread(buffer,1,MAT_PARSER_BUFFER_SIZE,file)) > 0) {
    if (CSV_PARSER_parse(csv,
			 buffer,
			 bytes_read,
			 feof(file),
			 ',',
			 '\n',
			 0,
			 MAT_PARSER_callback_field,
			 MAT_PARSER_callback_row,
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

void MAT_PARSER_show(MAT_Parser* parser) {

  // Local variables
  int len_bus_list;
  int len_gen_list;
  int len_branch_list;
  int len_cost_list;

  if (!parser)
    return;

  // List lengths
  LIST_len(MAT_Bus,parser->bus_list,next,len_bus_list);
  LIST_len(MAT_Gen,parser->gen_list,next,len_gen_list);
  LIST_len(MAT_Branch,parser->branch_list,next,len_branch_list);
  LIST_len(MAT_Cost,parser->cost_list,next,len_cost_list);

  // Show
  printf("\nParsed Data\n");
  printf("base power : %.2f\n",parser->base_power);
  printf("bus list   : %d\n",len_bus_list);
  printf("bus hash   : %d\n",HASH_COUNT(parser->bus_hash));
  printf("gen list   : %d\n",len_gen_list);
  printf("branch list: %d\n",len_branch_list);
  printf("cost list  : %d\n",len_cost_list);
}

void MAT_PARSER_load(MAT_Parser* parser, Net* net) {

  // Local variables
  int index;
  int num_buses;
  int num_loads;
  int num_shunts;
  int num_gens;
  int num_branches;
  MAT_Bus* mat_bus;
  MAT_Gen* mat_gen;
  MAT_Branch* mat_branch;
  MAT_Cost* mat_cost;
  Bus* bus;
  Bus* busA;
  Bus* busB;
  Load* load;
  Shunt* shunt;
  Gen* gen;
  Branch* branch;
  REAL r;
  REAL x;
  REAL den;
  REAL g;
  REAL b;
  REAL t;
  REAL z;

  if (!parser || !net)
    return;

  // Base
  NET_set_base_power(net,parser->base_power);

  // Buses
  index = 0;
  num_buses = 0;
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    if (mat_bus->type != MAT_BUS_TYPE_IS)
      num_buses++;
  }
  NET_set_bus_array(net,BUS_array_new(num_buses),num_buses); // allocate buses
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    bus = NET_get_bus(net,index);                          // get bus
    BUS_set_number(bus,mat_bus->number); 
    BUS_set_v_mag(bus,mat_bus->Vm);           // per unit
    BUS_set_v_ang(bus,mat_bus->Va*PI/180.); // radians
    if (mat_bus->type == MAT_BUS_TYPE_SL)
      BUS_set_slack(bus,TRUE);
    NET_bus_hash_add(net,bus);
    index++;
  }

  // Loads
  index = 0;
  num_loads = 0;
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    if (mat_bus->type != MAT_BUS_TYPE_IS && (mat_bus->Pd != 0 || mat_bus->Qd != 0))
      num_loads++;
  }
  NET_set_load_array(net,LOAD_array_new(num_loads),num_loads);
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    if (mat_bus->type != MAT_BUS_TYPE_IS && (mat_bus->Pd != 0 || mat_bus->Qd != 0)) {
      bus = BUS_hash_number_find(NET_get_bus_hash(net),mat_bus->number);
      load = NET_get_load(net,index);
      BUS_add_load(bus,load);                          // connect load to bus
      LOAD_set_bus(load,bus);                          // connect bus to load
      LOAD_set_P(load,mat_bus->Pd/parser->base_power); // per unit 
      LOAD_set_Q(load,mat_bus->Qd/parser->base_power); // per unit 
      index++;
    }
  }

  // Shunts
  index = 0;
  num_shunts = 0;
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    if (mat_bus->type != MAT_BUS_TYPE_IS && (mat_bus->Gs != 0 || mat_bus->Bs != 0))
      num_shunts++;
  }
  NET_set_shunt_array(net,SHUNT_array_new(num_shunts),num_shunts);
  for (mat_bus = parser->bus_list; mat_bus != NULL; mat_bus = mat_bus->next) {
    if (mat_bus->type != MAT_BUS_TYPE_IS && (mat_bus->Gs != 0 || mat_bus->Bs != 0)) {
      bus = BUS_hash_number_find(NET_get_bus_hash(net),mat_bus->number);
      shunt = NET_get_shunt(net,index);
      BUS_add_shunt(bus,shunt);                          // connect shunt to bus
      SHUNT_set_bus(shunt,bus);                          // connect bus to shunt
      SHUNT_set_type(shunt,SHUNT_TYPE_FIXED);            // set switchable flag
      SHUNT_set_g(shunt,mat_bus->Gs/parser->base_power); // per unit 
      SHUNT_set_b(shunt,mat_bus->Bs/parser->base_power); // per unit 
      SHUNT_set_b_max(shunt,SHUNT_get_b(shunt));         // per unit
      SHUNT_set_b_min(shunt,SHUNT_get_b(shunt));         // per unit
      index++;
    }
  }

  // Gens
  index = 0;
  num_gens = 0;
  for (mat_gen = parser->gen_list; mat_gen != NULL; mat_gen = mat_gen->next) {
    if (mat_gen->status > 0)
      num_gens++;
  }
  NET_set_gen_array(net,GEN_array_new(num_gens),num_gens);
  for (mat_gen = parser->gen_list; mat_gen != NULL; mat_gen = mat_gen->next) {
    if (mat_gen->status > 0) {
      bus = BUS_hash_number_find(NET_get_bus_hash(net),mat_gen->bus_number);
      gen = NET_get_gen(net,index);
      BUS_add_gen(bus,gen);                                // connect gen to bus
      GEN_set_bus(gen,bus);                                // connect bus to gen
      GEN_set_P(gen,mat_gen->Pg/parser->base_power);       // per unit
      GEN_set_P_max(gen,mat_gen->Pmax/parser->base_power); // per unit
      GEN_set_P_min(gen,mat_gen->Pmin/parser->base_power); // per unit
      GEN_set_Q(gen,mat_gen->Qg/parser->base_power);       // per unit
      GEN_set_Q_max(gen,mat_gen->Qmax/parser->base_power); // per unit
      GEN_set_Q_min(gen,mat_gen->Qmin/parser->base_power); // per unit
      if (BUS_is_slack(bus))  { // generator provides regulation
	GEN_set_regulator(gen,TRUE);        
	GEN_set_reg_bus(gen,bus);           
	BUS_add_reg_gen(bus,gen);
	BUS_set_v_set(bus,mat_gen->Vg); // p.u.
      }
      else if (GEN_get_Q_max(gen) > GEN_get_Q_min(gen)) { // generator provides regulation
	GEN_set_regulator(gen,TRUE);    
	GEN_set_reg_bus(gen,bus);       
	BUS_add_reg_gen(bus,gen);
	BUS_set_v_set(bus,mat_gen->Vg); // p.u.
      }
      else if (GEN_get_Q_max(gen) < GEN_get_Q_min(gen)) {
	GEN_set_Q_max(gen,mat_gen->Qmin/parser->base_power); // per unit
	GEN_set_Q_min(gen,mat_gen->Qmax/parser->base_power); // per unit
      }
      index++;
    }
  }

  // Branches
  index = 0;
  num_branches = 0;
  LIST_len(MAT_Branch,parser->branch_list,next,num_branches);
  NET_set_branch_array(net,BRANCH_array_new(num_branches),num_branches);
  for (mat_branch = parser->branch_list; mat_branch != NULL; mat_branch = mat_branch->next) {
    busA = BUS_hash_number_find(NET_get_bus_hash(net),mat_branch->bus_from_number);
    busB = BUS_hash_number_find(NET_get_bus_hash(net),mat_branch->bus_to_number);
    branch = NET_get_branch(net,index);
    r = mat_branch->r;
    x = mat_branch->x;
    den = pow(r,2.)+pow(x,2.);
    g = r/den;
    b = -x/den;
    if (mat_branch->ratio > 0)
      t = mat_branch->ratio;
    else
      t = 1.;
    z = mat_branch->angle*PI/180.;
    if (t == 1. && z == 0)
      BRANCH_set_type(branch,BRANCH_TYPE_LINE);
    else
      BRANCH_set_type(branch,BRANCH_TYPE_TRAN_FIXED);
    BRANCH_set_bus_from(branch,busA);
    BRANCH_set_bus_to(branch,busB);
    BUS_add_branch_from(busA,branch);
    BUS_add_branch_to(busB,branch);
    BRANCH_set_ratio(branch,1./t);                         // units of bus_from_base/bus_to_base
    BRANCH_set_ratio_max(branch,BRANCH_get_ratio(branch));
    BRANCH_set_ratio_min(branch,BRANCH_get_ratio(branch));
    BRANCH_set_phase(branch,z);                            // radians
    BRANCH_set_phase_max(branch,BRANCH_get_phase(branch));
    BRANCH_set_phase_min(branch,BRANCH_get_phase(branch));			   
    BRANCH_set_g(branch,g);                                // per unit
    BRANCH_set_b(branch,b);                                // per unit
    BRANCH_set_b_from(branch,mat_branch->b/2.);            // per unit
    BRANCH_set_b_to(branch,mat_branch->b/2.);              // per unit
    BRANCH_set_ratingA(branch,mat_branch->rateA/parser->base_power); // p.u.
    BRANCH_set_ratingB(branch,mat_branch->rateB/parser->base_power); // p.u.
    BRANCH_set_ratingC(branch,mat_branch->rateC/parser->base_power); // p.u.
    index++;
  }

  // Costs
  for (mat_cost = parser->cost_list; mat_cost != NULL; mat_cost = mat_cost->next) {
    if (mat_cost->gen_index < NET_get_num_gens(net)) {
      gen = NET_get_gen(net,mat_cost->gen_index);
      GEN_set_cost_coeff_Q2(gen,mat_cost->Q2*100.);
      GEN_set_cost_coeff_Q1(gen,mat_cost->Q1);
      GEN_set_cost_coeff_Q0(gen,mat_cost->Q0/100.);
    }
  }
}

void MAT_PARSER_del(MAT_Parser* parser) {

  if (!parser)
    return;

  // Buses
  while (parser->bus_hash)
    HASH_DEL(parser->bus_hash,parser->bus_hash);
  LIST_map(MAT_Bus,parser->bus_list,bus,next,{free(bus);});

  // Gens
  LIST_map(MAT_Gen,parser->gen_list,gen,next,{free(gen);});

  // Branches
  LIST_map(MAT_Branch,parser->branch_list,branch,next,{free(branch);});

  // Costs
  LIST_map(MAT_Cost,parser->cost_list,cost,next,{free(cost);});

  // Parser
  free(parser);  
}

BOOL MAT_PARSER_has_error(MAT_Parser* parser) {
  if (!parser)
    return FALSE;
  else
    return parser->error_flag;
}

char* MAT_PARSER_get_error_string(MAT_Parser* parser) {
  if (!parser)
    return NULL;
  else
    return parser->error_string;
}

void MAT_PARSER_callback_field(char* s, void* data) {

  // Local variables
  MAT_Parser* parser = (MAT_Parser*)data;

  // Field callback
  switch (parser->state) {
  case MAT_PARSER_STATE_TOKEN:
    MAT_PARSER_parse_token_field((char*)s,parser);
    break;
  case MAT_PARSER_STATE_TITLE:
    MAT_PARSER_parse_title_field((char*)s,parser);
    break;
  case MAT_PARSER_STATE_BUS:
    MAT_PARSER_parse_bus_field((char*)s,parser);
    break;
  case MAT_PARSER_STATE_GEN:
    MAT_PARSER_parse_gen_field((char*)s,parser);
    break;
  case MAT_PARSER_STATE_BRANCH:
    MAT_PARSER_parse_branch_field((char*)s,parser);
    break;
  case MAT_PARSER_STATE_COST:
    MAT_PARSER_parse_cost_field((char*)s,parser);
    break;
  }
    
  // Update field
  parser->field++;

}

void MAT_PARSER_callback_row(void *data) {

  // Local variables
  MAT_Parser* parser = (MAT_Parser*)data;

  // Row callback
  switch (parser->state) {
  case MAT_PARSER_STATE_TOKEN:
    MAT_PARSER_parse_token_row(parser);
    break;
  case MAT_PARSER_STATE_TITLE:
    MAT_PARSER_parse_title_row(parser);
    break;
  case MAT_PARSER_STATE_BUS:
    MAT_PARSER_parse_bus_row(parser);
    break;
  case MAT_PARSER_STATE_GEN:
    MAT_PARSER_parse_gen_row(parser);
    break;
  case MAT_PARSER_STATE_BRANCH:
    MAT_PARSER_parse_branch_row(parser);
    break;
  case MAT_PARSER_STATE_COST:
    MAT_PARSER_parse_cost_row(parser);
    break;
  }
}

void MAT_PARSER_parse_token_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;

  if (strcmp(s,MAT_GEN_TOKEN) == 0) {
    parser->state = MAT_PARSER_STATE_GEN;
    parser->field = 0;
    parser->record = -2;
  }
  else if (strcmp(s,MAT_BRANCH_TOKEN) == 0) {
    parser->state = MAT_PARSER_STATE_BRANCH;
    parser->field = 0;
    parser->record = -2;
  }
  else if (strcmp(s,MAT_COST_TOKEN) == 0) {
    parser->state = MAT_PARSER_STATE_COST;
    parser->field = 0;
    parser->record = -2;
  }
}

void MAT_PARSER_parse_token_row(MAT_Parser* parser) {

  if (!parser)
    return;

  parser->field = 0;
  parser->record = 0;
}

void MAT_PARSER_parse_title_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;

  // Token
  if (strstr(s,MAT_BUS_TOKEN) != NULL) {
    parser->state = MAT_PARSER_STATE_BUS;
    parser->field = 0;
    parser->record = -2;
    return;
  }

  // Base power
  if (parser->field == 0 && parser->record == 1) {
    parser->base_power = atof(s);
  }
}

void MAT_PARSER_parse_title_row(MAT_Parser* parser) {

  if (!parser)
    return;

  parser->field = 0;
  parser->record++;
}

void MAT_PARSER_parse_bus_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;
  
  // Token
  if (strstr(s,MAT_END_TOKEN) != NULL) {
    parser->state = MAT_PARSER_STATE_TOKEN;
    parser->field = 0;
    parser->record = 0;
    return;
  }

  // Labels
  if (parser->record < 0)
    return;

  // New bus
  if (parser->field == 0) {
    parser->bus = (MAT_Bus*)malloc(sizeof(MAT_Bus));
  }

  // Fields
  if (parser->bus) {
    switch (parser->field) {
    case 0: 
      parser->bus->number = atoi(s);
      break;
    case 1:
      parser->bus->type = atoi(s);
      break;
    case 2:
      parser->bus->Pd = atof(s);
      break;
    case 3:
      parser->bus->Qd = atof(s);
      break;
    case 4:
      parser->bus->Gs = atof(s);
      break;
    case 5:
      parser->bus->Bs = atof(s);
      break;
    case 6:
      parser->bus->area = atoi(s);
      break;
    case 7:
      parser->bus->Vm = atof(s);
      break;
    case 8:
      parser->bus->Va = atof(s);
      break;
    case 9:
      parser->bus->basekv = atof(s);
      break;
    case 10:
      parser->bus->zone = atoi(s);
      break;
    case 11:
      parser->bus->maxVm = atof(s);
      break;
    case 12:
      parser->bus->minVm = atof(s);
      break;
    }
  }
}

void MAT_PARSER_parse_bus_row(MAT_Parser* parser) {

  if (!parser)
    return;

  if (parser->bus && parser->record >= 0) {
    LIST_add(parser->bus_list,parser->bus,next);
    HASH_ADD_INT(parser->bus_hash,number,parser->bus);
  }
  parser->bus = NULL;
  parser->field = 0;
  parser->record++;
}

void MAT_PARSER_parse_gen_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;

  // Token
  if (strstr(s,MAT_END_TOKEN) != NULL) {
    parser->state = MAT_PARSER_STATE_TOKEN;
    parser->field = 0;
    parser->record = 0;
    return;
  }

  // Labels
  if (parser->record < 0)
    return;

  // New gen
  if (parser->field == 0) {
    parser->gen = (MAT_Gen*)malloc(sizeof(MAT_Gen));
  }

  // Fields
  if (parser->gen) {
    switch (parser->field) {
    case 0: 
      parser->gen->bus_number = atoi(s);
      break;
    case 1: 
      parser->gen->Pg = atof(s);
      break;
    case 2: 
      parser->gen->Qg = atof(s);
      break;
    case 3: 
      parser->gen->Qmax = atof(s);
      break;
    case 4: 
      parser->gen->Qmin = atof(s);
      break;
    case 5: 
      parser->gen->Vg = atof(s);
      break;
    case 6: 
      parser->gen->mBase = atof(s);
      break;
    case 7: 
      parser->gen->status = atoi(s);
      break;
    case 8: 
      parser->gen->Pmax = atof(s);
      break;
    case 9: 
      parser->gen->Pmin = atof(s);
      break;
    }
  }
}

void MAT_PARSER_parse_gen_row(MAT_Parser* parser) {

  if (!parser)
    return;
  
  if (parser->gen && parser->record >= 0)
    LIST_add(parser->gen_list,parser->gen,next);
  parser->gen = NULL;
  parser->field = 0;
  parser->record++;
}

void MAT_PARSER_parse_branch_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;

  // Token
  if (strstr(s,MAT_END_TOKEN) != NULL) {
    parser->state = MAT_PARSER_STATE_TOKEN;
    parser->field = 0;
    parser->record = 0;
    return;
  }

  // Labels
  if (parser->record < 0)
    return;

  // New branch
  if (parser->field == 0) {
    parser->branch = (MAT_Branch*)malloc(sizeof(MAT_Branch));
  }

  // Fields
  if (parser->branch) {
    switch (parser->field) {
    case 0: 
      parser->branch->bus_from_number = atoi(s);
      break;
    case 1: 
      parser->branch->bus_to_number = atoi(s);
      break;
    case 2: 
      parser->branch->r = atof(s);
      break;
    case 3: 
      parser->branch->x = atof(s);
      break;
    case 4: 
      parser->branch->b = atof(s);
      break;
    case 5: 
      parser->branch->rateA = atof(s);
      break;
    case 6: 
      parser->branch->rateB = atof(s);
      break;
    case 7: 
      parser->branch->rateC = atof(s);
      break;
    case 8: 
      parser->branch->ratio = atof(s);
      break;
    case 9: 
      parser->branch->angle = atof(s);
      break;
    }
  }
}

void MAT_PARSER_parse_branch_row(MAT_Parser* parser) {

  if (!parser)
    return;

  if (parser->branch && parser->record >= 0)
    LIST_add(parser->branch_list,parser->branch,next);
  parser->branch = NULL;
  parser->field = 0;
  parser->record++;  
}

void MAT_PARSER_parse_cost_field(char* s, MAT_Parser* parser) {

  if (!parser)
    return;

  // Token
  if (strstr(s,MAT_END_TOKEN) != NULL) {
    parser->state = MAT_PARSER_STATE_TOKEN;
    parser->field = 0;
    parser->record = 0;
    return;
  }

  // Labels
  if (parser->record < 0)
    return;

  // New cost
  if (parser->field == 0) {
    parser->cost = (MAT_Cost*)malloc(sizeof(MAT_Cost));
  }

  // Fields
  if (parser->cost) {
    switch (parser->field) {
    case 0: 
      parser->cost->gen_index = atoi(s);
      break;
    case 1: 
      parser->cost->Q2 = atof(s);
      break;
    case 2: 
      parser->cost->Q1 = atof(s);
      break;
    case 3: 
      parser->cost->Q0 = atof(s);
      break;
    }
  }
}

void MAT_PARSER_parse_cost_row(MAT_Parser* parser) {

  if (!parser)
    return;
  
  if (parser->cost && parser->record >= 0)
    LIST_add(parser->cost_list,parser->cost,next);
  parser->cost = NULL;
  parser->field = 0;
  parser->record++;
}
