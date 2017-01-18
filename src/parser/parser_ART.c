/** @file parser_ART.c
 *  @brief This file defines the ART_Parser and related data structures and their associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/utils.h>
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
  char k_bus[10];     // from/i bus
  char m_bus[10];     // to/j bus`
  REAL r;             // ohms
  REAL x;             // ohms
  REAL wc_half;       // micro siemens
  REAL snom;          // mva
  REAL br;            // breaker
  struct ART_Line* next;
};

struct ART_Transfo {
  char name[22];
  int index;
  char k_bus[10];     // from/i bus
  char m_bus[10];     // to/j bus
  REAL r;             // % on the Vb1,SNOM base
  REAL x;             // % on the Vb1,SNOM base
  REAL b1;            // % on the Vb1,SNOM base
  REAL b2;            // % on the Vb1,SNOM base
  REAL n;             // % on the Vb1,Vb2 base
  REAL phi;           // degrees
  REAL snom;          // mva
  REAL br;            // breaker
  struct ART_Transfo* next;
  UT_hash_handle hh;
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
  char k_bus[10];     // from/i bus number
  char m_bus[10];     // to/j bus number
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
  REAL br;            // breaker
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
  char con_bus[10]; // connected bus
  char mon_bus[10]; // not used
  REAL p;           // mw
  REAL q;           // mvar
  REAL vimp;        // per unit
  REAL snom;        // MVA
  REAL qmin;        // mvar
  REAL qmax;        // mvar
  REAL br;          // breaker
  struct ART_Gener* next;
};

struct ART_Slack {
  char at_bus[10];
  struct ART_Slack* next;
};

struct ART_Vargen {
  char name[22];
  char bus[10];
  REAL p;             // mw
  REAL q;             // mvar
  REAL pmax;          // mw
  REAL pmin;          // mw
  REAL qmax;          // mvar
  REAL qmin;          // mvar
  struct ART_Vargen* next;
};

struct ART_Bat {
  char bus[10];
  REAL p;             // mw
  REAL pmax;          // mw
  REAL pmin;          // mw
  REAL e;             // mwh
  REAL emax;          // mhw
  REAL eta_c;         // unitless
  REAL eta_d;         // unitless
  struct ART_Bat* next;
};

struct ART_Parser {

  // Error
  BOOL error_flag;
  char error_string[ART_PARSER_BUFFER_SIZE];

  // State
  int state;
  int field;
  int record;

  // Options
  int output_level;

  // Base
  REAL base_power; // MVA

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
  ART_Transfo* transfo_hash;

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

  // Vargens
  ART_Vargen* vargen;
  ART_Vargen* vargen_list;

  // Bats
  ART_Bat* bat;
  ART_Bat* bat_list;
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

  // Options
  parser->output_level = 0;

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
  parser->transfo_hash = NULL;

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

  // Vargen
  parser->vargen = NULL;
  parser->vargen_list = NULL;

  // Bats
  parser->bat = NULL;
  parser->bat_list = NULL;

  // Return
  return parser;
}

void ART_PARSER_read(ART_Parser* parser, char* filename) {

  // Local variables
  FILE* file;
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
  int len_vargen_list;
  int len_bat_list;

  if (!parser)
    return;

  // LEVEL 0
  if (parser->output_level <= 0)
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
  LIST_len(ART_Vargen,parser->vargen_list,next,len_vargen_list);
  LIST_len(ART_Bat,parser->bat_list,next,len_bat_list);

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
  printf("vargen list  : %d\n",len_vargen_list);
  printf("bat list     : %d\n",len_bat_list);

  // LEVEL 1
  if (parser->output_level <= 1)
    return;

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
	   line->k_bus,
	   line->m_bus,
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
	   transfo->k_bus,
	   transfo->m_bus,
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
	   trfo->k_bus,
	   trfo->m_bus,
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

  // Debugging VARGEN
  ART_Vargen* vargen;
  for (vargen = parser->vargen_list; vargen != NULL; vargen = vargen->next) {
    printf("Vargen %s %s %.5f %.5f %.5f %.5f %.5f %.5f\n",
	   vargen->name,
	   vargen->bus,
	   vargen->p,
	   vargen->q,
	   vargen->pmin,
	   vargen->pmax,
	   vargen->qmin,
	   vargen->qmax);
  }

  // Debugging BAT
  ART_Bat* bat;
  for (bat = parser->bat_list; bat != NULL; bat = bat->next) {
    printf("Bat %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
	   bat->bus,
	   bat->p,
	   bat->pmin,
	   bat->pmax,
	   bat->e,
	   bat->emax,
	   bat->eta_c,
	   bat->eta_d);
  }
}

void ART_PARSER_load(ART_Parser* parser, Net* net) {

  // Local variables
  int index;
  int num_buses;
  int num_loads;
  int num_shunts;
  int num_gens;
  int num_lines;
  int num_transfo;
  int num_branches;
  int num_vargens;
  int num_bats;
  ART_Bus* art_bus;
  ART_Bus* art_busA;
  ART_Bus* art_busB;
  ART_Slack* art_slack;
  ART_Gener* art_gen;
  ART_Line* art_line;
  ART_Transfo* art_transfo;
  ART_Ltcv* art_ltcv;
  ART_Vargen* art_vargen;
  ART_Bat* art_bat;
  Bus* bus;
  Bus* busA;
  Bus* busB;
  Load* load;
  Shunt* shunt;
  Gen* gen;
  Branch* branch;
  Vargen* vargen;
  Bat* bat;
  REAL r;
  REAL x;
  REAL den;
  REAL g;
  REAL b;
  int num_periods;

  // Check
  if (!parser || !net)
    return;

  // Base power
  NET_set_base_power(net,parser->base_power);

  // Num periods
  num_periods = NET_get_num_periods(net);

  // Buses
  index = 0;
  num_buses = 0;
  LIST_len(ART_Bus,parser->bus_list,next,num_buses);
  NET_set_bus_array(net,BUS_array_new(num_buses,num_periods),num_buses);
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    art_bus->index = index;
    bus = NET_get_bus(net,index);
    BUS_set_number(bus,art_bus->index+1);
    BUS_set_name(bus,art_bus->name);
    NET_bus_hash_number_add(net,bus);
    NET_bus_hash_name_add(net,bus);
    index++;
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
  num_loads = 0;
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->pload != 0. || art_bus->qload != 0. || art_bus->qshunt != 0.)
      num_loads++;
  }
  NET_set_load_array(net,LOAD_array_new(num_loads,num_periods),num_loads);
  index = 0;
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->pload != 0. || art_bus->qload != 0. || art_bus->qshunt != 0.) {
      bus = NET_get_bus(net,art_bus->index);
      load = NET_get_load(net,index);
      BUS_add_load(bus,load);                             // connect load to bus
      LOAD_set_bus(load,bus);                             // connect bus to load
      LOAD_set_P(load,art_bus->pload/parser->base_power,0); // per unit 
      LOAD_set_Q(load,(art_bus->qload-art_bus->qshunt)/parser->base_power,0); // per unit
      LOAD_set_P_min(load,LOAD_get_P(load,0));              // Pmin = P = Pmax
      LOAD_set_P_max(load,LOAD_get_P(load,0));              // Pmin = P = Pmax
      index++;
    }
  }

  // Shunts
  num_shunts = 0;
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->bshunt != 0.)
      num_shunts++;
  }
  NET_set_shunt_array(net,SHUNT_array_new(num_shunts,num_periods),num_shunts);
  index = 0;
  for (art_bus = parser->bus_list; art_bus != NULL; art_bus = art_bus->next) {
    if (art_bus->bshunt != 0.) {
      bus = NET_get_bus(net,art_bus->index);
      shunt = NET_get_shunt(net,index);
      BUS_add_shunt(bus,shunt);                              // connect shunt to bus
      SHUNT_set_bus(shunt,bus);                              // connect bus to shunt
      SHUNT_set_b(shunt,art_bus->bshunt/parser->base_power,0); // per unit
      SHUNT_set_b_max(shunt,SHUNT_get_b(shunt,0));             // per unit
      SHUNT_set_b_min(shunt,SHUNT_get_b(shunt,0));             // per unit
      index++;
    }
  }

  // Generators
  num_gens = 0;
  for (art_gen = parser->gener_list; art_gen != NULL; art_gen = art_gen->next) {
    if (art_gen->br != 0)
      num_gens++;
  }
  NET_set_gen_array(net,GEN_array_new(num_gens,num_periods),num_gens);
  index = 0;
  for (art_gen = parser->gener_list; art_gen != NULL; art_gen = art_gen->next) {
    if (art_gen->br != 0) {
      art_bus = NULL;
      HASH_FIND_STR(parser->bus_hash,art_gen->con_bus,art_bus);
      if (art_bus) {
	bus = NET_get_bus(net,art_bus->index);
	gen = NET_get_gen(net,index);
	BUS_add_gen(bus,gen);                                // connect gen to bus
	GEN_set_bus(gen,bus);                                // connect bus to gen
	GEN_set_P(gen,art_gen->p/parser->base_power,0);        // per unit
	GEN_set_P_max(gen,GEN_INF_P);                        // per unit
	GEN_set_P_min(gen,-GEN_INF_P);                       // per unit
	GEN_set_Q(gen,art_gen->q/parser->base_power,0);        // per unit
	GEN_set_Q_max(gen,art_gen->qmax/parser->base_power); // per unit
	GEN_set_Q_min(gen,art_gen->qmin/parser->base_power); // per unit
	if (art_gen->vimp != 0.) {
	  GEN_set_reg_bus(gen,bus);
	  BUS_add_reg_gen(bus,gen);
	  BUS_set_v_set(bus,art_gen->vimp,0); // p.u.
	  BUS_set_v_mag(bus,art_gen->vimp,0); // p.u.
	}
	else if (BUS_is_slack(bus)) {
	  sprintf(parser->error_string,"invalid imposed voltage of slack generator %s",art_gen->name);
	  parser->error_flag = TRUE;
	}
      }
      else {
	sprintf(parser->error_string,"unable to find generator bus %s",art_gen->con_bus);
	parser->error_flag = TRUE;
      }
      index++;
    }
  }

  // Branches
  num_lines = 0;
  num_transfo = 0;
  num_branches = 0;
  for (art_line = parser->line_list; art_line != NULL; art_line = art_line->next) {
    if (art_line->br != 0)
      num_lines++;
  }
  for (art_transfo = parser->transfo_list; art_transfo != NULL; art_transfo = art_transfo->next) {
    if (art_transfo->br != 0)
      num_transfo++;
  }
  num_branches = num_lines + num_transfo;
  NET_set_branch_array(net,BRANCH_array_new(num_branches,num_periods),num_branches);

  // Lines
  index = 0;
  for (art_line = parser->line_list; art_line != NULL; art_line = art_line->next) {
    if (art_line->br != 0) {

      art_busA = NULL;
      art_busB = NULL;
      HASH_FIND_STR(parser->bus_hash,art_line->k_bus,art_busA);
      HASH_FIND_STR(parser->bus_hash,art_line->m_bus,art_busB);

      if (art_busA && art_busB) {

	busA = NET_get_bus(net,art_busA->index);
	busB = NET_get_bus(net,art_busB->index);
	branch = NET_get_branch(net,index);

	BRANCH_set_type(branch,BRANCH_TYPE_LINE);

	BRANCH_set_bus_k(branch,busA);
	BRANCH_set_bus_m(branch,busB);
	BUS_add_branch_k(busA,branch);
	BUS_add_branch_m(busB,branch);

	r = art_line->r*(parser->base_power*1e6)/pow(art_busA->vnom*1e3,2.); // per unit
	x = art_line->x*(parser->base_power*1e6)/pow(art_busA->vnom*1e3,2.); // per unit

	den = pow(r,2.)+pow(x,2.);
	g = r/den;
	b = -x/den;

	BRANCH_set_g(branch,g);                                // per unit
	BRANCH_set_b(branch,b);                                // per unit

	b = (art_line->wc_half*1e-6)*pow(art_busA->vnom*1e3,2.)/(parser->base_power*1e6); // per unit

	BRANCH_set_b_k(branch,b);            // per unit
	BRANCH_set_b_m(branch,b);              // per unit

      }
      else {
	sprintf(parser->error_string,"unable to find buses of line %s",art_line->name);
	parser->error_flag = TRUE;
      }

      index++;
    }
  }

  // Transfo
  for (art_transfo = parser->transfo_list; art_transfo != NULL; art_transfo = art_transfo->next) {
    if (art_transfo->br != 0) {

      art_transfo->index = index;

      art_busA = NULL;
      art_busB = NULL;
      HASH_FIND_STR(parser->bus_hash,art_transfo->k_bus,art_busA);
      HASH_FIND_STR(parser->bus_hash,art_transfo->m_bus,art_busB);

      if (art_busA && art_busB) {

	busA = NET_get_bus(net,art_busA->index);
	busB = NET_get_bus(net,art_busB->index);
	branch = NET_get_branch(net,index);

	BRANCH_set_type(branch,BRANCH_TYPE_TRAN_FIXED);

	BRANCH_set_bus_k(branch,busB);  // reversed
	BRANCH_set_bus_m(branch,busA);    // reversed
	BUS_add_branch_k(busB,branch);
	BUS_add_branch_m(busA,branch);

	r = (art_transfo->r/100.)*(parser->base_power/art_transfo->snom); // per unit (VB1,SNOM)
	x = (art_transfo->x/100.)*(parser->base_power/art_transfo->snom); // per unit (VB1,SNOM)

	den = pow(r,2.)+pow(x,2.);
	g = r/den;
	b = -x/den;

	BRANCH_set_g(branch,g);                                // per unit
	BRANCH_set_b(branch,b);                                // per unit

	BRANCH_set_b_m(branch,(art_transfo->b1/100.)*(art_transfo->snom/parser->base_power));   // per unit (VB1,SNOM)
	BRANCH_set_b_k(branch,(art_transfo->b2/100.)*(art_transfo->snom/parser->base_power)); // per unit (VB1,SNOM)

	BRANCH_set_ratio(branch,100./art_transfo->n,0);          // units of bus_from_base/bus_to_base
	BRANCH_set_ratio_max(branch,BRANCH_get_ratio(branch,0));
	BRANCH_set_ratio_min(branch,BRANCH_get_ratio(branch,0));

	BRANCH_set_phase(branch,-art_transfo->phi*PI/180.,0);    // radians
	BRANCH_set_phase_max(branch,BRANCH_get_phase(branch,0));
	BRANCH_set_phase_min(branch,BRANCH_get_phase(branch,0));

      }
      else {
	sprintf(parser->error_string,"unable to find buses of transfo %s",art_transfo->name);
	parser->error_flag = TRUE;
      }

      index++;
    }
  }

  // LTC-V
  for (art_ltcv = parser->ltcv_list; art_ltcv != NULL; art_ltcv = art_ltcv->next) {

    art_transfo = NULL;
    HASH_FIND_STR(parser->transfo_hash,art_ltcv->name,art_transfo);
    if (art_transfo) {

      art_bus = NULL;
      HASH_FIND_STR(parser->bus_hash,art_ltcv->con_bus,art_bus);
      if (art_bus) {

	branch = NET_get_branch(net,art_transfo->index);
	bus = NET_get_bus(net,art_bus->index);
	busA = BRANCH_get_bus_k(branch);
	busB = BRANCH_get_bus_m(branch);

	BRANCH_set_type(branch,BRANCH_TYPE_TRAN_TAP_V); // tap changer tap that regulates voltage
	BRANCH_set_reg_bus(branch,bus);                 // branch regulates bus
	BUS_add_reg_tran(bus,branch);               // add regulating transformer to bus

	// Ratio limits
	BRANCH_set_ratio_max(branch,100./art_ltcv->nfirst);
	BRANCH_set_ratio_min(branch,100./art_ltcv->nlast);

	// Voltage limits
	BUS_set_v_set(bus,art_ltcv->vdes,0); // per unit
	BUS_set_v_max(bus,art_ltcv->vdes+art_ltcv->tolv);
	BUS_set_v_min(bus,art_ltcv->vdes-art_ltcv->tolv);

	// tap-voltage sensitivity
	if (busA == bus)
	  BRANCH_set_pos_ratio_v_sens(branch,FALSE); // negative sensitivity
	else if (busB == bus)
	  BRANCH_set_pos_ratio_v_sens(branch,TRUE);  // positive ratio-v sensitivity
	else {
	  sprintf(parser->error_string,"bus regulated by LTC-V transformer %s is not a terminal bus",
		  art_ltcv->name);
	  parser->error_flag = TRUE;
	}
      }
      else {
	sprintf(parser->error_string,"unable to find controlled bus %s of LTC-V transformer %s",
		art_ltcv->con_bus,
		art_ltcv->name);
	parser->error_flag = TRUE;
      }
    }
    else {
      sprintf(parser->error_string,"unable to find LTC-V transformer %s",art_ltcv->name);
      parser->error_flag = TRUE;
    }
  }

  // Vargens
  index = 0;
  LIST_len(ART_Vargen,parser->vargen_list,next,num_vargens);
  NET_set_vargen_array(net,VARGEN_array_new(num_vargens,num_periods),num_vargens);
  for (art_vargen = parser->vargen_list; art_vargen != NULL; art_vargen = art_vargen->next) {
    art_bus = NULL;
    HASH_FIND_STR(parser->bus_hash,art_vargen->bus,art_bus);
    if (art_bus) {
      bus = NET_get_bus(net,art_bus->index);
      vargen = NET_get_vargen(net,index);
      BUS_add_vargen(bus,vargen);                                   // connect vargen to bus
      VARGEN_set_bus(vargen,bus);                                   // connect bus to vargen
      VARGEN_set_name(vargen,art_vargen->name);
      VARGEN_set_P(vargen,art_vargen->p/parser->base_power,0);        // per unit
      VARGEN_set_P_max(vargen,art_vargen->pmax/parser->base_power); // per unit
      VARGEN_set_P_min(vargen,art_vargen->pmin/parser->base_power); // per unit
      VARGEN_set_Q(vargen,art_vargen->q/parser->base_power,0);        // per unit
      VARGEN_set_Q_max(vargen,art_vargen->qmax/parser->base_power); // per unit
      VARGEN_set_Q_min(vargen,art_vargen->qmin/parser->base_power); // per unit
      NET_vargen_hash_name_add(net,vargen);
    }
    else {
      sprintf(parser->error_string,"unable to find var-generator bus %s",art_vargen->bus);
      parser->error_flag = TRUE;
    }
    index++;
  }

  // Batteries
  index = 0;
  LIST_len(ART_Bat,parser->bat_list,next,num_bats);
  NET_set_bat_array(net,BAT_array_new(num_bats,num_periods),num_bats);
  for (art_bat = parser->bat_list; art_bat != NULL; art_bat = art_bat->next) {
    art_bus = NULL;
    HASH_FIND_STR(parser->bus_hash,art_bat->bus,art_bus);
    if (art_bus) {
      bus = NET_get_bus(net,art_bus->index);
      bat = NET_get_bat(net,index);
      BUS_add_bat(bus,bat);                                // connect bat to bus
      BAT_set_bus(bat,bus);                                // connect bus to bat
      BAT_set_P(bat,art_bat->p/parser->base_power,0);        // per unit
      BAT_set_P_max(bat,art_bat->pmax/parser->base_power); // per unit
      BAT_set_P_min(bat,art_bat->pmin/parser->base_power); // per unit
      BAT_set_E(bat,art_bat->e/parser->base_power,0);        // per unit times hour
      BAT_set_E_max(bat,art_bat->emax/parser->base_power); // per unit times hour
      BAT_set_eta_c(bat,art_bat->eta_c);                   // unitless
      BAT_set_eta_d(bat,art_bat->eta_d);                   // unitless
    }
    else {
      sprintf(parser->error_string,"unable to find battery bus %s",art_bat->bus);
      parser->error_flag = TRUE;
    }
    index++;
  }
}

void ART_PARSER_del(ART_Parser* parser) {

  // Check
  if (!parser)
    return;

  // Buses
  while (parser->bus_hash)
    HASH_DEL(parser->bus_hash,parser->bus_hash);
  LIST_map(ART_Bus,parser->bus_list,bus,next,{free(bus);});

  // Line
  LIST_map(ART_Line,parser->line_list,line,next,{free(line);});

  // Transformers
  while (parser->transfo_hash)
    HASH_DEL(parser->transfo_hash,parser->transfo_hash);
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

  // Vargens
  LIST_map(ART_Vargen,parser->vargen_list,vargen,next,{free(vargen);});

  // Bats
  LIST_map(ART_Bat,parser->bat_list,bat,next,{free(bat);});

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
	parser->state = ART_PARSER_STATE_BUS;
      }

      // Line
      else if (strstr(s,ART_LINE_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_LINE;
      }

      // Transformer
      else if (strstr(s,ART_TRANSFO_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_TRANSFO;
      }

      // LTC-V
      else if (strstr(s,ART_LTCV_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_LTCV;
      }

      // TRFOs
      else if (strstr(s,ART_TRFO_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_TRFO;
      }

      // PSHIFTPs
      else if (strstr(s,ART_PSHIFTP_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_PSHIFTP;
      }

      // Generators
      else if (strstr(s,ART_GENER_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_GENER;
      }

      // Slacks
      else if (strstr(s,ART_SLACK_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_SLACK;
      }

      // Vargen
      else if (strstr(s,ART_VARGEN_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_VARGEN;
      }

      // Bat
      else if (strstr(s,ART_BAT_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_BAT;
      }

      // Base
      else if (strstr(s,ART_BASE_TOKEN) != NULL) {
	parser->state = ART_PARSER_STATE_BASE;
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
  case ART_PARSER_STATE_VARGEN:
    ART_PARSER_parse_vargen_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_BAT:
    ART_PARSER_parse_bat_field((char*)s,parser);
    break;
  case ART_PARSER_STATE_BASE:
    ART_PARSER_parse_base_field((char*)s,parser);
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
  case ART_PARSER_STATE_VARGEN:
    ART_PARSER_parse_vargen_record(parser);
    break;
  case ART_PARSER_STATE_BAT:
    ART_PARSER_parse_bat_record(parser);
    break;
  case ART_PARSER_STATE_BASE:
    ART_PARSER_parse_base_record(parser);
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
      strcpy(parser->bus->name,trim(s));
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
    LIST_push(parser->bus_list,parser->bus,next);
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
      strcpy(parser->line->k_bus,s);
      break;
    case 3:
      strcpy(parser->line->m_bus,s);
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

  if (parser->line)
    LIST_push(parser->line_list,parser->line,next);
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
      strcpy(parser->transfo->k_bus,s);
      break;
    case 3:
      strcpy(parser->transfo->m_bus,s);
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
    LIST_push(parser->transfo_list,parser->transfo,next);
    HASH_ADD_STR(parser->transfo_hash,name,parser->transfo);
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

  if (parser->ltcv)
    LIST_push(parser->ltcv_list,parser->ltcv,next);
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
      strcpy(parser->trfo->k_bus,s);
      break;
    case 3:
      strcpy(parser->trfo->m_bus,s);
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

  if (parser->trfo)
    LIST_push(parser->trfo_list,parser->trfo,next);
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

  if (parser->pshiftp)
    LIST_push(parser->pshiftp_list,parser->pshiftp,next);
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

  if (parser->gener)
    LIST_push(parser->gener_list,parser->gener,next);
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

  if (parser->slack)
    LIST_push(parser->slack_list,parser->slack,next);
  parser->slack = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_vargen_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New vargen
  if (parser->field == 1) {
    parser->vargen = (ART_Vargen*)malloc(sizeof(ART_Vargen));
    parser->vargen->next = NULL;
  }

  // Fields
  if (parser->vargen) {
    switch (parser->field) {
    case 1:
      strcpy(parser->vargen->name,s);
      break;
    case 2:
      strcpy(parser->vargen->bus,s);
      break;
    case 3:
      parser->vargen->p = atof(s);
      break;
    case 4:
      parser->vargen->q = atof(s);
      break;
    case 5:
      parser->vargen->pmin = atof(s);
      break;
    case 6:
      parser->vargen->pmax = atof(s);
      break;
    case 7:
      parser->vargen->qmin = atof(s);
      break;
    case 8:
      parser->vargen->qmax = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_vargen_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->vargen)
    LIST_push(parser->vargen_list,parser->vargen,next);
  parser->vargen = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_bat_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // New bat
  if (parser->field == 1) {
    parser->bat = (ART_Bat*)malloc(sizeof(ART_Bat));
    parser->bat->next = NULL;
  }

  // Fields
  if (parser->bat) {
    switch (parser->field) {
    case 1:
      strcpy(parser->bat->bus,s);
      break;
    case 2:
      parser->bat->p = atof(s);
      break;
    case 3:
      parser->bat->pmin = atof(s);
      break;
    case 4:
      parser->bat->pmax = atof(s);
      break;
    case 5:
      parser->bat->e = atof(s);
      break;
    case 6:
      parser->bat->emax = atof(s);
      break;
    case 7:
      parser->bat->eta_c = atof(s);
      break;
    case 8:
      parser->bat->eta_d = atof(s);
      break;
    }
  }
}

void ART_PARSER_parse_bat_record(ART_Parser* parser) {

  if (!parser)
    return;

  if (parser->bat)
    LIST_push(parser->bat_list,parser->bat,next);
  parser->bat = NULL;
  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_parse_base_field(char* s, ART_Parser* parser) {

  if (!parser)
    return;

  // Fields
  switch (parser->field) {
  case 1:
    parser->base_power = atof(s);
    break;
  }
}

void ART_PARSER_parse_base_record(ART_Parser* parser) {

  if (!parser)
    return;

  parser->field = 0;
  parser->record = 0;
  parser->state = ART_PARSER_STATE_INIT;
}

void ART_PARSER_set(ART_Parser* parser, char* key, REAL value) {

  if (!parser)
    return;

  // Output level
  if (strcmp(key,"output_level") == 0)
    parser->output_level = (int)value;
}
