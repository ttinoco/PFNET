/** @file parser_JSON.c
 *  @brief This file defines the JSON_Parser data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_JSON.h>

Parser* JSON_PARSER_new(void) {
  Parser* p = PARSER_new();
  PARSER_set_func_init(p,&JSON_PARSER_init);
  PARSER_set_func_parse(p,&JSON_PARSER_parse);
  PARSER_set_func_set(p,&JSON_PARSER_set);
  PARSER_set_func_show(p,&JSON_PARSER_show);
  PARSER_set_func_write(p,&JSON_PARSER_write);
  PARSER_set_func_free(p,&JSON_PARSER_free);
  PARSER_init(p);
  return p;
}

void JSON_PARSER_init(Parser* p) {
  // pass
}

Net* JSON_PARSER_parse(Parser* p, char* filename, int num_periods) {

  // Local variables
  Net* net;
  char* ext;
  FILE* file;
  size_t file_size;
  char* file_contents;
  json_char* json_string;
  json_value* value;
  json_value* json_base_power = NULL;
  json_value* json_num_periods = NULL;
  json_value* json_bus_array = NULL;
  json_value* json_branch_array = NULL;
  json_value* json_gen_array = NULL;
  json_value* json_vargen_array = NULL;
  json_value* json_shunt_array = NULL;
  json_value* json_load_array = NULL;
  json_value* json_bat_array = NULL;
  json_value* val;
  char* key;
  REAL data_num_periods;
  int num_buses;
  int num_branches;
  int num_gens;
  int num_vargens;
  int num_shunts;
  int num_loads;
  int num_bats;
  int i;
	
  // Check extension
  ext = strrchr(filename,'.');
  ext = strtolower(ext);
  if (!ext || strcmp(ext+1,"json") != 0) {
    PARSER_set_error(p,"invalid file extension");
    return NULL;
  }

  // Open file
  file = fopen(filename,"rb");
  if (!file) {
    PARSER_set_error(p,"unable to open file");
    return NULL;
  }

  // File size
  fseek(file,0L,SEEK_END);
  file_size = ftell(file);
  rewind(file);

  // Allocate
  file_contents = (char*)malloc(file_size);

  // Read file contents
  if (fread(file_contents,file_size,1,file) != 1 ) {
    PARSER_set_error(p,"Unable to read file contents");
    free(file_contents);
    fclose(file);
    return NULL;
  }
  fclose(file);

  // Parse json string
  json_string = (json_char*)file_contents;
  value = json_parse(json_string,file_size);

  // Checks
  if (value == NULL) {
    PARSER_set_error(p,"Unable to parse json data");
    free(file_contents);
    return NULL;
  }
  if (value->type != json_object) {
    PARSER_set_error(p,"Bad json data");
    json_value_free(value);
    free(file_contents);
    return NULL;
  }

  // Get data
  for (i = 0; i < value->u.object.length; i++) {
    key = value->u.object.values[i].name;
    val = value->u.object.values[i].value;
    if (strcmp(key,"base_power") == 0)
      json_base_power = val;
    else if (strcmp(key,"num_periods") == 0)
      json_num_periods = val;
    else if (strcmp(key,"buses") == 0)
      json_bus_array = val;
    else if (strcmp(key,"branches") == 0)
      json_branch_array = val;
    else if (strcmp(key,"generators") == 0)
      json_gen_array = val;
    else if (strcmp(key,"var_generators") == 0)
      json_vargen_array = val;
    else if (strcmp(key,"shunts") == 0)
      json_shunt_array = val;
    else if (strcmp(key,"loads") == 0)
      json_load_array = val;
    else if (strcmp(key,"batteries") == 0)
      json_bat_array = val;   
  }

  // Check data
  if ((!json_base_power || json_base_power->type != json_double) ||
      (!json_num_periods || json_num_periods->type != json_integer) ||
      (!json_bus_array || json_bus_array->type != json_array) ||
      (!json_branch_array || json_branch_array->type != json_array) ||
      (!json_gen_array || json_gen_array->type != json_array) ||
      (!json_vargen_array || json_vargen_array->type != json_array) ||
      (!json_shunt_array || json_shunt_array->type != json_array) ||
      (!json_load_array || json_load_array->type != json_array) ||
      (!json_bat_array || json_bat_array->type != json_array)) {
    PARSER_set_error(p,"Bad json data");
    json_value_free(value);
    free(file_contents);
  }

  // Num periods
  data_num_periods = json_num_periods->u.integer;
  if (num_periods <= 0)
    num_periods = data_num_periods;

  // Network
  net = NET_new(num_periods);

  // Base power
  NET_set_base_power(net,json_base_power->u.dbl);

  // Set arrays
  num_buses = json_bus_array->u.array.length;
  num_branches = json_branch_array->u.array.length;
  num_gens = json_gen_array->u.array.length;
  num_vargens = json_vargen_array->u.array.length;
  num_shunts = json_shunt_array->u.array.length;
  num_loads = json_load_array->u.array.length;
  num_bats = json_bat_array->u.array.length;  
  NET_set_bus_array(net,BUS_array_new(num_buses,num_periods),num_buses);
  NET_set_branch_array(net,BRANCH_array_new(num_branches,num_periods),num_branches);
  NET_set_gen_array(net,GEN_array_new(num_gens,num_periods),num_gens);
  NET_set_vargen_array(net,VARGEN_array_new(num_vargens,num_periods),num_vargens);
  NET_set_shunt_array(net,SHUNT_array_new(num_shunts,num_periods),num_shunts);
  NET_set_load_array(net,LOAD_array_new(num_loads,num_periods),num_loads);
  NET_set_bat_array(net,BAT_array_new(num_bats,num_periods),num_bats);

  // Process arrays
  JSON_PARSER_process_json_bus_array(p,net,json_bus_array);
  JSON_PARSER_process_json_branch_array(p,net,json_branch_array);
  JSON_PARSER_process_json_gen_array(p,net,json_gen_array);
  JSON_PARSER_process_json_vargen_array(p,net,json_vargen_array);
  JSON_PARSER_process_json_shunt_array(p,net,json_shunt_array);
  JSON_PARSER_process_json_load_array(p,net,json_load_array);
  JSON_PARSER_process_json_bat_array(p,net,json_bat_array);
  
  // Set hash tables
  
  // Free
  json_value_free(value);
  free(file_contents);

  // DEBUG
  printf("a network would be returned here\n");
  
  // Return
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

void JSON_PARSER_process_json_bus_array(Parser* p, Net* net, json_value* json_bus_array) {

  // Local variables
  Bus* bus;
  json_value* json_bus;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs bus array
  for (i = 0; i < json_bus_array->u.array.length; i++) {

    // Json bus
    json_bus = json_bus_array->u.array.values[i];

    // Check
    if (!json_bus || json_bus->type != json_object) {
      PARSER_set_error(p,"Bad json bus array");
      continue;
    }
    
    // Get bus
    bus = NULL;
    for (j = 0; j < json_bus->u.object.length; j++) {
      key = json_bus->u.object.values[j].name;
      val = json_bus->u.object.values[j].value;
      if (strcmp(key,"index") == 0)
	bus = NET_get_bus(net,val->u.integer);
    }
    
    // Check
    if (!bus) {
      PARSER_set_error(p,"Bad json bus data");
      continue;
    }

    // Fill
    for (j = 0; j < json_bus->u.object.length; j++) {
      
      key = json_bus->u.object.values[j].name;
      val = json_bus->u.object.values[j].value;

      // number
      if (strcmp(key,"number") == 0)
	BUS_set_number(bus,val->u.integer);	

      // name
      else if (strcmp(key,"name") == 0)
	BUS_set_name(bus,val->u.string.ptr);
	
      // v_mag
      else if (strcmp(key,"v_mag") == 0) {
	for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
	  BUS_set_v_mag(bus,val->u.array.values[k]->u.dbl,k);
      }

      // v_ang
      else if (strcmp(key,"v_ang") == 0) {
	for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
	  BUS_set_v_ang(bus,val->u.array.values[k]->u.dbl,k);
      }

      // v_set
      else if (strcmp(key,"v_set") == 0) {
	for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
	  BUS_set_v_set(bus,val->u.array.values[k]->u.dbl,k);
      }

      // v_max_reg
      else if (strcmp(key,"v_max_reg") == 0)
	BUS_set_v_max_reg(bus,val->u.dbl);

      // v_min_reg
      else if (strcmp(key,"v_min_reg") == 0)
	BUS_set_v_min_reg(bus,val->u.dbl);

      // v_max_norm
      else if (strcmp(key,"v_max_norm") == 0)
	BUS_set_v_max_norm(bus,val->u.dbl);

      // v_min_norm
      else if (strcmp(key,"v_min_norm") == 0)
	BUS_set_v_min_norm(bus,val->u.dbl);

      // v_max_emer
      else if (strcmp(key,"v_max_emer") == 0)
	BUS_set_v_max_emer(bus,val->u.dbl);

      // v_min_emer
      else if (strcmp(key,"v_min_emer") == 0)
	BUS_set_v_min_emer(bus,val->u.dbl);
      
      // slack
      else if (strcmp(key,"slack") == 0)
	BUS_set_slack(bus,val->u.boolean);
	
      // price
      else if (strcmp(key,"price") == 0) {
	for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
	  BUS_set_price(bus,val->u.array.values[k]->u.dbl,k);
      }

      // gen
      else if (strcmp(key,"generators") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_gen(bus,NET_get_gen(net,val->u.array.values[k]->u.integer));
      }

      // reg_gen
      else if (strcmp(key,"reg_generators") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_reg_gen(bus,NET_get_gen(net,val->u.array.values[k]->u.integer));
      }

      // load
      else if (strcmp(key,"loads") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_load(bus,NET_get_load(net,val->u.array.values[k]->u.integer));
      }

      // shunt
      else if (strcmp(key,"shunts") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_shunt(bus,NET_get_shunt(net,val->u.array.values[k]->u.integer));
      }

      // reg_shunt
      else if (strcmp(key,"reg_shunts") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_reg_shunt(bus,NET_get_shunt(net,val->u.array.values[k]->u.integer));
      }

      // branch_k
      else if (strcmp(key,"branches_k") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_branch_k(bus,NET_get_branch(net,val->u.array.values[k]->u.integer));
      }

      // branch_m
      else if (strcmp(key,"branches_m") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_branch_m(bus,NET_get_branch(net,val->u.array.values[k]->u.integer));
      }

      // reg_tran
      else if (strcmp(key,"reg_transformers") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_reg_tran(bus,NET_get_branch(net,val->u.array.values[k]->u.integer));
      }

      // vargen
      else if (strcmp(key,"var_generators") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_vargen(bus,NET_get_vargen(net,val->u.array.values[k]->u.integer));
      }

      // bat
      else if (strcmp(key,"batteries") == 0) {
	for (k = 0; k < val->u.array.length; k++)
	  BUS_add_bat(bus,NET_get_bat(net,val->u.array.values[k]->u.integer));
      }
    }
  }
}
											
void JSON_PARSER_process_json_branch_array(Parser* p, Net* net, json_value* json_branch_array) {

  // Local variables
  Branch* branch;
  json_value* json_branch;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs branch array
  for (i = 0; i < json_branch_array->u.array.length; i++) {

    // Json branch
    json_branch = json_branch_array->u.array.values[i];
    
    // Check
    if (!json_branch || json_branch->type != json_object) {
      PARSER_set_error(p,"Bad json branch array");
      continue;
    }
    
    // Get branch
    branch = NULL;
    for (j = 0; j < json_branch->u.object.length; j++) {
      key = json_branch->u.object.values[j].name;
      val = json_branch->u.object.values[j].value;
      if (strcmp(key,"index") == 0)
	branch = NET_get_branch(net,val->u.integer);
    }
    
    // Check
    if (!branch) {
      PARSER_set_error(p,"Bad json branch data");
      continue;
    }

    // Fill
    for (j = 0; j < json_branch->u.object.length; j++) {
      
      key = json_branch->u.object.values[j].name;
      val = json_branch->u.object.values[j].value;

      // type
      if (strcmp(key,"type") == 0)
	BRANCH_set_type(branch,val->u.integer);
      
      // bus_k
      else if (strcmp(key,"bus_k") == 0) {
	if (val->type == json_integer)
	  BRANCH_set_bus_k(branch,NET_get_bus(net,val->u.integer));
      }
	    
      // bus_m
      else if (strcmp(key,"bus_m") == 0) {
	if (val->type == json_integer)
	  BRANCH_set_bus_m(branch,NET_get_bus(net,val->u.integer));
      }      
	    
      // reg_bus
      else if (strcmp(key,"reg_bus") == 0) {
	if (val->type == json_integer)
	  BRANCH_set_reg_bus(branch,NET_get_bus(net,val->u.integer));
      }
      
      // g
      else if (strcmp(key,"g") == 0)
	BRANCH_set_g(branch,val->u.dbl);

      // g_k
      else if (strcmp(key,"g_k") == 0)
	BRANCH_set_g_k(branch,val->u.dbl);

      // g_m
      else if (strcmp(key,"g_m") == 0)
	BRANCH_set_g_m(branch,val->u.dbl);
      
      // b
      else if (strcmp(key,"b") == 0)
	BRANCH_set_b(branch,val->u.dbl);

      // b_k
      else if (strcmp(key,"b_k") == 0)
	BRANCH_set_b_k(branch,val->u.dbl);

      // b_m
      else if (strcmp(key,"b_m") == 0)
	BRANCH_set_b_m(branch,val->u.dbl);

      // ratio
      else if (strcmp(key,"ratio") == 0) {
	for (k = 0; k < imin(BRANCH_get_num_periods(branch),val->u.array.length); k++)
	  BRANCH_set_ratio(branch,val->u.array.values[k]->u.dbl,k);
      }
      
      // ratio_max
      else if (strcmp(key,"ratio_max") == 0)
	BRANCH_set_ratio_max(branch,val->u.dbl);

      // ratio_min
      else if (strcmp(key,"ratio_min") == 0)
	BRANCH_set_ratio_min(branch,val->u.dbl);

      // phase
      else if (strcmp(key,"phase") == 0) {
	for (k = 0; k < imin(BRANCH_get_num_periods(branch),val->u.array.length); k++)
	  BRANCH_set_phase(branch,val->u.array.values[k]->u.dbl,k);
      }

      // phase_max
      else if (strcmp(key,"phase_max") == 0)
	BRANCH_set_phase_max(branch,val->u.dbl);

      // phase_min
      else if (strcmp(key,"phase_min") == 0)
	BRANCH_set_phase_min(branch,val->u.dbl);

      // ratingA
      else if (strcmp(key,"ratingA") == 0)
	BRANCH_set_ratingA(branch,val->u.dbl);

      // ratingB
      else if (strcmp(key,"ratingB") == 0)
	BRANCH_set_ratingB(branch,val->u.dbl);

      // ratingC
      else if (strcmp(key,"ratingC") == 0)
	BRANCH_set_ratingC(branch,val->u.dbl);
      
      // outage
      else if (strcmp(key,"outage") == 0)
	BRANCH_set_outage(branch,val->u.boolean);

      // pos_ratio_v_sens
      else if (strcmp(key,"pos_ratio_v_sens") == 0)
	BRANCH_set_pos_ratio_v_sens(branch,val->u.boolean);
    }
  }
}

void JSON_PARSER_process_json_gen_array(Parser* p, Net* net, json_value* json_gen_array) {

  // Local variables
  Gen* gen;
  json_value* json_gen;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs gen array
  for (i = 0; i < json_gen_array->u.array.length; i++) {

    // Json gen
    json_gen = json_gen_array->u.array.values[i];
    
    // Check
    if (!json_gen || json_gen->type != json_object) {
      PARSER_set_error(p,"Bad json gen array");
      continue;
    }
    
    // Get gen
    gen = NULL;
    for (j = 0; j < json_gen->u.object.length; j++) {
      key = json_gen->u.object.values[j].name;
      val = json_gen->u.object.values[j].value;
      if (strcmp(key,"index") == 0)
	gen = NET_get_gen(net,val->u.integer);
    }
    
    // Check
    if (!gen) {
      PARSER_set_error(p,"Bad json gen data");
      continue;
    }

    // Fill
    for (j = 0; j < json_gen->u.object.length; j++) {
      
      key = json_gen->u.object.values[j].name;
      val = json_gen->u.object.values[j].value;

      // bus
      if (strcmp(key,"bus") == 0) {
	if (val->type == json_integer)
	  GEN_set_bus(gen,NET_get_bus(net,val->u.integer));
      }

      // reg_bus
      else if (strcmp(key,"reg_bus") == 0) {
	if (val->type == json_integer)
	  GEN_set_reg_bus(gen,NET_get_bus(net,val->u.integer));
      }

      // outage
      else if (strcmp(key,"outage") == 0)
	GEN_set_outage(gen,val->u.boolean);

      // P
      else if (strcmp(key,"P") == 0) {
	for (k = 0; k < imin(GEN_get_num_periods(gen),val->u.array.length); k++)
	  GEN_set_P(gen,val->u.array.values[k]->u.dbl,k);
      }

      // P_max
      else if (strcmp(key,"P_max") == 0)
	GEN_set_P_max(gen,val->u.dbl);

      // P_min
      else if (strcmp(key,"P_min") == 0)
	GEN_set_P_min(gen,val->u.dbl);

      // dP_max
      else if (strcmp(key,"dP_max") == 0)
	GEN_set_dP_max(gen,val->u.dbl);

      // P_prev
      else if (strcmp(key,"P_prev") == 0)
	GEN_set_P_prev(gen,val->u.dbl);

      // Q
      else if (strcmp(key,"Q") == 0) {
	for (k = 0; k < imin(GEN_get_num_periods(gen),val->u.array.length); k++)
	  GEN_set_Q(gen,val->u.array.values[k]->u.dbl,k);
      }

      // Q_max
      else if (strcmp(key,"Q_max") == 0)
	GEN_set_Q_max(gen,val->u.dbl);

      // Q_min
      else if (strcmp(key,"Q_min") == 0)
	GEN_set_Q_min(gen,val->u.dbl);

      // cost_coeff_Q0
      else if (strcmp(key,"cost_coeff_Q0") == 0)
	GEN_set_cost_coeff_Q0(gen,val->u.dbl);

      // cost_coeff_Q1
      else if (strcmp(key,"cost_coeff_Q1") == 0)
	GEN_set_cost_coeff_Q1(gen,val->u.dbl);
      
      // cost_coeff_Q2
      else if (strcmp(key,"cost_coeff_Q2") == 0)
	GEN_set_cost_coeff_Q2(gen,val->u.dbl);
    }
  }
}

void JSON_PARSER_process_json_vargen_array(Parser* p, Net* net, json_value* json_vargen_array) {

  // Local variables
  Vargen* vargen;
  json_value* json_vargen;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs vargen array
  for (i = 0; i < json_vargen_array->u.array.length; i++) {

    // Json vargen
    json_vargen = json_vargen_array->u.array.values[i];
    
    // Check
    if (!json_vargen || json_vargen->type != json_object) {
      PARSER_set_error(p,"Bad json vargen array");
      continue;
    }
    
    // Get vargen
    vargen = NULL;
    for (j = 0; j < json_vargen->u.object.length; j++) {
      key = json_vargen->u.object.values[j].name;
      val = json_vargen->u.object.values[j].value;
      if (strcmp(key,"index") == 0)
	vargen = NET_get_vargen(net,val->u.integer);
    }
    
    // Check
    if (!vargen) {
      PARSER_set_error(p,"Bad json vargen data");
      continue;
    }

    // Fill
    for (j = 0; j < json_vargen->u.object.length; j++) {
      
      key = json_vargen->u.object.values[j].name;
      val = json_vargen->u.object.values[j].value;

      // bus
      if (strcmp(key,"bus") == 0) {
	if (val->type == json_integer)
	  VARGEN_set_bus(vargen,NET_get_bus(net,val->u.integer));
      }

      // name
      else if (strcmp(key,"name") == 0)
	VARGEN_set_name(vargen,val->u.string.ptr);

      // P
      else if (strcmp(key,"P") == 0) {
	for (k = 0; k < imin(VARGEN_get_num_periods(vargen),val->u.array.length); k++)
	  VARGEN_set_P(vargen,val->u.array.values[k]->u.dbl,k);
      }

      // P_ava
      else if (strcmp(key,"P_ava") == 0) {
	for (k = 0; k < imin(VARGEN_get_num_periods(vargen),val->u.array.length); k++)
	  VARGEN_set_P_ava(vargen,val->u.array.values[k]->u.dbl,k);
      }

      // P_max
      else if (strcmp(key,"P_max") == 0)
	VARGEN_set_P_max(vargen,val->u.dbl);

      // P_min
      else if (strcmp(key,"P_min") == 0)
	VARGEN_set_P_min(vargen,val->u.dbl);

      // P_std
      else if (strcmp(key,"P_std") == 0) {
	for (k = 0; k < imin(VARGEN_get_num_periods(vargen),val->u.array.length); k++)
	  VARGEN_set_P_std(vargen,val->u.array.values[k]->u.dbl,k);
      }

      // Q
      else if (strcmp(key,"Q") == 0) {
	for (k = 0; k < imin(VARGEN_get_num_periods(vargen),val->u.array.length); k++)
	  VARGEN_set_Q(vargen,val->u.array.values[k]->u.dbl,k);
      }

      // Q_max
      else if (strcmp(key,"Q_max") == 0)
	VARGEN_set_Q_max(vargen,val->u.dbl);

      // Q_min
      else if (strcmp(key,"Q_min") == 0)
	VARGEN_set_Q_min(vargen,val->u.dbl);
    }
  }
}

void JSON_PARSER_process_json_shunt_array(Parser* p, Net* net, json_value* json_shunt_array) {

}

void JSON_PARSER_process_json_load_array(Parser* p, Net* net, json_value* json_load_array) {

}

void JSON_PARSER_process_json_bat_array(Parser* p, Net* net, json_value* json_bat_array) {

}
