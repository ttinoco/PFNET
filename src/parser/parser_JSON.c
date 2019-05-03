/** @file parser_JSON.c
 *  @brief This file defines the JSON_Parser data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
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
  PARSER_init(p,TRUE);
  return p;
}

void JSON_PARSER_init(Parser* p, BOOL init_params) {
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
  json_value* json_bus_dc_array = NULL;
  json_value* json_branch_dc_array = NULL;
  json_value* json_conv_vsc_array = NULL;
  json_value* json_conv_csc_array = NULL;
  json_value* json_facts_array = NULL;
  json_value* json_red_bus_array = NULL;
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
  int num_dc_buses;
  int num_dc_branches;
  int num_csc_convs;
  int num_vsc_convs;
  int num_facts;
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
    else if (strcmp(key,"dc_buses") == 0)
      json_bus_dc_array = val;
   else if (strcmp(key,"dc_branches") == 0)
     json_branch_dc_array = val;
   else if (strcmp(key,"csc_converters") == 0)
     json_conv_csc_array = val;
   else if (strcmp(key,"vsc_converters") == 0)
     json_conv_vsc_array = val;
   else if (strcmp(key,"facts") == 0)
     json_facts_array = val;
   else if (strcmp(key,"redundant_buses") == 0)
     json_red_bus_array = val;
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
      (!json_bat_array || json_bat_array->type != json_array) ||
      (!json_bus_dc_array || json_bus_dc_array->type != json_array) ||
      (!json_branch_dc_array || json_branch_dc_array->type != json_array) ||
      (!json_conv_vsc_array || json_conv_vsc_array->type != json_array) ||
      (!json_conv_csc_array || json_conv_csc_array->type != json_array) ||
      (!json_facts_array || json_facts_array->type != json_array) ||
      (!json_red_bus_array || json_red_bus_array->type != json_array)) {
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
  num_dc_buses = json_bus_dc_array->u.array.length;
  num_dc_branches = json_branch_dc_array->u.array.length;
  num_vsc_convs = json_conv_vsc_array->u.array.length;
  num_csc_convs = json_conv_csc_array->u.array.length;
  num_facts = json_facts_array->u.array.length;
  
  NET_set_bus_array(net,BUS_array_new(num_buses,num_periods),num_buses);
  NET_set_branch_array(net,BRANCH_array_new(num_branches,num_periods),num_branches);
  NET_set_gen_array(net,GEN_array_new(num_gens,num_periods),num_gens);
  NET_set_vargen_array(net,VARGEN_array_new(num_vargens,num_periods),num_vargens);
  NET_set_shunt_array(net,SHUNT_array_new(num_shunts,num_periods),num_shunts);
  NET_set_load_array(net,LOAD_array_new(num_loads,num_periods),num_loads);
  NET_set_bat_array(net,BAT_array_new(num_bats,num_periods),num_bats);
  NET_set_dc_bus_array(net,BUSDC_array_new(num_dc_buses,num_periods),num_dc_buses);
  NET_set_dc_branch_array(net,BRANCHDC_array_new(num_dc_branches,num_periods),num_dc_branches);
  NET_set_vsc_conv_array(net,CONVVSC_array_new(num_vsc_convs,num_periods),num_vsc_convs);
  NET_set_csc_conv_array(net,CONVCSC_array_new(num_csc_convs,num_periods),num_csc_convs);
  NET_set_facts_array(net,FACTS_array_new(num_facts,num_periods),num_facts);

  // Process arrays
  JSON_PARSER_process_json_bus_array(p,net,json_bus_array);
  JSON_PARSER_process_json_branch_array(p,net,json_branch_array);
  JSON_PARSER_process_json_gen_array(p,net,json_gen_array);
  JSON_PARSER_process_json_vargen_array(p,net,json_vargen_array);
  JSON_PARSER_process_json_shunt_array(p,net,json_shunt_array);
  JSON_PARSER_process_json_load_array(p,net,json_load_array);
  JSON_PARSER_process_json_bat_array(p,net,json_bat_array);
  JSON_PARSER_process_json_bus_dc_array(p,net,json_bus_dc_array);
  JSON_PARSER_process_json_branch_dc_array(p,net,json_branch_dc_array);
  JSON_PARSER_process_json_conv_csc_array(p,net,json_conv_csc_array);
  JSON_PARSER_process_json_conv_vsc_array(p,net,json_conv_vsc_array);
  JSON_PARSER_process_json_facts_array(p,net,json_facts_array);
  JSON_PARSER_process_json_red_bus_array(p,net,json_red_bus_array);

  // Propagate in time
  NET_propagate_data_in_time(net,data_num_periods-1,num_periods);
  
  // Free
  json_value_free(value);
  free(file_contents);
  
  // Return
  return net;
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
  if (json_string)
    fprintf(file,"%s",json_string);

  // Clean up
  free(json_string);
  fclose(file);
}

void JSON_PARSER_free(Parser* p, BOOL del_parser) {
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
      if (strcmp(key,"index") == 0) {
        bus = NET_get_bus(net,val->u.integer);
        break;
      }
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
      if (strcmp(key,"number") == 0) {
        BUS_set_number(bus,val->u.integer);
        NET_bus_hash_number_add(net,bus);
      }

      // oindex
      else if (strcmp(key,"oindex") == 0)
        BUS_set_oindex(bus,val->u.integer);

      // area
      else if (strcmp(key,"area") == 0)
        BUS_set_area(bus,val->u.integer);

      // zone
      else if (strcmp(key,"zone") == 0)
        BUS_set_zone(bus,val->u.integer);

      // name
      else if (strcmp(key,"name") == 0) {
        BUS_set_name(bus,val->u.string.ptr);
        NET_bus_hash_name_add(net,bus);
      }

      // v_base
      else if (strcmp(key,"v_base") == 0)
        BUS_set_v_base(bus,val->u.dbl);
	
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
        BUS_set_slack_flag(bus,val->u.boolean);
      
      // star
      else if (strcmp(key,"star") == 0)
        BUS_set_star_flag(bus,val->u.boolean);
      
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

      // sens_P_balance
      else if (strcmp(key,"sens_P_balance") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_P_balance(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_Q_balance
      else if (strcmp(key,"sens_Q_balance") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_Q_balance(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_mag_u_bound
      else if (strcmp(key,"sens_v_mag_u_bound") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_mag_u_bound(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_mag_l_bound
      else if (strcmp(key,"sens_v_mag_l_bound") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_mag_l_bound(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_ang_u_bound
      else if (strcmp(key,"sens_v_ang_u_bound") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_ang_u_bound(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_ang_l_bound
      else if (strcmp(key,"sens_v_ang_l_bound") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_ang_l_bound(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_set_reg
      else if (strcmp(key,"sens_v_set_reg") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_set_reg(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_reg_by_tran
      else if (strcmp(key,"sens_v_reg_by_tran") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_reg_by_tran(bus,val->u.array.values[k]->u.dbl,k);
      }

      // sens_v_reg_by_shunt
      else if (strcmp(key,"sens_v_reg_by_shunt") == 0) {
        for (k = 0; k < imin(BUS_get_num_periods(bus),val->u.array.length); k++)
          BUS_set_sens_v_reg_by_shunt(bus,val->u.array.values[k]->u.dbl,k);
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
      if (strcmp(key,"index") == 0) {
        branch = NET_get_branch(net,val->u.integer);
        break;
      }
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

      // oindex
      else if (strcmp(key,"oindex") == 0)
        BRANCH_set_oindex(branch,val->u.integer);

      // name
      else if (strcmp(key,"name") == 0)
        BRANCH_set_name(branch,val->u.string.ptr);
      
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

      // num_ratios
      else if (strcmp(key,"num_ratios") == 0) {
        if (val->type == json_integer)
          BRANCH_set_num_ratios(branch,val->u.integer);
      }

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
      
      // in_service
      else if (strcmp(key,"in_service") == 0)
        BRANCH_set_in_service(branch,val->u.boolean);

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
      if (strcmp(key,"index") == 0) {
        gen = NET_get_gen(net,val->u.integer);
        break;
      }
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

      // name
      else if (strcmp(key,"name") == 0)
        GEN_set_name(gen,val->u.string.ptr);

      // in_service
      else if (strcmp(key,"in_service") == 0)
        GEN_set_in_service(gen,val->u.boolean);

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

      // Q_par
      else if (strcmp(key,"Q_par") == 0)
        GEN_set_Q_par(gen,val->u.dbl);

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
      if (strcmp(key,"index") == 0) {
        vargen = NET_get_vargen(net,val->u.integer);
        break;
      }
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

  // Local variables
  Shunt* shunt;
  json_value* json_shunt;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;
  REAL* b_values;
  int num_b;

  // Processs shunt array
  for (i = 0; i < json_shunt_array->u.array.length; i++) {

    // Json shunt
    json_shunt = json_shunt_array->u.array.values[i];
    
    // Check
    if (!json_shunt || json_shunt->type != json_object) {
      PARSER_set_error(p,"Bad json shunt array");
      continue;
    }
    
    // Get shunt
    shunt = NULL;
    for (j = 0; j < json_shunt->u.object.length; j++) {
      key = json_shunt->u.object.values[j].name;
      val = json_shunt->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        shunt = NET_get_shunt(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!shunt) {
      PARSER_set_error(p,"Bad json shunt data");
      continue;
    }

    // Fill
    for (j = 0; j < json_shunt->u.object.length; j++) {
      
      key = json_shunt->u.object.values[j].name;
      val = json_shunt->u.object.values[j].value;
      
      // bus
      if (strcmp(key,"bus") == 0) {
        if (val->type == json_integer)
          SHUNT_set_bus(shunt,NET_get_bus(net,val->u.integer));
      }

      // reg_bus
      else if (strcmp(key,"reg_bus") == 0) {
        if (val->type == json_integer)
          SHUNT_set_reg_bus(shunt,NET_get_bus(net,val->u.integer));
      }

      // type
      else if (strcmp(key,"type") == 0)
        SHUNT_set_type(shunt,val->u.integer);

      // mode
      else if (strcmp(key,"mode") == 0)
        SHUNT_set_mode(shunt,val->u.integer);

      // name
      else if (strcmp(key,"name") == 0)
        SHUNT_set_name(shunt,val->u.string.ptr);
      
      // g
      else if (strcmp(key,"g") == 0)
        SHUNT_set_g(shunt,val->u.dbl); 

      // b
      else if (strcmp(key,"b") == 0) {
        for (k = 0; k < imin(SHUNT_get_num_periods(shunt),val->u.array.length); k++)
          SHUNT_set_b(shunt,val->u.array.values[k]->u.dbl,k);
      }

      // b_max
      else if (strcmp(key,"b_max") == 0)
        SHUNT_set_b_max(shunt,val->u.dbl);

      // b_min
      else if (strcmp(key,"b_min") == 0)
        SHUNT_set_b_min(shunt,val->u.dbl);      

      // b_values
      else if (strcmp(key,"b_values") == 0) {
        num_b = val->u.array.length;
        b_values = (REAL*)malloc(sizeof(REAL)*num_b);
        for (k = 0; k < num_b; k++)
          b_values[k] = val->u.array.values[k]->u.dbl;
        SHUNT_set_b_values(shunt,b_values,num_b); // sets num_b as well
      }      
    }
  }
}

void JSON_PARSER_process_json_load_array(Parser* p, Net* net, json_value* json_load_array) {

  // Local variables
  Load* load;
  json_value* json_load;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs load array
  for (i = 0; i < json_load_array->u.array.length; i++) {

    // Json load
    json_load = json_load_array->u.array.values[i];
    
    // Check
    if (!json_load || json_load->type != json_object) {
      PARSER_set_error(p,"Bad json load array");
      continue;
    }
    
    // Get load
    load = NULL;
    for (j = 0; j < json_load->u.object.length; j++) {
      key = json_load->u.object.values[j].name;
      val = json_load->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        load = NET_get_load(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!load) {
      PARSER_set_error(p,"Bad json load data");
      continue;
    }

    // Fill
    for (j = 0; j < json_load->u.object.length; j++) {
      
      key = json_load->u.object.values[j].name;
      val = json_load->u.object.values[j].value;

      // bus
      if (strcmp(key,"bus") == 0) {
        if (val->type == json_integer)
          LOAD_set_bus(load,NET_get_bus(net,val->u.integer));
      }

      // name
      else if (strcmp(key,"name") == 0)
        LOAD_set_name(load,val->u.string.ptr);

      // in_service
      else if (strcmp(key,"in_service") == 0)
        LOAD_set_in_service(load,val->u.boolean);

      // P
      else if (strcmp(key,"P") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_P(load,val->u.array.values[k]->u.dbl,k);
      }

      // P_max
      else if (strcmp(key,"P_max") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_P_max(load,val->u.array.values[k]->u.dbl,k);
      }


      // P_min
      else if (strcmp(key,"P_min") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_P_min(load,val->u.array.values[k]->u.dbl,k);
      }

      // Q
      else if (strcmp(key,"Q") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_Q(load,val->u.array.values[k]->u.dbl,k);
      }

      // Q max
      else if (strcmp(key,"Q_max") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_Q_max(load,val->u.array.values[k]->u.dbl,k);
      }

      // Q min
      else if (strcmp(key,"Q_min") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_Q_min(load,val->u.array.values[k]->u.dbl,k);
      }

      // comp_cp
      else if (strcmp(key,"comp_cp") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_comp_cp(load,val->u.array.values[k]->u.dbl,k);
      }

      // comp_cq
      else if (strcmp(key,"comp_cq") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_comp_cq(load,val->u.array.values[k]->u.dbl,k);
      }

      // comp_ci
      else if (strcmp(key,"comp_ci") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_comp_ci(load,val->u.array.values[k]->u.dbl,k);
      }

      // comp_cj
      else if (strcmp(key,"comp_cj") == 0) {
        for (k = 0; k < imin(LOAD_get_num_periods(load),val->u.array.length); k++)
          LOAD_set_comp_cj(load,val->u.array.values[k]->u.dbl,k);
      }

      // comp_cg
      else if (strcmp(key,"comp_cg") == 0)
        LOAD_set_comp_cg(load,val->u.dbl);

      // comp_cb
      else if (strcmp(key,"comp_cb") == 0)
        LOAD_set_comp_cb(load,val->u.dbl);

      // target power factor
      else if (strcmp(key,"target_power_factor") == 0)
        LOAD_set_target_power_factor(load,val->u.dbl);

      // util_coeff_Q0
      else if (strcmp(key,"util_coeff_Q0") == 0)
        LOAD_set_util_coeff_Q0(load,val->u.dbl);

      // util_coeff_Q1
      else if (strcmp(key,"util_coeff_Q1") == 0)
        LOAD_set_util_coeff_Q1(load,val->u.dbl);
      
      // util_coeff_Q2
      else if (strcmp(key,"util_coeff_Q2") == 0)
        LOAD_set_util_coeff_Q2(load,val->u.dbl);
    }
  }  
}

void JSON_PARSER_process_json_bat_array(Parser* p, Net* net, json_value* json_bat_array) {

  // Local variables
  Bat* bat;
  json_value* json_bat;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs bat array
  for (i = 0; i < json_bat_array->u.array.length; i++) {

    // Json bat
    json_bat = json_bat_array->u.array.values[i];
    
    // Check
    if (!json_bat || json_bat->type != json_object) {
      PARSER_set_error(p,"Bad json battery array");
      continue;
    }
    
    // Get bat
    bat = NULL;
    for (j = 0; j < json_bat->u.object.length; j++) {
      key = json_bat->u.object.values[j].name;
      val = json_bat->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        bat = NET_get_bat(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!bat) {
      PARSER_set_error(p,"Bad json battery data");
      continue;
    }

    // Fill
    for (j = 0; j < json_bat->u.object.length; j++) {
      
      key = json_bat->u.object.values[j].name;
      val = json_bat->u.object.values[j].value;

      // bus
      if (strcmp(key,"bus") == 0) {
        if (val->type == json_integer)
          BAT_set_bus(bat,NET_get_bus(net,val->u.integer));
      }

      // name
      else if (strcmp(key,"name") == 0)
        BAT_set_name(bat,val->u.string.ptr);

      // P
      else if (strcmp(key,"P") == 0) {
        for (k = 0; k < imin(BAT_get_num_periods(bat),val->u.array.length); k++)
          BAT_set_P(bat,val->u.array.values[k]->u.dbl,k);
      }

      // P_max
      else if (strcmp(key,"P_max") == 0)
        BAT_set_P_max(bat,val->u.dbl);

      // P_min
      else if (strcmp(key,"P_min") == 0)
        BAT_set_P_min(bat,val->u.dbl);

      // eta_c
      else if (strcmp(key,"eta_c") == 0)
        BAT_set_eta_c(bat,val->u.dbl);

      // eta_d
      else if (strcmp(key,"eta_d") == 0)
        BAT_set_eta_d(bat,val->u.dbl);

      // E
      else if (strcmp(key,"E") == 0) {
        for (k = 0; k < imin(BAT_get_num_periods(bat),val->u.array.length); k++)
          BAT_set_E(bat,val->u.array.values[k]->u.dbl,k);
      }

      // E_init
      else if (strcmp(key,"E_init") == 0)
        BAT_set_E_init(bat,val->u.dbl);

      // E_final
      else if (strcmp(key,"E_final") == 0)
        BAT_set_E_final(bat,val->u.dbl);

      // E_max
      else if (strcmp(key,"E_max") == 0)
        BAT_set_E_max(bat,val->u.dbl);
    }
  }
}


void JSON_PARSER_process_json_bus_dc_array(Parser* p, Net* net, json_value* json_bus_dc_array) {

  // Local variables
  BusDC* bus;
  json_value* json_bus;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs bus array
  for (i = 0; i < json_bus_dc_array->u.array.length; i++) {

    // Json bus
    json_bus = json_bus_dc_array->u.array.values[i];
    
    // Check
    if (!json_bus || json_bus->type != json_object) {
      PARSER_set_error(p,"Bad json dc bus array");
      continue;
    }
    
    // Get dc bus
    bus = NULL;
    for (j = 0; j < json_bus->u.object.length; j++) {
      key = json_bus->u.object.values[j].name;
      val = json_bus->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        bus = NET_get_dc_bus(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!bus) {
      PARSER_set_error(p,"Bad json dc bus data");
      continue;
    }

    // Fill
    for (j = 0; j < json_bus->u.object.length; j++) {

      key = json_bus->u.object.values[j].name;
      val = json_bus->u.object.values[j].value;
      
      // number
      if (strcmp(key,"number") == 0) {
        BUSDC_set_number(bus,val->u.integer);
        NET_dc_bus_hash_number_add(net,bus);
      }
      
      // name
      else if (strcmp(key,"name") == 0) {
        BUSDC_set_name(bus,val->u.string.ptr);
        NET_dc_bus_hash_name_add(net,bus);
      }

      // v_base
      else if (strcmp(key,"v_base") == 0)
        BUSDC_set_v_base(bus,val->u.dbl);

      // v
      else if (strcmp(key,"v") == 0) {
        for (k = 0; k < imin(BUSDC_get_num_periods(bus),val->u.array.length); k++)
          BUSDC_set_v(bus,val->u.array.values[k]->u.dbl,k);
      }

      // csc_converters
      else if (strcmp(key,"csc_converters") == 0) {
        for (k = 0; k < val->u.array.length; k++)
          BUSDC_add_csc_conv(bus,NET_get_csc_conv(net,val->u.array.values[k]->u.integer));
      }
      
      // vsc_converters
      else if (strcmp(key,"vsc_converters") == 0) {
        for (k = 0; k < val->u.array.length; k++)
          BUSDC_add_vsc_conv(bus,NET_get_vsc_conv(net,val->u.array.values[k]->u.integer));
      }

      // branches_k
      else if (strcmp(key,"branches_k") == 0) {
        for (k = 0; k < val->u.array.length; k++)
          BUSDC_add_branch_k(bus,NET_get_dc_branch(net,val->u.array.values[k]->u.integer));
      }

      // branches_m
      else if (strcmp(key,"branches_m") == 0) {
        for (k = 0; k < val->u.array.length; k++)
          BUSDC_add_branch_m(bus,NET_get_dc_branch(net,val->u.array.values[k]->u.integer));
      }
    }
  }
}

void JSON_PARSER_process_json_branch_dc_array(Parser* p, Net* net, json_value* json_branch_dc_array) {

  // Local variables
  BranchDC* br;
  json_value* json_br;
  json_value* val;
  char* key;  
  int i;
  int j;

  // Processs branch array
  for (i = 0; i < json_branch_dc_array->u.array.length; i++) {

    // Json br
    json_br = json_branch_dc_array->u.array.values[i];
    
    // Check
    if (!json_br || json_br->type != json_object) {
      PARSER_set_error(p,"Bad json dc branch array");
      continue;
    }
    
    // Get dc branch
    br = NULL;
    for (j = 0; j < json_br->u.object.length; j++) {
      key = json_br->u.object.values[j].name;
      val = json_br->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        br = NET_get_dc_branch(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!br) {
      PARSER_set_error(p,"Bad json dc branch data");
      continue;
    }

    // Fill
    for (j = 0; j < json_br->u.object.length; j++) {

      key = json_br->u.object.values[j].name;
      val = json_br->u.object.values[j].value;

      // name
      if (strcmp(key,"name") == 0)
        BRANCHDC_set_name(br,val->u.string.ptr);

      // bus_k
      else if (strcmp(key,"bus_k") == 0) {
        if (val->type == json_integer)
          BRANCHDC_set_bus_k(br,NET_get_dc_bus(net,val->u.integer));
      }

      // bus_m
      else if (strcmp(key,"bus_m") == 0) {
        if (val->type == json_integer)
          BRANCHDC_set_bus_m(br,NET_get_dc_bus(net,val->u.integer));
      }

      // r
      else if (strcmp(key,"r") == 0)
        BRANCHDC_set_r(br,val->u.dbl);
    }
  }
}

void JSON_PARSER_process_json_conv_csc_array(Parser* p, Net* net, json_value* json_conv_csc_array) {

  // Local variables
  ConvCSC* conv;
  json_value* json_conv;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs conv array
  for (i = 0; i < json_conv_csc_array->u.array.length; i++) {

    // Json conv
    json_conv = json_conv_csc_array->u.array.values[i];
    
    // Check
    if (!json_conv || json_conv->type != json_object) {
      PARSER_set_error(p,"Bad json csc converter array");
      continue;
    }
    
    // Get conv csc
    conv = NULL;
    for (j = 0; j < json_conv->u.object.length; j++) {
      key = json_conv->u.object.values[j].name;
      val = json_conv->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        conv = NET_get_csc_conv(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!conv) {
      PARSER_set_error(p,"Bad json csc converter data");
      continue;
    }

    // Fill
    for (j = 0; j < json_conv->u.object.length; j++) {

      key = json_conv->u.object.values[j].name;
      val = json_conv->u.object.values[j].value;

      // type
      if (strcmp(key,"type") == 0)
        CONVCSC_set_type(conv,val->u.integer);

      // name
      else if (strcmp(key,"name") == 0)
        CONVCSC_set_name(conv,val->u.string.ptr);

      // ac_bus
      else if (strcmp(key,"ac_bus") == 0) {
        if (val->type == json_integer)
          CONVCSC_set_ac_bus(conv,NET_get_bus(net,val->u.integer));
      }

      // dc_bus
      else if (strcmp(key,"dc_bus") == 0) {
        if (val->type == json_integer)
          CONVCSC_set_dc_bus(conv,NET_get_dc_bus(net,val->u.integer));
      }

      // mode_dc
      else if (strcmp(key,"mode_dc") == 0)
        CONVCSC_set_mode_dc(conv,val->u.integer);

      // P_dc_set
      else if (strcmp(key,"P_dc_set") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_P_dc_set(conv,val->u.array.values[k]->u.dbl,k);
      }

      // i_dc_set
      else if (strcmp(key,"i_dc_set") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_i_dc_set(conv,val->u.array.values[k]->u.dbl,k);
      }

      // v_dc_set
      else if (strcmp(key,"v_dc_set") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_v_dc_set(conv,val->u.array.values[k]->u.dbl,k);
      }

      // P
      else if (strcmp(key,"P") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_P(conv,val->u.array.values[k]->u.dbl,k);
      }

      // Q
      else if (strcmp(key,"Q") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_Q(conv,val->u.array.values[k]->u.dbl,k);
      }

      // P_dc
      else if (strcmp(key,"P_dc") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_P_dc(conv,val->u.array.values[k]->u.dbl,k);
      }

      // num_bridges
      if (strcmp(key,"num_bridges") == 0)
        CONVCSC_set_num_bridges(conv,val->u.integer);

      // x_cap
      else if (strcmp(key,"x_cap") == 0)
        CONVCSC_set_x_cap(conv,val->u.dbl);

      // x
      else if (strcmp(key,"x") == 0)
        CONVCSC_set_x(conv,val->u.dbl);

      // r
      else if (strcmp(key,"r") == 0)
        CONVCSC_set_r(conv,val->u.dbl);

      // ratio
      else if (strcmp(key,"ratio") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_ratio(conv,val->u.array.values[k]->u.dbl,k);
      }

      // ratio_max
      else if (strcmp(key,"ratio_max") == 0)
        CONVCSC_set_ratio_max(conv,val->u.dbl);

      // ratio_min
      else if (strcmp(key,"ratio_min") == 0)
        CONVCSC_set_ratio_min(conv,val->u.dbl);

      // angle
      else if (strcmp(key,"angle") == 0) {
        for (k = 0; k < imin(CONVCSC_get_num_periods(conv),val->u.array.length); k++)
          CONVCSC_set_angle(conv,val->u.array.values[k]->u.dbl,k);
      }

      // angle_max
      else if (strcmp(key,"angle_max") == 0)
        CONVCSC_set_angle_max(conv,val->u.dbl);

      // angle_min
      else if (strcmp(key,"angle_min") == 0)
        CONVCSC_set_angle_min(conv,val->u.dbl);

      // v_base_p
      else if (strcmp(key,"v_base_p") == 0)
        CONVCSC_set_v_base_p(conv,val->u.dbl);

      // v_base_s
      else if (strcmp(key,"v_base_s") == 0)
        CONVCSC_set_v_base_s(conv,val->u.dbl);
    }
  }
}

void JSON_PARSER_process_json_conv_vsc_array(Parser* p, Net* net, json_value* json_conv_vsc_array) {

  // Local variables
  ConvVSC* conv;
  json_value* json_conv;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs conv array
  for (i = 0; i < json_conv_vsc_array->u.array.length; i++) {

    // Json conv
    json_conv = json_conv_vsc_array->u.array.values[i];
    
    // Check
    if (!json_conv || json_conv->type != json_object) {
      PARSER_set_error(p,"Bad json vsc converter array");
      continue;
    }
    
    // Get conv vsc
    conv = NULL;
    for (j = 0; j < json_conv->u.object.length; j++) {
      key = json_conv->u.object.values[j].name;
      val = json_conv->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        conv = NET_get_vsc_conv(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!conv) {
      PARSER_set_error(p,"Bad json vsc converter data");
      continue;
    }

    // Fill
    for (j = 0; j < json_conv->u.object.length; j++) {

      key = json_conv->u.object.values[j].name;
      val = json_conv->u.object.values[j].value;

      // name
      if (strcmp(key,"name") == 0)
        CONVVSC_set_name(conv,val->u.string.ptr);
      
      // ac_bus
      else if (strcmp(key,"ac_bus") == 0) {
        if (val->type == json_integer)
          CONVVSC_set_ac_bus(conv,NET_get_bus(net,val->u.integer));
      }

      // dc_bus
      else if (strcmp(key,"dc_bus") == 0) {
        if (val->type == json_integer)
          CONVVSC_set_dc_bus(conv,NET_get_dc_bus(net,val->u.integer));
      }

      // reg_bus
      else if (strcmp(key,"reg_bus") == 0) {
        if (val->type == json_integer)
          CONVVSC_set_reg_bus(conv,NET_get_bus(net,val->u.integer));
      }

      // mode_ac
      else if (strcmp(key,"mode_ac") == 0)
        CONVVSC_set_mode_ac(conv,val->u.integer);

      // mode_dc
      else if (strcmp(key,"mode_dc") == 0)
        CONVVSC_set_mode_dc(conv,val->u.integer);

      // P_dc_set
      else if (strcmp(key,"P_dc_set") == 0) {
        for (k = 0; k < imin(CONVVSC_get_num_periods(conv),val->u.array.length); k++)
          CONVVSC_set_P_dc_set(conv,val->u.array.values[k]->u.dbl,k);
      }

      // v_dc_set
      else if (strcmp(key,"v_dc_set") == 0) {
        for (k = 0; k < imin(CONVVSC_get_num_periods(conv),val->u.array.length); k++)
          CONVVSC_set_v_dc_set(conv,val->u.array.values[k]->u.dbl,k);
      }
      
      // P
      else if (strcmp(key,"P") == 0) {
        for (k = 0; k < imin(CONVVSC_get_num_periods(conv),val->u.array.length); k++)
          CONVVSC_set_P(conv,val->u.array.values[k]->u.dbl,k);
      }

      // Q
      else if (strcmp(key,"Q") == 0) {
        for (k = 0; k < imin(CONVVSC_get_num_periods(conv),val->u.array.length); k++)
          CONVVSC_set_Q(conv,val->u.array.values[k]->u.dbl,k);
      }

      // loss_coeff_A
      else if (strcmp(key,"loss_coeff_A") == 0)
        CONVVSC_set_loss_coeff_A(conv,val->u.dbl);

      // loss_coeff_B
      else if (strcmp(key,"loss_coeff_B") == 0)
        CONVVSC_set_loss_coeff_B(conv,val->u.dbl);

      // P_max
      else if (strcmp(key,"P_max") == 0)
        CONVVSC_set_P_max(conv,val->u.dbl);

      // P_min
      else if (strcmp(key,"P_min") == 0)
        CONVVSC_set_P_min(conv,val->u.dbl);

      // Q_max
      else if (strcmp(key,"Q_max") == 0)
        CONVVSC_set_Q_max(conv,val->u.dbl);

      // Q_min
      else if (strcmp(key,"Q_min") == 0)
        CONVVSC_set_Q_min(conv,val->u.dbl);

      // Q_par
      else if (strcmp(key,"Q_par") == 0)
        CONVVSC_set_Q_par(conv,val->u.dbl);

      // target_power_factor
      else if (strcmp(key,"target_power_factor") == 0)
        CONVVSC_set_target_power_factor(conv,val->u.dbl);

      // P_dc
      else if (strcmp(key,"P_dc") == 0) {
        for (k = 0; k < imin(CONVVSC_get_num_periods(conv),val->u.array.length); k++)
          CONVVSC_set_P_dc(conv,val->u.array.values[k]->u.dbl,k);
      }
    }
  }
}

void JSON_PARSER_process_json_facts_array(Parser* p, Net* net, json_value* json_facts_array) {

  // Local variables
  Facts* facts;
  json_value* json_facts;
  json_value* val;
  char* key;  
  int i;
  int j;
  int k;

  // Processs facts array
  for (i = 0; i < json_facts_array->u.array.length; i++) {

    // Json facts
    json_facts = json_facts_array->u.array.values[i];
    
    // Check
    if (!json_facts || json_facts->type != json_object) {
      PARSER_set_error(p,"Bad json facts array");
      continue;
    }
    
    // Get facts
    facts = NULL;
    for (j = 0; j < json_facts->u.object.length; j++) {
      key = json_facts->u.object.values[j].name;
      val = json_facts->u.object.values[j].value;
      if (strcmp(key,"index") == 0) {
        facts = NET_get_facts(net,val->u.integer);
        break;
      }
    }
    
    // Check
    if (!facts) {
      PARSER_set_error(p,"Bad json facts data");
      continue;
    }

    // Fill
    for (j = 0; j < json_facts->u.object.length; j++) {

      key = json_facts->u.object.values[j].name;
      val = json_facts->u.object.values[j].value;

      // bus_k
      if (strcmp(key,"bus_k") == 0) {
        if (val->type == json_integer)
          FACTS_set_bus_k(facts,NET_get_bus(net,val->u.integer));
      }
      
      // bus_m
      else if (strcmp(key,"bus_m") == 0) {
        if (val->type == json_integer)
          FACTS_set_bus_m(facts,NET_get_bus(net,val->u.integer));
      }

      // reg_bus
      else if (strcmp(key,"reg_bus") == 0) {
        if (val->type == json_integer)
          FACTS_set_reg_bus(facts,NET_get_bus(net,val->u.integer));
      }

      // name
      else if (strcmp(key,"name") == 0)
        FACTS_set_name(facts,val->u.string.ptr);

      // mode_s
      else if (strcmp(key,"mode_s") == 0)
        FACTS_set_mode_s(facts,val->u.integer);

      // P_k
      else if (strcmp(key,"P_k") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_P_k(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_k
      else if (strcmp(key,"Q_k") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_Q_k(facts,val->u.array.values[k]->u.dbl,k);
      }

      // P_m
      else if (strcmp(key,"P_m") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_P_m(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_m
      else if (strcmp(key,"Q_m") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_Q_m(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_sh
      else if (strcmp(key,"Q_sh") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_Q_sh(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_s
      else if (strcmp(key,"Q_s") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_Q_s(facts,val->u.array.values[k]->u.dbl,k);
      }

      // P_dc
      else if (strcmp(key,"P_dc") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_P_dc(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_par
      else if (strcmp(key,"Q_par") == 0)
        FACTS_set_Q_par(facts,val->u.dbl);

      // P_set
      else if (strcmp(key,"P_set") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_P_set(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_set
      else if (strcmp(key,"Q_set") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_Q_set(facts,val->u.array.values[k]->u.dbl,k);
      }

      // Q_max_s
      else if (strcmp(key,"Q_max_s") == 0)
        FACTS_set_Q_max_s(facts,val->u.dbl);

      // Q_min_s
      else if (strcmp(key,"Q_min_s") == 0)
        FACTS_set_Q_min_s(facts,val->u.dbl);

      // Q_max_sh
      else if (strcmp(key,"Q_max_sh") == 0)
        FACTS_set_Q_max_sh(facts,val->u.dbl);

      // Q_min_sh
      else if (strcmp(key,"Q_min_sh") == 0)
        FACTS_set_Q_min_sh(facts,val->u.dbl);

      // i_max_s
      else if (strcmp(key,"i_max_s") == 0)
        FACTS_set_i_max_s(facts,val->u.dbl);

      // i_max_sh
      else if (strcmp(key,"i_max_sh") == 0)
        FACTS_set_i_max_sh(facts,val->u.dbl);

      // P_max_dc
      else if (strcmp(key,"P_max_dc") == 0)
        FACTS_set_P_max_dc(facts,val->u.dbl);

      // v_min_m
      else if (strcmp(key,"v_min_m") == 0)
        FACTS_set_v_min_m(facts,val->u.dbl);

      // v_max_m
      else if (strcmp(key,"v_max_m") == 0)
        FACTS_set_v_max_m(facts,val->u.dbl);

      // v_mag_s
      else if (strcmp(key,"v_mag_s") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_v_mag_s(facts,val->u.array.values[k]->u.dbl,k);
      }

      // v_ang_s
      else if (strcmp(key,"v_ang_s") == 0) {
        for (k = 0; k < imin(FACTS_get_num_periods(facts),val->u.array.length); k++)
          FACTS_set_v_ang_s(facts,val->u.array.values[k]->u.dbl,k);
      }

      // v_max_s
      else if (strcmp(key,"v_max_s") == 0)
        FACTS_set_v_max_s(facts,val->u.dbl);

      // g
      else if (strcmp(key,"g") == 0)
        FACTS_set_g(facts,val->u.dbl);

      // b
      else if (strcmp(key,"b") == 0)
        FACTS_set_b(facts,val->u.dbl);
    }
  }
}

void JSON_PARSER_process_json_red_bus_array(Parser* p, Net* net, json_value* json_red_bus_array) {

  // Local variables
  Bus* bus;
  json_value* json_bus;
  json_value* val;
  char* key;  
  int i;
  int j;

  // Processs red bus array
  for (i = 0; i < json_red_bus_array->u.array.length; i++) {

    // Json bus
    json_bus = json_red_bus_array->u.array.values[i];

    // Check
    if (!json_bus || json_bus->type != json_object) {
      PARSER_set_error(p,"Bad json red bus array");
      continue;
    }
    
    // Create bus
    bus = BUS_new(NET_get_num_periods(net));
   
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

      // alt_number
      if (strcmp(key,"alt_number") == 0)
        BUS_set_alt_number(bus,val->u.integer);

      // alt_name
      else if (strcmp(key,"alt_name") == 0)
        BUS_set_alt_name(bus,val->u.string.ptr);
    }

    // Add bus
    NET_add_red_bus(net,bus);
  }
}
