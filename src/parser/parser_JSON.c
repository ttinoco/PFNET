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
  //  Net* net;
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
  char* name;
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
    free(file_contents);
    return NULL;
  }

  // Get data
  for (i = 0; i < value->u.object.length; i++) {
    name = value->u.object.values[i].name;
    val = value->u.object.values[i].value;
    if (strcmp(name,"base_power") == 0)
      json_base_power = val;
    else if (strcmp(name,"num_periods") == 0)
      json_num_periods = val;
    else if (strcmp(name,"buses") == 0)
      json_bus_array = val;
    else if (strcmp(name,"branches") == 0)
      json_branch_array = val;
    else if (strcmp(name,"generators") == 0)
      json_gen_array = val;
    else if (strcmp(name,"var_generators") == 0)
      json_vargen_array = val;
    else if (strcmp(name,"shunts") == 0)
      json_shunt_array = val;
    else if (strcmp(name,"loads") == 0)
      json_load_array = val;
    else if (strcmp(name,"batteries") == 0)
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

  // DEBUG
  printf("things going pretty well so far\n");

  // Free
  json_value_free(value);
  free(file_contents);
  
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
