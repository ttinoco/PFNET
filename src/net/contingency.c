/** @file contingency.c
 *  @brief This file defines the Cont data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/contingency.h>
#include <pfnet/json_macros.h>

// Gen outage
struct Gen_outage {
  int gen_index;
  struct Gen_outage* next;
};

// Branch outage
struct Branch_outage {
  int br_index;
  struct Branch_outage* next;
};

// Types
typedef struct Gen_outage Gen_outage;
typedef struct Branch_outage Branch_outage;

// Contingency
struct Cont {

  // Name
  char name[CONT_BUFFER_SIZE]; /**< @brief Contingency name */

  // Output
  char output_string[CONT_BUFFER_SIZE]; /**< @brief Output string */

  // Generator outages
  Gen_outage* gen_outage;           /**< @brief List of generator outages */

  // Branch outages
  Branch_outage* br_outage;         /**< @brief List of branch outages */
};


void CONT_apply(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Gen* gen;
  Branch* br;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {
      gen = NET_get_gen(net,go->gen_index);
      GEN_set_outage(gen,TRUE);
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
      br = NET_get_branch(net,bo->br_index);
      BRANCH_set_outage(br,TRUE);
    }
  }
}

void CONT_clear(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Gen* gen;
  Branch* br;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {
      gen = NET_get_gen(net,go->gen_index);
      GEN_set_outage(gen,FALSE);
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
      br = NET_get_branch(net,bo->br_index);
      BRANCH_set_outage(br,FALSE);
    }
  }
}

void CONT_init(Cont* cont) {
  if (cont) {
    strcpy(cont->name,"");
    strcpy(cont->output_string,"");
    cont->gen_outage = NULL;
    cont->br_outage = NULL;
  }
}

void CONT_del(Cont* cont) {
  if (cont){
    LIST_map(Gen_outage,cont->gen_outage,go,next,{free(go);});
    LIST_map(Branch_outage,cont->br_outage,bo,next,{free(bo);});
    free(cont);
  }
}

char* CONT_get_name(Cont* cont) {
  if (cont)
    return cont->name;
  else
    return NULL;
}

int CONT_get_num_gen_outages(Cont* cont) {
  int len;
  if (cont) {
    LIST_len(Gen_outage,cont->gen_outage,next,len);
    return len;
  }
  else
    return 0;
}

int CONT_get_num_branch_outages(Cont* cont) {
  int len;
  if (cont) {
    LIST_len(Branch_outage,cont->br_outage,next,len);
    return len;
  }
  else
    return 0;
}

int* CONT_get_branch_outages(Cont* cont) {

  // Local variables
  int* outages;
  Branch_outage* bo;
  int i;
  
  // Allocate
  outages = (int*)malloc(sizeof(int)*CONT_get_num_branch_outages(cont));
  
  // Fill
  i = 0;
  for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
    outages[i] = bo->br_index;
    i++;
  }

  // Return
  return outages;
}

int* CONT_get_gen_outages(Cont* cont) {

  // Local variables
  int* outages;
  Gen_outage* go;
  int i;
  
  // Allocate
  outages = (int*)malloc(sizeof(int)*CONT_get_num_gen_outages(cont));
  
  // Fill
  i = 0;
  for (go = cont->gen_outage; go != NULL; go = go->next) {
    outages[i] = go->gen_index;
    i++;
  }

  // Return
  return outages;
}

void CONT_add_gen_outage(Cont* cont, int gen_index) {
  Gen_outage* go;
  if (cont) {
    for (go = cont->gen_outage; go != NULL; go = go->next) {
      if (go->gen_index == gen_index)
    	return;
    }
    go = (Gen_outage*)malloc(sizeof(Gen_outage));
    go->gen_index = gen_index;
    go->next = NULL;
    LIST_add(Gen_outage,cont->gen_outage,go,next);
  }
}

void CONT_add_branch_outage(Cont* cont, int br_index) {
  Branch_outage* bo;
  if (cont) {
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
      if (bo->br_index == br_index)
	return;
    }
    bo = (Branch_outage*)malloc(sizeof(Branch_outage));
    bo->br_index = br_index;
    bo->next = NULL;
    LIST_add(Branch_outage,cont->br_outage,bo,next);
  }
}

BOOL CONT_has_gen_outage(Cont* cont, int gen_index) {
  Gen_outage* go;
  if (!cont)
    return FALSE;
  for (go = cont->gen_outage; go != NULL; go = go->next) {
    if (gen_index == go->gen_index)
      return TRUE;
  }
  return FALSE;
}

BOOL CONT_has_branch_outage(Cont* cont, int br_index) {
  Branch_outage* bo;
  if (!cont)
    return FALSE;
  for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
    if (br_index == bo->br_index)
      return TRUE;
  }
  return FALSE;
}

Cont* CONT_new(void) {
  Cont* cont = (Cont*)malloc(sizeof(Cont));
  CONT_init(cont);
  return cont;
}

char* CONT_get_show_str(Cont* cont) {

  Gen_outage* go;
  Branch_outage* bo;
  char* out;

  if (!cont)
    return NULL;

  out = cont->output_string;
  strcpy(out,"");
  sprintf(out+strlen(out),"\nName: %s\n",cont->name);
  sprintf(out+strlen(out),"\nGenerator outages\n");
  for (go = cont->gen_outage; go != NULL; go = go->next)
    sprintf(out+strlen(out),"index %d\n",go->gen_index);
  sprintf(out+strlen(out),"\nBranch outages\n");
  for (bo = cont->br_outage; bo != NULL; bo = bo->next)
    sprintf(out+strlen(out),"index %d\n",bo->br_index);

  return out;
}

void CONT_show(Cont* cont) {

  printf("%s",CONT_get_show_str(cont));
}

char* CONT_get_json_string(Cont* cont) {

  // Local vars
  Gen_outage* go;
  Branch_outage* bo;
  char* output;
  char* output_start;
  int buffer_size = CONT_BUFFER_SIZE+JSON_STR_BUFFER_EXTRA;
  char temp[buffer_size];
  int* indices;
  int num;

  // No contingency
  if (!cont)
    return NULL;

  // Output
  output = (char*)malloc(sizeof(char)*buffer_size*5);
  output_start = output;

  // Write
  JSON_start(output);

  // Name
  JSON_str(temp,output,"name",cont->name,FALSE);

  // Gen outages
  num = 0;
  for (go = cont->gen_outage; go != NULL; go = go->next)
    num++;
  indices = (int*)malloc(num*sizeof(int));
  num = 0;
  for (go = cont->gen_outage; go != NULL; go = go->next) {
    indices[num] = go->gen_index;
    num++;
  }
  JSON_array_int(temp,output,"generator_outages",indices,num,FALSE);
  free(indices);

  // Branch outages
  num = 0;
  for (bo = cont->br_outage; bo != NULL; bo = bo->next)
    num++;
  indices = (int*)malloc(num*sizeof(int));
  num = 0;
  for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
    indices[num] = bo->br_index;
    num++;
  }
  JSON_array_int(temp,output,"branch_outages",indices,num,TRUE);
  free(indices);

  // End
  JSON_end(output);

  // Resize
  output = (char*)realloc(output_start,sizeof(char)*(strlen(output_start)+1)); // +1 important!

  // Return
  return output;
}

void CONT_set_name(Cont* cont, char* name) {
  if (cont)
    strncpy(cont->name,name,(size_t)(CONT_BUFFER_SIZE-1));
}
