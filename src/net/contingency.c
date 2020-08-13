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

// Load outage
struct Load_outage {
  int load_index;
  struct Load_outage* next;
};

// Shunt outage
struct Shunt_outage {
  int shunt_index;
  struct Shunt_outage* next;
};

// Types
typedef struct Gen_outage Gen_outage;
typedef struct Branch_outage Branch_outage;
typedef struct Load_outage Load_outage;
typedef struct Shunt_outage Shunt_outage;

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

  // Load outages
  Load_outage* load_outage;         /**< @brief List of load outages */

  // Shunt outages
  Shunt_outage* shunt_outage;       /**< @brief List of shunt outages */
};

void CONT_apply(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Load_outage* lo;
  Shunt_outage* so;
  Gen* gen;
  Branch* br;
  Load* load;
  Shunt* shunt;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {
      gen = NET_get_gen(net,go->gen_index);
      GEN_set_in_service(gen,FALSE);
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
      br = NET_get_branch(net,bo->br_index);
      BRANCH_set_in_service(br,FALSE);
    }

    // Loads
    for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
      load = NET_get_load(net,lo->load_index);
      LOAD_set_in_service(load,FALSE);
    }

    // Shunts
    for (so = cont->shunt_outage; so != NULL; so = so->next) {
      shunt = NET_get_shunt(net,so->shunt_index);
      SHUNT_set_in_service(shunt,FALSE);
    }
  }
}

void CONT_clear(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Load_outage* lo;
  Shunt_outage* so;
  Gen* gen;
  Branch* br;
  Load* load;
  Shunt* shunt;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {
      gen = NET_get_gen(net,go->gen_index);
      GEN_set_in_service(gen,TRUE);
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
      br = NET_get_branch(net,bo->br_index);
      BRANCH_set_in_service(br,TRUE);
    }

    // Loads
    for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
      load = NET_get_load(net,lo->load_index);
      LOAD_set_in_service(load,TRUE);
    }

    // Shunts
    for (so = cont->shunt_outage; so != NULL; so = so->next) {
      shunt = NET_get_shunt(net,so->shunt_index);
      SHUNT_set_in_service(shunt,TRUE);
    }
  }
}

void CONT_init(Cont* cont) {
  if (cont) {
    strcpy(cont->name,"");
    strcpy(cont->output_string,"");
    cont->gen_outage = NULL;
    cont->br_outage = NULL;
    cont->load_outage = NULL;
    cont->shunt_outage = NULL;
  }
}

void CONT_del(Cont* cont) {
  if (cont){
    LIST_map(Gen_outage,cont->gen_outage,go,next,{free(go);});
    LIST_map(Branch_outage,cont->br_outage,bo,next,{free(bo);});
    LIST_map(Load_outage,cont->load_outage,lo,next,{free(lo);});
    LIST_map(Shunt_outage,cont->shunt_outage,so,next,{free(so);});
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

int CONT_get_num_load_outages(Cont* cont) {
  int len;
  if (cont) {
    LIST_len(Load_outage,cont->load_outage,next,len);
    return len;
  }
  else
    return 0;
}

int CONT_get_num_shunt_outages(Cont* cont) {
  int len;
  if (cont) {
    LIST_len(Shunt_outage,cont->shunt_outage,next,len);
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

int* CONT_get_load_outages(Cont* cont) {

  // Local variables
  int* outages;
  Load_outage* lo;
  int i;

  // Allocate
  outages = (int*)malloc(sizeof(int)*CONT_get_num_load_outages(cont));

  // Fill
  i = 0;
  for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
    outages[i] = lo->load_index;
    i++;
  }

  // Return
  return outages;
}

int* CONT_get_shunt_outages(Cont* cont) {

  // Local variables
  int* outages;
  Shunt_outage* so;
  int i;

  // Allocate
  outages = (int*)malloc(sizeof(int)*CONT_get_num_shunt_outages(cont));

  // Fill
  i = 0;
  for (so = cont->shunt_outage; so != NULL; so = so->next) {
    outages[i] = so->shunt_index;
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

void CONT_add_load_outage(Cont* cont, int load_index) {
  Load_outage* lo;
  if (cont) {
    for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
      if (lo->load_index == load_index)
        return;
    }
    lo = (Load_outage*)malloc(sizeof(Load_outage));
    lo->load_index = load_index;
    lo->next = NULL;
    LIST_add(Load_outage,cont->load_outage,lo,next);
  }
}

void CONT_add_shunt_outage(Cont* cont, int shunt_index) {
  Shunt_outage* so;
  if (cont) {
    for (so = cont->shunt_outage; so != NULL; so = so->next) {
      if (so->shunt_index == shunt_index)
        return;
    }
    so = (Shunt_outage*)malloc(sizeof(Shunt_outage));
    so->shunt_index = shunt_index;
    so->next = NULL;
    LIST_add(Shunt_outage,cont->shunt_outage,so,next);
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

BOOL CONT_has_load_outage(Cont* cont, int load_index) {
  Load_outage* lo;
  if (!cont)
    return FALSE;
  for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
    if (load_index == lo->load_index)
      return TRUE;
  }
  return FALSE;
}

BOOL CONT_has_shunt_outage(Cont* cont, int shunt_index) {
  Shunt_outage* so;
  if (!cont)
    return FALSE;
  for (so = cont->shunt_outage; so != NULL; so = so->next) {
    if (shunt_index == so->shunt_index)
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
  Load_outage* lo;
  Shunt_outage* so;
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
  sprintf(out+strlen(out),"\nLoad outages\n");
  for (lo = cont->load_outage; lo != NULL; lo = lo->next)
    sprintf(out+strlen(out),"index %d\n",lo->load_index);
  sprintf(out+strlen(out),"\nShunt outages\n");
  for (so = cont->shunt_outage; so != NULL; so = so->next)
    sprintf(out+strlen(out),"index %d\n",so->shunt_index);

  return out;
}

void CONT_show(Cont* cont) {

  printf("%s",CONT_get_show_str(cont));
}

char* CONT_get_json_string(Cont* cont) {

  // Local vars
  Gen_outage* go;
  Branch_outage* bo;
  Load_outage* lo;
  Shunt_outage* so;
  char* output;
  char* output_start;
  char temp[CONT_JSON_BUFFER_SIZE];
  int* indices;
  int num;

  // No contingency
  if (!cont)
    return NULL;

  // Output
  output = (char*)malloc(sizeof(char)*CONT_BUFFER_SIZE*10);
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
  JSON_array_int(temp,output,"branch_outages",indices,num,FALSE);
  free(indices);

  // Load outages
  num = 0;
  for (lo = cont->load_outage; lo != NULL; lo = lo->next)
    num++;
  indices = (int*)malloc(num*sizeof(int));
  num = 0;
  for (lo = cont->load_outage; lo != NULL; lo = lo->next) {
    indices[num] = lo->load_index;
    num++;
  }
  JSON_array_int(temp,output,"load_outages",indices,num,FALSE);
  free(indices);

  // Shunt outages
  num = 0;
  for (so = cont->shunt_outage; so != NULL; so = so->next)
    num++;
  indices = (int*)malloc(num*sizeof(int));
  num = 0;
  for (so = cont->shunt_outage; so != NULL; so = so->next) {
    indices[num] = so->shunt_index;
    num++;
  }
  JSON_array_int(temp,output,"shunt_outages",indices,num,TRUE);
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
