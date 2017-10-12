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
  int bus_index;
  char bus_slack_flag;
  int reg_bus_index;
  BOOL applied;
  struct Gen_outage* next;
};

// Branch outage
struct Branch_outage {
  int br_index;
  int bus_k_index;
  int bus_m_index;
  int reg_bus_index;
  char br_type;
  BOOL applied;
  struct Branch_outage* next;
};

// Types
typedef struct Gen_outage Gen_outage;
typedef struct Branch_outage Branch_outage;

// Contingency
struct Cont {

  // Output
  char output_string[CONT_BUFFER_SIZE]; /**< @brief Output string */

  // Generator outages
  Gen_outage* gen_outage;           /**< @brief List of generator outages */

  // Branch outages
  Branch_outage* br_outage;           /**< @brief List of branch outages */
};


void CONT_apply(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Gen* gen;
  Bus* bus;
  Branch* br;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {

      // Check
      if (go->applied)
	continue;

      // Get data
      gen = NET_get_gen(net,go->gen_index);
      bus = GEN_get_bus(gen);
      reg_bus = GEN_get_reg_bus(gen);

      // Save data
      go->bus_index = BUS_get_index(bus);
      go->bus_slack_flag = BUS_is_slack(bus);
      go->reg_bus_index = BUS_get_index(reg_bus);
      go->applied = TRUE;

      // Outage flag
      GEN_set_outage(gen,TRUE);

      // Connection
      GEN_set_bus(gen,NULL); // disconnect bus from gen
      BUS_del_gen(bus,gen);  // disconnect gen from bus

      // Regulation
      GEN_set_reg_bus(gen,NULL);    // gen does not regulate reg_bus
      BUS_del_reg_gen(reg_bus,gen); // reg_bus is not regulated by gen

      // Slack flag
      if (!BUS_get_gen(bus))
	BUS_set_slack_flag(bus,FALSE);
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {

      // Check
      if (bo->applied)
	continue;

      // Get data
      br = NET_get_branch(net,bo->br_index);
      bus_k = BRANCH_get_bus_k(br);
      bus_m = BRANCH_get_bus_m(br);
      reg_bus = BRANCH_get_reg_bus(br);

      // Save data
      bo->bus_k_index = BUS_get_index(bus_k);
      bo->bus_m_index = BUS_get_index(bus_m);
      bo->reg_bus_index = BUS_get_index(reg_bus);
      bo->br_type = BRANCH_get_type(br);
      bo->applied = TRUE;

      // Outage flag
      BRANCH_set_outage(br,TRUE);

      // Connection
      BRANCH_set_bus_k(br,NULL);  // disconnect bus_k from branch
      BRANCH_set_bus_m(br,NULL);  // disconnect bus_m from branch
      BUS_del_branch_k(bus_k,br); // disconnect branch from bus_k
      BUS_del_branch_m(bus_m,br); // disconnect branch from bus_m

      // Regulation
      BRANCH_set_reg_bus(br,NULL);  // branch does not regulate reg_bus
      BUS_del_reg_tran(reg_bus,br); // reg_bus is not regulated by branch

      // Type
      if (BRANCH_get_type(br) != BRANCH_TYPE_LINE)
	BRANCH_set_type(br,BRANCH_TYPE_TRAN_FIXED);
    }
  }
}

void CONT_clear(Cont* cont, Net* net) {

  // Local variables
  Gen_outage* go;
  Branch_outage* bo;
  Gen* gen;
  Bus* bus;
  Branch* br;
  Bus* bus_k;
  Bus* bus_m;
  Bus* reg_bus;

  if (cont) {

    // Generators
    for (go = cont->gen_outage; go != NULL; go = go->next) {

      // Check
      if (!go->applied)
	continue;

      // Get data
      gen = NET_get_gen(net,go->gen_index);
      bus = NET_get_bus(net,go->bus_index);
      reg_bus = NET_get_bus(net,go->reg_bus_index);
      
      // Outage flag
      GEN_set_outage(gen,FALSE);

      // Connections
      GEN_set_bus(gen,bus); // connect bus to gen
      BUS_add_gen(bus,gen); // connect gen to bus

      // Regulation
      GEN_set_reg_bus(gen,reg_bus); // gen does regulates reg_bus
      BUS_add_reg_gen(reg_bus,gen); // reg_bus is regulated by gen

      // Slack flag
      BUS_set_slack_flag(bus,go->bus_slack_flag);

      // Clear flag
      go->applied = FALSE;
    }

    // Branches
    for (bo = cont->br_outage; bo != NULL; bo = bo->next) {

      // Check
      if (!bo->applied)
	continue;

      // Get data
      br = NET_get_branch(net,bo->br_index);
      bus_k = NET_get_bus(net,bo->bus_k_index);
      bus_m = NET_get_bus(net,bo->bus_m_index);
      reg_bus = NET_get_bus(net,bo->reg_bus_index);

      // Outage flag
      BRANCH_set_outage(br,FALSE);

      // Connection
      BRANCH_set_bus_k(br,bus_k);     // connect bus_k from branch
      BRANCH_set_bus_m(br,bus_m);     // connect bus_m from branch
      BUS_add_branch_k(bus_k,br);     // connect branch from bus_k
      BUS_add_branch_m(bus_m,br);     // connect branch from bus_m

      // Regulation
      BRANCH_set_reg_bus(br,reg_bus); // branch regulates reg_bus
      BUS_add_reg_tran(reg_bus,br);   // reg_bus is regulated by branch

      // Type
      BRANCH_set_type(br,bo->br_type);

      // Clear flag
      bo->applied = FALSE;
    }
  }
}

void CONT_init(Cont* cont) {
  if (cont) {
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
    go->applied = FALSE;
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
    bo->applied = FALSE;
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
  char temp[CONT_BUFFER_SIZE];
  int* indices;
  int num;

  // No contingency
  if (!cont)
    return NULL;

  // Output
  output = (char*)malloc(sizeof(char)*CONT_BUFFER_SIZE*5);
  output_start = output;

  // Write
  JSON_start(output);

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
