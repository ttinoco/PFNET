/** @file contingency.c
 *  @brief This file defines the Cont data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/contingency.h>
#include <pfnet/branch.h>
#include <pfnet/gen.h>

// Gen outage
struct Gen_outage {
  Gen* gen;
  struct Gen_outage* next;
};
typedef struct Gen_outage Gen_outage;

// Branch outage
struct Branch_outage {
  Branch* br;
  struct Branch_outage* next;
};
typedef struct Branch_outage Branch_outage;

// Contingency
struct Cont {

  // Generator outages
  Gen_outage* gen_outage;           /**< @brief List of generator outages */

  // Branch outages
  Branch_outage* br_outage;           /**< @brief List of branch outages */
};


void CONT_init(Cont* cont) { 
  if (cont) {
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

void CONT_add_gen_outage(Cont* cont, Gen* gen) {
  if (cont) {
    Gen_outage* go = (Gen_outage*)malloc(sizeof(Gen_outage));
    go->gen = gen;
    LIST_add(cont->gen_outage,go,next);
  }
}

void CONT_add_branch_outage(Cont* cont, Branch* br) {
  if (cont) {
    Branch_outage* bo = (Branch_outage*)malloc(sizeof(Branch_outage));
    bo->br = br;
    LIST_add(cont->br_outage,bo,next);
  }
}

BOOL CONT_has_gen_outage(Cont* cont, Gen* gen) {
  Gen_outage* go;
  if (!cont)
    return FALSE;
  for (go = cont->gen_outage; go != NULL; go = go->next) {
    if (gen == go->gen)
      return TRUE;
  }
  return FALSE; 
}

BOOL CONT_has_branch_outage(Cont* cont, Branch* br) {
  Branch_outage* bo;
  if (!cont)
    return FALSE;
  for (bo = cont->br_outage; bo != NULL; bo = bo->next) {
    if (br == bo->br)
      return TRUE;
  }
  return FALSE; 
}

Cont* CONT_new(void) { 
  Cont* cont = (Cont*)malloc(sizeof(Cont));
  CONT_init(cont);
  return cont;
}

void CONT_show(Cont* cont) {
  Gen_outage* go;
  Branch_outage* bo;
  if (cont) {
    printf("\nGenerator outages\n");
    for (go = cont->gen_outage; go != NULL; go = go->next)
      printf("index %d\n",GEN_get_index(go->gen));
    printf("\nBranch outages\n");
    for (bo = cont->br_outage; bo != NULL; bo = bo->next)
      printf("index %d\n",BRANCH_get_index(bo->br));
  }
}

