/** @file func.c
 *  @brief This file defines the Func data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func.h>
#include <pfnet/func_REG_VMAG.h>
#include <pfnet/func_REG_VANG.h>
#include <pfnet/func_REG_PQ.h>
#include <pfnet/func_REG_RATIO.h>
#include <pfnet/func_REG_PHASE.h>
#include <pfnet/func_REG_SUSC.h>
#include <pfnet/func_GEN_COST.h>
#include <pfnet/func_SP_CONTROLS.h>
#include <pfnet/func_SLIM_VMAG.h>

struct Func {
  
  // Error
  BOOL error_flag;                     /**< @brief Error flag */
  char error_string[FUNC_BUFFER_SIZE]; /**< @brief Error string */

  // Network
  Net* net;    /**< @brief Power network */
  
  // Type
  int type;    /**< @brief Function type */

  // Weight
  REAL weight; /**< @brief Function weight for forming an objective function */

  // Value
  REAL phi;    /**< @brief Function value */

  // Gradient
  Vec* gphi;   /**< @brief Function gradient vector */

  // Hessian
  Mat* Hphi;   /**< @brief Function Hessian matrix (only lower triangular part stored) */

  // Utils
  char* bus_counted;    /**< @brief Flags for processing buses */
  int bus_counted_size; /**< @brief Size of array of flags for processing buses */
  int Hcounter;         /**< @brief Counter of number of nonzero elements of the Hessian */
  int branch_counter;   /**< @brief Counter of network branches */

  // Functions
  void (*func_init)(Func* f); /**< @brief Initialization function */
  void (*func_count_branch)(Func* f, Branch *branch); /**< @brief Function for countinng nonzero entries */
  void (*func_allocate)(Func* f); /**< @brief Function for allocating required arrays */
  void (*func_clear)(Func* f); /**< @brief Function for clearing flags, counters, and function values */
  void (*func_analyze_branch)(Func* f, Branch *branch); /**< @brief Function for analyzing sparsity pattern */
  void (*func_eval_branch)(Func* f, Branch* branch, Vec* var_values); /**< @brief Function for evaluating function */
  void (*func_free)(Func* f); /**< @brief Function for de-allocating any data used */

  // List
  Func* next; /**< @brief List of functions for forming objective function */
};

void FUNC_clear_bus_counted(Func* f) {
  int i;
  if (f) {
    if (f->bus_counted) {
      for (i = 0; i < f->bus_counted_size; i++)
	f->bus_counted[i] = 0;
    }
  } 
}

void FUNC_del(Func* f) {
  if (f) {

    // Mat and vec
    VEC_del(f->gphi);
    MAT_del(f->Hphi);

    // Utils
    if (f->bus_counted)
      free(f->bus_counted);

    // Data
    if (f->func_free)
      (*(f->func_free))(f);

    free(f);
  }
}

int FUNC_get_type(Func* f) {
  if (f)
    return f->type;
  else
    return FUNC_TYPE_UNKNOWN;
}

char* FUNC_get_type_str(Func* f) {
  if (f)
    switch (f->type) {
    case FUNC_TYPE_REG_VMAG:
      return FUNC_TYPE_REG_VMAG_STR;
    case FUNC_TYPE_REG_VANG:
      return FUNC_TYPE_REG_VANG_STR;
    case FUNC_TYPE_REG_PQ:
      return FUNC_TYPE_REG_PQ_STR;
    case FUNC_TYPE_REG_RATIO:
      return FUNC_TYPE_REG_RATIO_STR;
    case FUNC_TYPE_REG_PHASE:
      return FUNC_TYPE_REG_PHASE_STR;
    case FUNC_TYPE_REG_SUSC:
      return FUNC_TYPE_REG_SUSC_STR;
    case FUNC_TYPE_GEN_COST:
      return FUNC_TYPE_GEN_COST_STR;
    case FUNC_TYPE_SP_CONTROLS:
      return FUNC_TYPE_SP_CONTROLS_STR;
    case FUNC_TYPE_SLIM_VMAG:
      return FUNC_TYPE_SLIM_VMAG_STR;
    default:
      return FUNC_TYPE_UNKNOWN_STR;
    }
  else
    return FUNC_TYPE_UNKNOWN_STR;
}

REAL FUNC_get_weight(Func* f) {
  if (f)
    return f->weight;
  else
    return 0;
}

REAL FUNC_get_phi(Func* f) {
  if (f)
    return f->phi;
  else
    return 0;
}

REAL* FUNC_get_phi_ptr(Func* f) {
  if (f)
    return &(f->phi);
  else
    return NULL;
}

Vec* FUNC_get_gphi(Func* f) {
  if (f)
    return f->gphi;
  else
    return NULL;
}

Mat* FUNC_get_Hphi(Func* f) {
  if (f)
    return f->Hphi;
  else
    return NULL;
}

int FUNC_get_Hcounter(Func* f) {
  if (f)
    return f->Hcounter;
  else
    return 0;
}

int* FUNC_get_Hcounter_ptr(Func* f) {
  if (f)
    return &(f->Hcounter);
  else
    return NULL;
}

char* FUNC_get_bus_counted(Func* f) {
  if (f)
    return f->bus_counted;
  else
    return NULL;
}

int FUNC_get_bus_counted_size(Func* f) {
  if (f)
    return f->bus_counted_size;
  else
    return 0;
}

Func* FUNC_get_next(Func* f) {
  if (f)
    return f->next;
  else
    return NULL;
}

int FUNC_get_branch_counter(Func* f) {
  if (f)
    return f->branch_counter;
  else
    return 0;
}

void FUNC_inc_branch_counter(Func* f) {
  if (f)
    f->branch_counter++;
}

Func* FUNC_list_add(Func* flist, Func* nf) {
  LIST_add(flist,nf,next);
  return flist;
}

int FUNC_list_len(Func* flist) {
  int len;
  LIST_len(Func,flist,next,len);
  return len;
}

void FUNC_list_del(Func* flist) {
  LIST_map(Func,flist,f,next,{FUNC_del(f);});
}

void FUNC_list_count_branch(Func* f, Branch* b) {
  Func* ff;
  for (ff = f; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_count_branch(ff,b);
}

void FUNC_list_allocate(Func* f) {
  Func* ff;
  for (ff = f; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_allocate(ff);
}

void FUNC_list_clear(Func* f) {
  Func* ff;
  for (ff = f; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_clear(ff);
}

void FUNC_list_analyze_branch(Func* f, Branch* b) {
  Func* ff;
  for (ff = f; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_analyze_branch(ff,b);
}

void FUNC_list_eval_branch(Func* f, Branch* b, Vec* values) {
  Func* ff;
  for (ff = f; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_eval_branch(ff,b,values);
}

Func* FUNC_new(int type, REAL weight, Net* net) {

  Func* f = (Func*)malloc(sizeof(Func));

  // Error
  f->error_flag = FALSE;
  strcpy(f->error_string,"");

  // Network
  f->net = net;

  // Fields
  f->type = type;
  f->weight = weight;
  f->phi = 0;
  f->gphi = NULL;
  f->Hphi = NULL;
  f->Hcounter = 0;
  f->branch_counter = 0;
  f->next = NULL;

  // Bus counted
  f->bus_counted_size = 0;
  f->bus_counted = NULL;
  
  // Functions
  switch (type) {

  case FUNC_TYPE_REG_VMAG: // Volage magnitude regularization
    f->func_init = FUNC_REG_VMAG_init;
    f->func_count_branch = FUNC_REG_VMAG_count_branch;
    f->func_allocate = FUNC_REG_VMAG_allocate;
    f->func_clear = FUNC_REG_VMAG_clear;
    f->func_analyze_branch = FUNC_REG_VMAG_analyze_branch;
    f->func_eval_branch = FUNC_REG_VMAG_eval_branch;
    f->func_free = FUNC_REG_VMAG_free;
    break;

  case FUNC_TYPE_REG_VANG: // Volage angle regularization
    f->func_init = FUNC_REG_VANG_init;
    f->func_count_branch = FUNC_REG_VANG_count_branch;
    f->func_allocate = FUNC_REG_VANG_allocate;
    f->func_clear = FUNC_REG_VANG_clear;
    f->func_analyze_branch = FUNC_REG_VANG_analyze_branch;
    f->func_eval_branch = FUNC_REG_VANG_eval_branch;
    f->func_free = FUNC_REG_VANG_free;
    break;
    
  case FUNC_TYPE_REG_PQ: // Generator reactive power regularization
    f->func_init = FUNC_REG_PQ_init;
    f->func_count_branch = FUNC_REG_PQ_count_branch;
    f->func_allocate = FUNC_REG_PQ_allocate;
    f->func_clear = FUNC_REG_PQ_clear;
    f->func_analyze_branch = FUNC_REG_PQ_analyze_branch;
    f->func_eval_branch = FUNC_REG_PQ_eval_branch;
    f->func_free = FUNC_REG_PQ_free;
    break;

  case FUNC_TYPE_REG_RATIO: // Transformer tap ratio regularization
    f->func_init = FUNC_REG_RATIO_init;
    f->func_count_branch = FUNC_REG_RATIO_count_branch;
    f->func_allocate = FUNC_REG_RATIO_allocate;
    f->func_clear = FUNC_REG_RATIO_clear;
    f->func_analyze_branch = FUNC_REG_RATIO_analyze_branch;
    f->func_eval_branch = FUNC_REG_RATIO_eval_branch;
    f->func_free = FUNC_REG_RATIO_free;
    break;

  case FUNC_TYPE_REG_PHASE: // Transformer phase shift regularization
    f->func_init = FUNC_REG_PHASE_init;
    f->func_count_branch = FUNC_REG_PHASE_count_branch;
    f->func_allocate = FUNC_REG_PHASE_allocate;
    f->func_clear = FUNC_REG_PHASE_clear;
    f->func_analyze_branch = FUNC_REG_PHASE_analyze_branch;
    f->func_eval_branch = FUNC_REG_PHASE_eval_branch;
    f->func_free = FUNC_REG_PHASE_free;
    break;

  case FUNC_TYPE_REG_SUSC: // Transformer tap ratio regularization
    f->func_init = FUNC_REG_SUSC_init;
    f->func_count_branch = FUNC_REG_SUSC_count_branch;
    f->func_allocate = FUNC_REG_SUSC_allocate;
    f->func_clear = FUNC_REG_SUSC_clear;
    f->func_analyze_branch = FUNC_REG_SUSC_analyze_branch;
    f->func_eval_branch = FUNC_REG_SUSC_eval_branch;
    f->func_free = FUNC_REG_SUSC_free;
    break;

  case FUNC_TYPE_GEN_COST: // Power generation cost
    f->func_init = FUNC_GEN_COST_init;
    f->func_count_branch = FUNC_GEN_COST_count_branch;
    f->func_allocate = FUNC_GEN_COST_allocate;
    f->func_clear = FUNC_GEN_COST_clear;
    f->func_analyze_branch = FUNC_GEN_COST_analyze_branch;
    f->func_eval_branch = FUNC_GEN_COST_eval_branch;
    f->func_free = FUNC_GEN_COST_free;
    break;

  case FUNC_TYPE_SP_CONTROLS: // Sparse controls
    f->func_init = FUNC_SP_CONTROLS_init;
    f->func_count_branch = FUNC_SP_CONTROLS_count_branch;
    f->func_allocate = FUNC_SP_CONTROLS_allocate;
    f->func_clear = FUNC_SP_CONTROLS_clear;
    f->func_analyze_branch = FUNC_SP_CONTROLS_analyze_branch;
    f->func_eval_branch = FUNC_SP_CONTROLS_eval_branch;
    f->func_free = FUNC_SP_CONTROLS_free;
    break;

  case FUNC_TYPE_SLIM_VMAG: // Voltage magnitude soft limits
    f->func_init = FUNC_SLIM_VMAG_init;
    f->func_count_branch = FUNC_SLIM_VMAG_count_branch;
    f->func_allocate = FUNC_SLIM_VMAG_allocate;
    f->func_clear = FUNC_SLIM_VMAG_clear;
    f->func_analyze_branch = FUNC_SLIM_VMAG_analyze_branch;
    f->func_eval_branch = FUNC_SLIM_VMAG_eval_branch;
    f->func_free = FUNC_SLIM_VMAG_free;
    break;
 
  default:
    f->func_init = NULL;
    f->func_count_branch = NULL;
    f->func_allocate = NULL;
    f->func_clear = NULL;
    f->func_analyze_branch = NULL;
    f->func_eval_branch = NULL;
    f->func_free = NULL;
    break;
  }

  // Update network
  FUNC_update_network(f);
  
  return f;
}

void FUNC_set_phi(Func* f, REAL phi) {
  if (f)
    f->phi = phi;
}

void FUNC_set_gphi(Func* f, Vec* gphi) {
  if (f)
    f->gphi = gphi;
}

void FUNC_set_Hphi(Func* f, Mat* Hphi) {
  if (f)
    f->Hphi = Hphi;
}

void FUNC_set_Hcounter(Func* f, int counter) {
  if (f)
    f->Hcounter = counter;
}

void FUNC_set_bus_counted(Func* f, char* counted, int size) {
  if (f) {
    f->bus_counted = counted;
    f->bus_counted_size = size;
  }  
}

void FUNC_set_branch_counter(Func* f, int counter) {
  if (f)
    f->branch_counter = counter;
}

void FUNC_count(Func* f) {
  int i;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    FUNC_count_branch(f,NET_get_branch(net,i));
}

void FUNC_count_branch(Func* f, Branch* b) {
  if (f && f->func_count_branch && FUNC_is_safe_to_count(f))
    (*(f->func_count_branch))(f,b);
}

void FUNC_allocate(Func* f) {
  if (f && f->func_allocate && FUNC_is_safe_to_count(f))
    (*(f->func_allocate))(f);
}

void FUNC_clear(Func* f) {
  if (f && f->func_clear)
    (*(f->func_clear))(f);
}

void FUNC_analyze(Func* f) {
  int i;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    FUNC_analyze_branch(f,NET_get_branch(net,i));
}

void FUNC_analyze_branch(Func* f, Branch* b) {
  if (f && f->func_analyze_branch && FUNC_is_safe_to_analyze(f))
    (*(f->func_analyze_branch))(f,b);
}

void FUNC_eval(Func* f, Vec* values) {
  int i;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    FUNC_eval_branch(f,NET_get_branch(net,i),values);
}

void FUNC_eval_branch(Func* f, Branch* b, Vec* values) {
  if (f && f->func_eval_branch && FUNC_is_safe_to_eval(f,values))
    (*(f->func_eval_branch))(f,b,values);
}

BOOL FUNC_is_safe_to_count(Func* f) {
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(FUNC_get_network(f)))
    return TRUE;
  else {
    sprintf(f->error_string,"function is not safe to count");
    f->error_flag = TRUE;
    return FALSE;
  }  
}

BOOL FUNC_is_safe_to_analyze(Func* f) {
  Net* net = FUNC_get_network(f);
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(net) &&
      MAT_get_size1(f->Hphi) == NET_get_num_vars(net) &&
      MAT_get_size2(f->Hphi) == NET_get_num_vars(net))
    return TRUE;
  else {
    sprintf(f->error_string,"function is not safe to analyze");
    f->error_flag = TRUE;
    return FALSE;
  }
}

BOOL FUNC_is_safe_to_eval(Func* f, Vec* values) {
  Net* net = FUNC_get_network(f);
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(net) &&
      MAT_get_size1(f->Hphi) == NET_get_num_vars(net) &&
      MAT_get_size2(f->Hphi) == NET_get_num_vars(net) &&
      VEC_get_size(f->gphi) == NET_get_num_vars(net) &&
      VEC_get_size(values) == NET_get_num_vars(net))
    return TRUE;
  else {
    sprintf(f->error_string,"function is not safe to eval");
    f->error_flag = TRUE;
    return FALSE;
  }
}

BOOL FUNC_has_error(Func* f) {
  if (!f)
    return FALSE;
  else
    return f->error_flag;
}

void FUNC_clear_error(Func * f) {
  if (f) {
    f->error_flag = FALSE;
    strcpy(f->error_string,"");
  }
}

char* FUNC_get_error_string(Func* f) {
  if (!f)
    return NULL;
  else
    return f->error_string;
}

void FUNC_update_network(Func* f) {
  if (!f)
    return;
 
  // Bus counted
  if (f->bus_counted)
    free(f->bus_counted);
  f->bus_counted_size = NET_get_num_buses(f->net);
  f->bus_counted = (char*)calloc(NET_get_num_buses(f->net),sizeof(char));

  // Type-specific data
  if (f->func_free)
    (*(f->func_free))(f);
  if (f->func_init)
    (*(f->func_init))(f);
}

Net* FUNC_get_network(Func* f) {
  if (f)
    return f->net;
  else
    return NULL;
}
