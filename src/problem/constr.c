/** @file constr.c
 *  @brief This file defines the Constr data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr.h>
#include <pfnet/constr_PF.h>
#include <pfnet/constr_FIX.h>
#include <pfnet/constr_BOUND.h>
#include <pfnet/constr_PAR_GEN.h>
#include <pfnet/constr_REG_GEN.h>
#include <pfnet/constr_REG_TRAN.h>
#include <pfnet/constr_REG_SHUNT.h>

struct Constr {

  // Error
  BOOL error_flag;                       /**< @brief Error flag */
  char error_string[CONSTR_BUFFER_SIZE]; /**< @brief Error string */

  // Network
  Net* net;    /**< @brief Power network */

  // Type
  int type; /**< @brief Constraint type */

  // Nonlinear
  Vec* f;           /**< @brief Vector of nonlinear constraint violations */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints */
  Mat* H_array;     /**< @brief Array of Hessian matrices of nonlinear constraints */
  int H_array_size; /**< @brief Size of Hessian array */
  Mat* H_combined;  /**< @brief Linear combination of Hessians of the nonlinear constraints */

  // Linear
  Mat* A; /**< @brief Matrix of constraint normals of linear cosntraints */
  Vec* b; /**< @brief Right-hand side vector of linear constraints */

  // Utils
  int Acounter;         /**< @brief Counter for nonzeros of matrix A */
  int Jcounter;         /**< @brief Counter for nonzeros matrix J */
  int* Hcounter;        /**< @brief Array of counters of nonzeros of nonlinear constraint Hessians */
  int Hcounter_size;    /**< @brief Size of array of counter of Hessian nonzeros */
  int Aconstr_index;    /**< @brief Index for linear constraints */
  int Jconstr_index;    /**< @brief Index for nonlinear constraints */
  char* bus_counted;    /**< @brief Flag for processing buses */
  int bus_counted_size; /**< @brief Size of array of flags for processing buses */
  int branch_counter;   /**< @brief Counter for processing branches */

  // Type functions
  void (*func_init)(Constr* c); /**< @brief Initialization function */
  void (*func_count_branch)(Constr* c, Branch* br); /**< @brief Function for counting nonzero entries */
  void (*func_allocate)(Constr* c); /**< @brief Function for allocating required arrays */
  void (*func_clear)(Constr* c); /**< @brief Function for clearing flags, counters, and function values */
  void (*func_analyze_branch)(Constr* c, Branch* br); /**< @brief Function for analyzing sparsity pattern */
  void (*func_eval_branch)(Constr* c, Branch* br, Vec* var_values); /**< @brief Function for evaluating constraint */
  void (*func_store_sens_branch)(Constr* c, Branch* br, Vec* sens); /**< @brief Function for storing sensitivities */
  void (*func_free)(Constr* c); /**< @brief Function for de-allocating any data used */

  // Type data
  void* data; /**< @brief Type-dependent constraint data structure */

  // List
  Constr* next; /**< @brief List of constraints */
};

void CONSTR_clear_Hcounter(Constr* c) {
  int i;
  if (c) {
    if (c->Hcounter) {
      for (i = 0; i < c->Hcounter_size; i++)
	c->Hcounter[i] = 0;
    } 
  }
}

void CONSTR_clear_bus_counted(Constr* c) {
  int i;
  if (c) {
    if (c->bus_counted) {
      for (i = 0; i < c->bus_counted_size; i++)
	c->bus_counted[i] = 0;
    }
  } 
}

void CONSTR_combine_H(Constr* c, Vec* coeff, BOOL ensure_psd) {
  
  // Local variabels
  REAL* Hd;
  REAL* Hd_comb;
  REAL* coeffd;
  REAL coeffk;
  int Hcounter_comb;
  int k;
  int m;

  if (!c)
    return;

  // Check dimensions
  if (VEC_get_size(coeff) != c->H_array_size) {
    sprintf(c->error_string,"invalid dimensions");
    c->error_flag = TRUE;
  }
  
  // Combine
  Hcounter_comb = 0;
  coeffd = VEC_get_data(coeff);
  Hd_comb = MAT_get_data_array(c->H_combined);
  for (k = 0; k < c->H_array_size; k++) {
    Hd = MAT_get_data_array(MAT_array_get(c->H_array,k));
    if (ensure_psd)
      coeffk = 0;
    else
      coeffk = coeffd[k];
    for (m = 0; m < MAT_get_nnz(MAT_array_get(c->H_array,k)); m++) {
      Hd_comb[Hcounter_comb] = coeffk*Hd[m];
      Hcounter_comb++;
    }
  }
}

void CONSTR_del(Constr* c) {
  if (c) {

    // Mat and vec
    VEC_del(c->b);
    MAT_del(c->A);
    VEC_del(c->f);
    MAT_del(c->J);
    MAT_array_del(c->H_array,c->H_array_size);
    MAT_del(c->H_combined);

    // Utils
    if (c->bus_counted)
      free(c->bus_counted);

    // Data
    if (c->func_free)
      (*(c->func_free))(c);

    free(c);
  }
}

int CONSTR_get_type(Constr* c) {
  if (c)
    return c->type;
  else
    return CONSTR_TYPE_UNKNOWN;
}

Vec* CONSTR_get_b(Constr* c) {
  if (c)
    return c->b;
  else
    return NULL;
}

Mat* CONSTR_get_A(Constr* c) {
  if (c)
    return c->A;
  else
    return NULL;
}

Vec* CONSTR_get_f(Constr* c) {
  if (c)
    return c->f;
  else
    return NULL;
}

Mat* CONSTR_get_J(Constr* c) {
  if (c)
    return c->J;
  else
    return NULL;
}

Mat* CONSTR_get_H_array(Constr* c) {
  if (c)
    return c->H_array;
  else
    return NULL;
}

int CONSTR_get_H_array_size(Constr* c) {
  if (c)
    return c->H_array_size;
  else
    return 0;
}

Mat* CONSTR_get_H_single(Constr* c, int i) {
  if (c && 0 <= i && i < c->H_array_size)
    return MAT_array_get(c->H_array,i);
  else
    return NULL;
}

Mat* CONSTR_get_H_combined(Constr* c) {
  if (c)
    return c->H_combined;
  else
    return NULL;
}

int CONSTR_get_Acounter(Constr* c) {
  if (c)
    return c->Acounter;
  else
    return 0;
}

int* CONSTR_get_Acounter_ptr(Constr* c) {
  if (c)
    return &(c->Acounter);
  else
    return NULL;
}

int CONSTR_get_Jcounter(Constr* c) {
  if (c)
    return c->Jcounter;
  else
    return 0;
}

int* CONSTR_get_Jcounter_ptr(Constr* c) {
  if (c)
    return &(c->Jcounter);
  else
    return 0;
}

int* CONSTR_get_Hcounter(Constr* c) {
  if (c)
    return c->Hcounter;
  else
    return NULL;
}

int CONSTR_get_Hcounter_size(Constr* c) {
  if (c)
    return c->Hcounter_size;
  else
    return 0;
}

int CONSTR_get_Aconstr_index(Constr* c) {
  if (c)
    return c->Aconstr_index;
  else
    return 0;
}

int CONSTR_get_Jconstr_index(Constr* c) {
  if (c)
    return c->Jconstr_index;
  else
    return 0;
}

int* CONSTR_get_Aconstr_index_ptr(Constr* c) {
  if (c)
    return &(c->Aconstr_index);
  else
    return NULL;
}

int* CONSTR_get_Jconstr_index_ptr(Constr* c) {
  if (c)
    return &(c->Jconstr_index);
  else
    return NULL;
}

char* CONSTR_get_bus_counted(Constr* c) {
  if (c)
    return c->bus_counted;
  else
    return NULL;
}

int CONSTR_get_bus_counted_size(Constr* c) {
  if (c)
    return c->bus_counted_size;
  else
    return 0;
}

void* CONSTR_get_data(Constr* c) {
  if (c)
    return c->data;
  else
    return NULL;
}

Constr* CONSTR_get_next(Constr* c) {
  if (c)
    return c->next;
  else
    return NULL;
}

int CONSTR_get_branch_counter(Constr* c) {
  if (c)
    return c->branch_counter;
  else
    return 0;
}

void CONSTR_inc_branch_counter(Constr* c) {
  if (c)
    c->branch_counter++;
}

Constr* CONSTR_list_add(Constr* clist, Constr* nc) {
  LIST_add(clist,nc,next);
  return clist;
}

int CONSTR_list_len(Constr* clist) {
  int len;
  LIST_len(Constr,clist,next,len);
  return len;
}

void CONSTR_list_del(Constr* clist) {
  LIST_map(Constr,clist,c,next,{CONSTR_del(c);});
}

void CONSTR_list_combine_H(Constr* clist, Vec* coeff, BOOL ensure_psd) {
  Constr* cc;
  Vec* v;
  int size = 0;
  int offset = 0;
  REAL* coeffd = VEC_get_data(coeff);

  // Size
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    size += VEC_get_size(CONSTR_get_f(cc));
    
  // Map
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (offset + VEC_get_size(CONSTR_get_f(cc)) <= VEC_get_size(coeff))
      v = VEC_new_from_array(&(coeffd[offset]),VEC_get_size(CONSTR_get_f(cc)));
    else
      v = NULL;       
    CONSTR_combine_H(cc,v,ensure_psd);
    offset += VEC_get_size(CONSTR_get_f(cc));
  }
}

void CONSTR_list_count_branch(Constr* clist, Branch* b) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_count_branch(cc,b);
}

void CONSTR_list_allocate(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_allocate(cc);
}

void CONSTR_list_clear(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_clear(cc);
}

void CONSTR_list_analyze_branch(Constr* clist, Branch* b) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_analyze_branch(cc,b);
}

void CONSTR_list_eval_branch(Constr* clist, Branch* br, Vec* values) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_eval_branch(cc,br,values);
}

void CONSTR_list_store_sens_branch(Constr* clist, Branch* br, Vec* sens) {
  Constr* cc;
  Vec* v;
  int size = 0;
  int offset = 0;
  REAL* sensd = VEC_get_data(sens);
  
  // Size
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    size += VEC_get_size(CONSTR_get_f(cc));
  
  // Map
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (offset + VEC_get_size(CONSTR_get_f(cc)) <= VEC_get_size(sens))
      v = VEC_new_from_array(&(sensd[offset]),VEC_get_size(CONSTR_get_f(cc)));
    else
      v = NULL;
    CONSTR_store_sens_branch(cc,br,v);
    offset += VEC_get_size(CONSTR_get_f(cc));
  }
}

Constr* CONSTR_new(int type, Net* net) {

  Constr* c = (Constr*)malloc(sizeof(Constr));

  // Error
  c->error_flag = FALSE;
  strcpy(c->error_string,"");

  // Network
  c->net = net;

  // Fields
  c->type = type;
  c->f = NULL;
  c->J = NULL;
  c->H_array = NULL;  
  c->H_array_size = 0;
  c->H_combined = NULL;
  c->A = NULL;
  c->b = NULL;
  c->Acounter = 0;
  c->Jcounter = 0;
  c->Hcounter = NULL;
  c->Hcounter_size = 0;
  c->Aconstr_index = 0;
  c->Jconstr_index = 0;
  c->branch_counter = 0;
  c->data = NULL;
  c->next = NULL;

  // Bus counted flags
  c->bus_counted_size = 0;
  c->bus_counted = NULL;

  // Constraint functions
  if (type == CONSTR_TYPE_PF) { // power flow
    c->func_init = &CONSTR_PF_init;
    c->func_count_branch = &CONSTR_PF_count_branch;
    c->func_allocate = &CONSTR_PF_allocate;
    c->func_clear = &CONSTR_PF_clear;
    c->func_analyze_branch = &CONSTR_PF_analyze_branch;
    c->func_eval_branch = &CONSTR_PF_eval_branch;
    c->func_store_sens_branch = &CONSTR_PF_store_sens_branch;
    c->func_free = &CONSTR_PF_free;
  }
  else if (type == CONSTR_TYPE_PAR_GEN) { // generator participation
    c->func_init = &CONSTR_PAR_GEN_init;
    c->func_count_branch = &CONSTR_PAR_GEN_count_branch;
    c->func_allocate = &CONSTR_PAR_GEN_allocate;
    c->func_clear = &CONSTR_PAR_GEN_clear;
    c->func_analyze_branch = &CONSTR_PAR_GEN_analyze_branch;
    c->func_eval_branch = &CONSTR_PAR_GEN_eval_branch;
    c->func_store_sens_branch = &CONSTR_PAR_GEN_store_sens_branch;
    c->func_free = &CONSTR_PAR_GEN_free;
  }
  else if (type == CONSTR_TYPE_FIX) { // variable fixing
    c->func_init = &CONSTR_FIX_init;
    c->func_count_branch = &CONSTR_FIX_count_branch;
    c->func_allocate = &CONSTR_FIX_allocate;
    c->func_clear = &CONSTR_FIX_clear;
    c->func_analyze_branch = &CONSTR_FIX_analyze_branch;
    c->func_eval_branch = &CONSTR_FIX_eval_branch;
    c->func_store_sens_branch = &CONSTR_FIX_store_sens_branch;
    c->func_free = &CONSTR_FIX_free;
  }
  else if (type == CONSTR_TYPE_REG_GEN) { // voltage regulation by generator
    c->func_init = &CONSTR_REG_GEN_init;
    c->func_count_branch = &CONSTR_REG_GEN_count_branch;
    c->func_allocate = &CONSTR_REG_GEN_allocate;
    c->func_clear = &CONSTR_REG_GEN_clear;
    c->func_analyze_branch = &CONSTR_REG_GEN_analyze_branch;
    c->func_eval_branch = &CONSTR_REG_GEN_eval_branch;
    c->func_store_sens_branch = &CONSTR_REG_GEN_store_sens_branch;
    c->func_free = &CONSTR_REG_GEN_free;
  }
  else if (type == CONSTR_TYPE_REG_TRAN) { // voltage regulation by transformer
    c->func_init = &CONSTR_REG_TRAN_init;
    c->func_count_branch = &CONSTR_REG_TRAN_count_branch;
    c->func_allocate = &CONSTR_REG_TRAN_allocate;
    c->func_clear = &CONSTR_REG_TRAN_clear;
    c->func_analyze_branch = &CONSTR_REG_TRAN_analyze_branch;
    c->func_eval_branch = &CONSTR_REG_TRAN_eval_branch;
    c->func_store_sens_branch = &CONSTR_REG_TRAN_store_sens_branch;
    c->func_free = &CONSTR_REG_TRAN_free;
  }
  else if (type == CONSTR_TYPE_REG_SHUNT) { // voltage regulation by shunt device
    c->func_init = &CONSTR_REG_SHUNT_init;
    c->func_count_branch = &CONSTR_REG_SHUNT_count_branch;
    c->func_allocate = &CONSTR_REG_SHUNT_allocate;
    c->func_clear = &CONSTR_REG_SHUNT_clear;
    c->func_analyze_branch = &CONSTR_REG_SHUNT_analyze_branch;
    c->func_eval_branch = &CONSTR_REG_SHUNT_eval_branch;
    c->func_store_sens_branch = &CONSTR_REG_SHUNT_store_sens_branch;
    c->func_free = &CONSTR_REG_SHUNT_free;
  }
  else if (type == CONSTR_TYPE_BOUND) { // variable bounds
    c->func_init = &CONSTR_BOUND_init;
    c->func_count_branch = &CONSTR_BOUND_count_branch;
    c->func_allocate = &CONSTR_BOUND_allocate;
    c->func_clear = &CONSTR_BOUND_clear;
    c->func_analyze_branch = &CONSTR_BOUND_analyze_branch;
    c->func_eval_branch = &CONSTR_BOUND_eval_branch;
    c->func_store_sens_branch = &CONSTR_BOUND_store_sens_branch;
    c->func_free = &CONSTR_BOUND_free;
  }
  else { // unknown 
    c->func_init = NULL;
    c->func_count_branch = NULL;
    c->func_allocate = NULL;
    c->func_clear = NULL;
    c->func_analyze_branch = NULL;
    c->func_eval_branch = NULL;
    c->func_store_sens_branch = NULL;
    c->func_free = NULL;
  }

  // Update network
  CONSTR_update_network(c);
  
  return c;
}

void CONSTR_set_b(Constr* c, Vec* b) {
  if (c)
    c->b = b;
}

void CONSTR_set_A(Constr* c, Mat* A) {
  if (c)
    c->A = A;
}

void CONSTR_set_f(Constr* c, Vec* f) {
  if (c)
    c->f = f;
}

void CONSTR_set_J(Constr* c, Mat* J) {
  if (c)
    c->J = J;
}

void CONSTR_set_H_array(Constr* c, Mat* array, int size) {
  if (c) {
    c->H_array = array;
    c->H_array_size = size;
  }  
}

void CONSTR_set_H_combined(Constr* c, Mat* H_combined) {
  if (c)
    c->H_combined = H_combined;
}

void CONSTR_set_Acounter(Constr* c, int counter) {
  if (c)
    c->Acounter = counter;
}

void CONSTR_set_Jcounter(Constr* c, int counter) {
  if (c)
    c->Jcounter = counter;
}

void CONSTR_set_Hcounter(Constr* c, int* counter, int size) {
  if (c) {
    c->Hcounter = counter;
    c->Hcounter_size = size;
  }
}

void CONSTR_set_Aconstr_index(Constr* c, int index) {
  if (c)
    c->Aconstr_index = index;
}

void CONSTR_set_Jconstr_index(Constr* c, int index) {
  if (c)
    c->Jconstr_index = index;
}

void CONSTR_set_bus_counted(Constr* c, char* counted, int size) {
  if (c) {
    c->bus_counted = counted;
    c->bus_counted_size = size;
  }
}

void CONSTR_set_branch_counter(Constr* c, int counter) {
  if (c)
    c->branch_counter = counter;
}

void CONSTR_set_data(Constr* c, void* data) {
  if (c)
    c->data = data;
}

void CONSTR_count(Constr* c) {
  int i;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    CONSTR_count_branch(c,NET_get_branch(net,i));
}

void CONSTR_count_branch(Constr* c, Branch* br) {
  if (c && c->func_count_branch && CONSTR_is_safe_to_count(c))
    (*(c->func_count_branch))(c,br);
}

void CONSTR_allocate(Constr* c) {
  if (c && c->func_allocate && CONSTR_is_safe_to_count(c))
    (*(c->func_allocate))(c);
}

void CONSTR_clear(Constr* c) {
  if (c && c->func_clear)
    (*(c->func_clear))(c);
}

void CONSTR_analyze(Constr* c) {
  int i;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    CONSTR_analyze_branch(c,NET_get_branch(net,i));
}

void CONSTR_analyze_branch(Constr* c, Branch* br) {
  if (c && c->func_analyze_branch && CONSTR_is_safe_to_analyze(c))
    (*(c->func_analyze_branch))(c,br);
}

void CONSTR_eval(Constr* c, Vec* values) {
  int i;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    CONSTR_eval_branch(c,NET_get_branch(net,i),values);
}

void CONSTR_eval_branch(Constr* c, Branch *br, Vec* values) {
  if (c && c->func_eval_branch && CONSTR_is_safe_to_eval(c,values))
    (*(c->func_eval_branch))(c,br,values);
}

void CONSTR_store_sens(Constr* c, Vec* sens) {
  int i;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (i = 0; i < NET_get_num_branches(net); i++) 
    CONSTR_store_sens_branch(c,NET_get_branch(net,i),sens);
}

void CONSTR_store_sens_branch(Constr* c, Branch *br, Vec* sens) {
  if (c && c->func_store_sens_branch && CONSTR_is_safe_to_count(c))
    (*(c->func_store_sens_branch))(c,br,sens);
}

BOOL CONSTR_is_safe_to_count(Constr* c) {
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(CONSTR_get_network(c)))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to count");
    c->error_flag = TRUE;
    return FALSE;
  }  
}

BOOL CONSTR_is_safe_to_analyze(Constr* c) {
  Net* net = CONSTR_get_network(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net) &&
      MAT_get_size2(c->A) == NET_get_num_vars(net) &&
      MAT_get_size2(c->J) == NET_get_num_vars(net)) {
    return TRUE;
  }
  else {
    sprintf(c->error_string,"constraint is not safe to analyze");
    c->error_flag = TRUE;
    return FALSE;
  }
}

BOOL CONSTR_is_safe_to_eval(Constr* c, Vec* values) {
  Net* net = CONSTR_get_network(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net) &&
      MAT_get_size2(c->A) == NET_get_num_vars(net) &&
      MAT_get_size2(c->J) == NET_get_num_vars(net) &&
      VEC_get_size(values) == NET_get_num_vars(net))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to eval");
    c->error_flag = TRUE;
    return FALSE;
  }
}

BOOL CONSTR_has_error(Constr* c) {
  if (!c)
    return FALSE;
  else
    return c->error_flag;
}

void CONSTR_set_error(Constr* c, char* string) {
  if (c) {
    c->error_flag = TRUE;
    strcpy(c->error_string,string);
  }
}

void CONSTR_clear_error(Constr * c) {
  if (c) {
    c->error_flag = FALSE;
    strcpy(c->error_string,"");
  }
}

char* CONSTR_get_error_string(Constr* c) {
  if (!c)
    return NULL;
  else
    return c->error_string;
}

void CONSTR_update_network(Constr* c) {
  if (!c)
    return;

  // Bus counted
  if (c->bus_counted)
    free(c->bus_counted);
  c->bus_counted_size = NET_get_num_buses(c->net);
  c->bus_counted = (char*)calloc(NET_get_num_buses(c->net),sizeof(char));

  // Type-specific data
  if (c->func_free)
    (*(c->func_free))(c);
  if (c->func_init)
    (*(c->func_init))(c);
}

Net* CONSTR_get_network(Constr* c) {
  if (c)
    return c->net;
  else
    return NULL;
}

