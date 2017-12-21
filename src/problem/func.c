/** @file func.c
 *  @brief This file defines the Func data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/func.h>

struct Func {
  
  // Error
  BOOL error_flag;                     /**< @brief Error flag */
  char error_string[FUNC_BUFFER_SIZE]; /**< @brief Error string */

  // Name
  char name[FUNC_BUFFER_SIZE]; /**< @brief Name string */

  // Network
  Net* net;    /**< @brief Power network */
  
  // Weight
  REAL weight; /**< @brief Function weight for forming an objective function */

  // Value
  REAL phi;    /**< @brief Function value */

  // Gradient
  Vec* gphi;   /**< @brief Function gradient vector */

  // Hessian
  Mat* Hphi;   /**< @brief Function Hessian matrix (only lower triangular part stored) */

  // Counters and flags
  int Hphi_nnz;         /**< @brief Counter of number of nonzero elements of the Hessian matrix */
  char* bus_counted;    /**< @brief Flags for processing buses */
  int bus_counted_size; /**< @brief Size of array of flags for processing buses */
  
  // Functions
  void (*func_init)(Func* f);                                    /**< @brief Initialization function */
  void (*func_count_step)(Func* f, Branch* br, int t);           /**< @brief Function for countinng nonzero entries */
  void (*func_allocate)(Func* f);                                /**< @brief Function for allocating required arrays */
  void (*func_clear)(Func* f);                                   /**< @brief Function for clearing flags, counters, and function values */
  void (*func_analyze_step)(Func* f, Branch* br, int t);         /**< @brief Function for analyzing sparsity pattern */
  void (*func_eval_step)(Func* f, Branch* br, int t, Vec* v);    /**< @brief Function for evaluating function */
  void (*func_free)(Func* f);                                    /**< @brief Function for de-allocating any data used */
  void (*func_set_parameter)(Func* c, char* key, void* value);   /**< @brief Function for setting function parameter */

  // Custom data
  void* data;  /**< @brief Type-dependent function data */

  // List
  Func* next; /**< @brief List of functions for forming objective function */
};

void FUNC_clear_bus_counted(Func* f) {
  if (f)
    ARRAY_clear(f->bus_counted,char,f->bus_counted_size);
}

void FUNC_del_matvec(Func* f) {
  if (f) {

    // Mat and vec
    VEC_del(f->gphi);
    MAT_del(f->Hphi);
    f->gphi = NULL;
    f->Hphi = NULL;
  }
}

void FUNC_del(Func* f) {
  if (f) {
    
    // Mat and vec
    FUNC_del_matvec(f);

    // Utils
    if (f->bus_counted)
      free(f->bus_counted);

    // Data
    if (f->func_free)
      (*(f->func_free))(f);

    free(f);
  }
}

void FUNC_finalize_structure_of_Hessian(Func* f) {

  // Local variables
  int* Hi;
  int* Hj;
  int temp;
  int m;
  
  // Check
  if (!f)
    return;

  // Ensure lower triangular
  Hi = MAT_get_row_array(f->Hphi);
  Hj = MAT_get_col_array(f->Hphi);
  for (m = 0; m < MAT_get_nnz(f->Hphi); m++) {
    if (Hi[m] < Hj[m]) {
      temp = Hi[m];
      Hi[m] = Hj[m];
      Hj[m] = temp;
    }
  }
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

int FUNC_get_Hphi_nnz(Func* f) {
  if (f)
    return f->Hphi_nnz;
  else
    return 0;
}

int* FUNC_get_Hphi_nnz_ptr(Func* f) {
  if (f)
    return &(f->Hphi_nnz);
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

void* FUNC_get_data(Func* f) {
  if (f)
    return f->data;
  else
    return NULL;
}

void FUNC_list_clear_error(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_clear_error(ff);
}

BOOL FUNC_list_has_error(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff)) {
    if (FUNC_has_error(ff))
      return TRUE;
  }
  return FALSE;
}

char* FUNC_list_get_error_string(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff)) {
    if (FUNC_has_error(ff))
      return FUNC_get_error_string(ff);
  }
  return "";
}

Func* FUNC_list_add(Func* flist, Func* nf) {
  LIST_add(Func,flist,nf,next);
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

void FUNC_list_count_step(Func* flist, Branch* br, int t) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_count_step(ff,br,t);
}

void FUNC_list_allocate(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_allocate(ff);
}

void FUNC_list_clear(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_clear(ff);
}

void FUNC_list_analyze_step(Func* flist, Branch* br, int t) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_analyze_step(ff,br,t);
}

void FUNC_list_eval_step(Func* flist, Branch* br, int t, Vec* values) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_eval_step(ff,br,t,values);
}

void FUNC_list_finalize_structure_of_Hessian(Func* flist) {
  Func* ff;
  for (ff = flist; ff != NULL; ff = FUNC_get_next(ff))
    FUNC_finalize_structure_of_Hessian(ff);
}

Func* FUNC_new(REAL weight, Net* net) {

  Func* f = (Func*)malloc(sizeof(Func));

  // Error
  f->error_flag = FALSE;
  strcpy(f->error_string,"");

  // Network
  f->net = net;

  // Name
  strcpy(f->name,"unknown");

  // Fields
  f->weight = weight;
  f->phi = 0;
  f->gphi = NULL;
  f->Hphi = NULL;
  f->Hphi_nnz = 0;
  f->next = NULL;

  // Bus counted
  f->bus_counted_size = 0;
  f->bus_counted = NULL;
  
  // Methods
  f->func_init = NULL;
  f->func_count_step = NULL;
  f->func_allocate = NULL;
  f->func_clear = NULL;
  f->func_analyze_step = NULL;
  f->func_eval_step = NULL;
  f->func_free = NULL;
  f->func_set_parameter = NULL;

  // Data
  f->data = NULL;

  // Update network
  FUNC_update_network(f);
  
  return f;
}

void FUNC_set_parameter(Func* f, char* key, void* value) {
  if (f) {
    if (f->func_set_parameter)
      (*(f->func_set_parameter))(f, key, value);
    else {
      sprintf(f->error_string,"function does not support setting parameters");
      f->error_flag = TRUE;
      return;
    }
  }
}

void FUNC_set_name(Func* f, char* name) {
  if (f)
    strcpy(f->name,name);
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

void FUNC_set_Hphi_nnz(Func* f, int nnz) {
  if (f)
    f->Hphi_nnz = nnz;
}

void FUNC_set_bus_counted(Func* f, char* bus_counted, int size) {
  if (f) {
    if (f->bus_counted)
      free(f->bus_counted);
    f->bus_counted = bus_counted;
    f->bus_counted_size = size;
  }  
}

void FUNC_set_data(Func* f, void* data) {
  if (f)
    f->data = data;
}

void FUNC_count(Func* f) {
  int i;
  int t;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      FUNC_count_step(f,NET_get_branch(net,i),t);
  }
}

void FUNC_count_step(Func* f, Branch* br, int t) {
  if (f && f->func_count_step && FUNC_is_safe_to_count(f))
    (*(f->func_count_step))(f,br,t);
}

void FUNC_allocate(Func* f) {
  if (f && f->func_allocate && FUNC_is_safe_to_count(f)) {
    FUNC_del_matvec(f);
    (*(f->func_allocate))(f);
  }
}

void FUNC_init(Func* f) {
  if (f && f->func_free)
    (*(f->func_free))(f);
  if (f && f->func_init)
    (*(f->func_init))(f);
}

void FUNC_clear(Func* f) {
  if (f && f->func_clear)
    (*(f->func_clear))(f);
}

void FUNC_analyze(Func* f) {
  int i;
  int t;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      FUNC_analyze_step(f,NET_get_branch(net,i),t);
  }
  FUNC_finalize_structure_of_Hessian(f);
}

void FUNC_analyze_step(Func* f, Branch* br, int t) {
  if (f && f->func_analyze_step && FUNC_is_safe_to_analyze(f))
    (*(f->func_analyze_step))(f,br,t);
}

void FUNC_eval(Func* f, Vec* values) {
  int i;
  int t;
  Net* net = FUNC_get_network(f);
  FUNC_clear(f);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      FUNC_eval_step(f,NET_get_branch(net,i),t,values);
  }
}

void FUNC_eval_step(Func* f, Branch* br, int t, Vec* values) {
  if (f && f->func_eval_step && FUNC_is_safe_to_eval(f,values))
    (*(f->func_eval_step))(f,br,t,values);
}

BOOL FUNC_is_safe_to_count(Func* f) {
  Net* net = FUNC_get_network(f);
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(net)*NET_get_num_periods(net))
    return TRUE;
  else {
    sprintf(f->error_string,"function is not safe to count");
    f->error_flag = TRUE;
    return FALSE;
  }  
}

BOOL FUNC_is_safe_to_analyze(Func* f) {
  Net* net = FUNC_get_network(f);
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
      VEC_get_size(f->gphi) == NET_get_num_vars(net) && 
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
  if (FUNC_get_bus_counted_size(f) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
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
  if (f)
    return f->error_flag;
  else
    return FALSE;
}

void FUNC_clear_error(Func * f) {
  if (f) {
    f->error_flag = FALSE;
    strcpy(f->error_string,"");
  }
}

char* FUNC_get_error_string(Func* f) {
  if (f)
    return f->error_string;
  else
    return NULL;
}

char* FUNC_get_name(Func* f) {
  if (f)
    return f->name;
  else
    return NULL;
}

void FUNC_update_network(Func* f) {
  
  // No f
  if (!f)
    return;
 
  // Bus counted
  if (f->bus_counted)
    free(f->bus_counted);
  f->bus_counted_size = NET_get_num_buses(f->net)*NET_get_num_periods(f->net);
  ARRAY_zalloc(f->bus_counted,char,f->bus_counted_size);

  // Init
  FUNC_init(f);
}

Net* FUNC_get_network(Func* f) {
  if (f)
    return f->net;
  else
    return NULL;
}

void FUNC_set_func_init(Func* f, void (*func)(Func* f)) {
  if (f)
    f->func_init = func;
}

void FUNC_set_func_count_step(Func* f, void (*func)(Func* f, Branch* br, int t)) {
  if (f)
    f->func_count_step = func;
}

void FUNC_set_func_allocate(Func* f, void (*func)(Func* f)) {
  if (f)
    f->func_allocate = func;
}

void FUNC_set_func_clear(Func* f, void (*func)(Func* f)) {
  if (f)
    f->func_clear = func;
}

void FUNC_set_func_analyze_step(Func* f, void (*func)(Func* f, Branch* br, int t)) {
  if (f)
    f->func_analyze_step = func;
}

void FUNC_set_func_eval_step(Func* f, void (*func)(Func* f, Branch* br, int t, Vec* v)) {
  if (f)
    f->func_eval_step = func;
}

void FUNC_set_func_free(Func* f, void (*func)(Func* f)) {
  if (f)
    f->func_free = func;
}

void FUNC_set_func_set_parameter(Func* f, void (*func)(Func* f, char* key, void* value)) {
  if (f)
    f->func_set_parameter = func;
}
