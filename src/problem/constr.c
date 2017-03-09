/** @file constr.c
 *  @brief This file defines the Constr data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr.h>
#include <pfnet/constr_PF.h>
#include <pfnet/constr_DCPF.h>
#include <pfnet/constr_LINPF.h>
#include <pfnet/constr_FIX.h>
#include <pfnet/constr_BOUND.h>
#include <pfnet/constr_PAR_GEN_P.h>
#include <pfnet/constr_PAR_GEN_Q.h>
#include <pfnet/constr_REG_GEN.h>
#include <pfnet/constr_REG_TRAN.h>
#include <pfnet/constr_REG_SHUNT.h>
#include <pfnet/constr_DC_FLOW_LIM.h>
#include <pfnet/constr_AC_FLOW_LIM.h>
#include <pfnet/constr_LBOUND.h>
#include <pfnet/constr_GEN_RAMP.h>

struct Constr {

  // Error
  BOOL error_flag;                       /**< @brief Error flag */
  char error_string[CONSTR_BUFFER_SIZE]; /**< @brief Error string */
  
  // Network
  Net* net;    /**< @brief Power network */

  // Type
  int type;                           /**< @brief Constraint type */
  char type_str[CONSTR_BUFFER_SIZE];  /**< @brief Constraint type string */
  
  // Nonlinear (f(x) + Jbar y = 0)
  Vec* f;           /**< @brief Vector of nonlinear constraint violations */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints wrt. variables */
  Mat* Jbar;        /**< @brief Jacobian matrix of nonlinear constraints wrt. extra variables */
  Mat* H_array;     /**< @brief Array of Hessian matrices of nonlinear constraints */
  int H_array_size; /**< @brief Size of Hessian array */
  Mat* H_combined;  /**< @brief Linear combination of Hessians of the nonlinear constraints */
  
  // Linear equality (Ax = b)
  Mat* A;           /**< @brief Matrix of constraint normals of linear equality constraints */
  Vec* b;           /**< @brief Right-hand side vector of linear equality constraints */

  // Linear inequalities (l <= Gx + Gbar y <= h)
  Mat* G;           /** @brief Matrix of constraint normals of linear inequality constraints wrt. variables */
  Mat* Gbar;        /** @brief Matrix of constraint normals of linear inequality constraints wrt. extra variables */
  Vec* l;           /** @brief Lower bound for linear inequality contraints */
  Vec* u;           /** @brief Upper bound for linear inequality contraints */
  
  // Extra variables
  int num_extra_vars;          /** @brief Number of extra variables (set during count) */
  
  // Counters and flags
  int A_nnz;             /**< @brief Counter for nonzeros of matrix A */
  int J_nnz;             /**< @brief Counter for nonzeros of matrix J */
  int Jbar_nnz;          /**< @brief Counter for nonzeros of matrix Jbar */
  int G_nnz;             /**< @brief Counter for nonzeros of matrix G */
  int Gbar_nnz;          /**< @brief Counter for nonzeros of matrix Gbar */
  int* H_nnz;            /**< @brief Array of counters of nonzeros of nonlinear constraint Hessians */
  int H_nnz_size;        /**< @brief Size of array of counter of Hessian nonzeros */
  int A_row;             /**< @brief Counter for linear equality constraints */
  int J_row;             /**< @brief Counter for nonlinear constraints */
  int G_row;             /**< @brief Counter for linear inequality constraints */
  char* bus_counted;     /**< @brief Flag for processing buses */
  int bus_counted_size;  /**< @brief Size of array of flags for processing buses */
  
  // Type functions
  void (*func_init)(Constr* c);                                       /**< @brief Initialization function */
  void (*func_count_step)(Constr* c, Branch* br, int t);              /**< @brief Function for counting nonzero entries */
  void (*func_allocate)(Constr* c);                                   /**< @brief Function for allocating required arrays */
  void (*func_clear)(Constr* c);                                      /**< @brief Function for clearing flags, counters, and function values */
  void (*func_analyze_step)(Constr* c, Branch* br, int t);            /**< @brief Function for analyzing sparsity pattern */
  void (*func_eval_step)(Constr* c, Branch* br, int t, Vec* v);       /**< @brief Function for evaluating constraint */
  void (*func_store_sens_step)(Constr* c, Branch* br, int t,
			       Vec* sA, Vec* sf, Vec* sGu, Vec* sGl); /**< @brief Func. for storing sensitivities */
  void (*func_free)(Constr* c);                                       /**< @brief Function for de-allocating any data used */

  // Type data
  void* data; /**< @brief Type-dependent constraint data structure */

  // List
  Constr* next; /**< @brief List of constraints */
};

int CONSTR_get_num_extra_vars(Constr* c) {
  if (c)
    return c->num_extra_vars;
  else
    return 0;
}

void CONSTR_set_num_extra_vars(Constr* c, int num) {
  if (c)
    c->num_extra_vars = num;
}

void CONSTR_clear_H_nnz(Constr* c) {
  if (c)
    ARRAY_clear(c->H_nnz,int,c->H_nnz_size);
}

void CONSTR_clear_bus_counted(Constr* c) {
  if (c)
    ARRAY_clear(c->bus_counted,char,c->bus_counted_size);
}

void CONSTR_combine_H(Constr* c, Vec* coeff, BOOL ensure_psd) {
  
  // Local variabels
  REAL* Hd;
  REAL* Hd_comb;
  REAL* coeffd;
  REAL coeffk;
  int H_nnz_comb;
  int k;
  int m;

  // No c
  if (!c)
    return;

  // Check dimensions
  if (VEC_get_size(coeff) != c->H_array_size) {
    sprintf(c->error_string,"invalid dimensions");
    c->error_flag = TRUE;
    return;
  }
  
  // Combine
  H_nnz_comb = 0;
  coeffd = VEC_get_data(coeff);
  Hd_comb = MAT_get_data_array(c->H_combined);
  for (k = 0; k < c->H_array_size; k++) {
    Hd = MAT_get_data_array(MAT_array_get(c->H_array,k));
    if (ensure_psd)
      coeffk = 0;
    else
      coeffk = coeffd[k];
    for (m = 0; m < MAT_get_nnz(MAT_array_get(c->H_array,k)); m++) {
      Hd_comb[H_nnz_comb] = coeffk*Hd[m];
      H_nnz_comb++;
    }
  }
}

void CONSTR_del_matvec(Constr* c) {
  if (c) {

    // Mat and vec
    VEC_del(c->b);
    MAT_del(c->A);
    VEC_del(c->f);
    MAT_del(c->J);
    MAT_del(c->Jbar);
    MAT_del(c->G);
    MAT_del(c->Gbar);
    VEC_del(c->l);
    VEC_del(c->u);
    MAT_array_del(c->H_array,c->H_array_size);
    MAT_del(c->H_combined);
    c->b = NULL;
    c->A = NULL;
    c->f = NULL;
    c->J = NULL;
    c->Jbar = NULL;
    c->G = NULL;
    c->Gbar = NULL;
    c->l = NULL;
    c->u = NULL;
    c->H_array = NULL;
    c->H_array_size = 0;
    c->H_combined = NULL;
  }
}

void CONSTR_del(Constr* c) {
  if (c) {

    // Mat and vec
    CONSTR_del_matvec(c);

    // Utils
    if (c->bus_counted)
      free(c->bus_counted);
    if (c->H_nnz)
      free(c->H_nnz);

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

char* CONSTR_get_type_str(Constr* c) {
  if (c)
    return c->type_str;
  else
    return "unknown";
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

Vec* CONSTR_get_l(Constr* c) {
  if (c)
    return c->l;
  else
    return NULL;
}

Vec* CONSTR_get_u(Constr* c) {
  if (c)
    return c->u;
  else
    return NULL;
}

Mat* CONSTR_get_G(Constr* c) {
  if (c)
    return c->G;
  else
    return NULL;
}

Mat* CONSTR_get_Gbar(Constr* c) {
  if (c)
    return c->Gbar;
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

Mat* CONSTR_get_Jbar(Constr* c) {
  if (c)
    return c->Jbar;
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

int CONSTR_get_A_nnz(Constr* c) {
  if (c)
    return c->A_nnz;
  else
    return 0;
}

int* CONSTR_get_A_nnz_ptr(Constr* c) {
  if (c)
    return &(c->A_nnz);
  else
    return NULL;
}

int CONSTR_get_G_nnz(Constr* c) {
  if (c)
    return c->G_nnz;
  else
    return 0;
}

int CONSTR_get_Gbar_nnz(Constr* c) {
  if (c)
    return c->Gbar_nnz;
  else
    return 0;
}

int* CONSTR_get_G_nnz_ptr(Constr* c) {
  if (c)
    return &(c->G_nnz);
  else
    return NULL;
}

int* CONSTR_get_Gbar_nnz_ptr(Constr* c) {
  if (c)
    return &(c->Gbar_nnz);
  else
    return NULL;
}

int CONSTR_get_J_nnz(Constr* c) {
  if (c)
    return c->J_nnz;
  else
    return 0;
}

int CONSTR_get_Jbar_nnz(Constr* c) {
  if (c)
    return c->Jbar_nnz;
  else
    return 0;
}

int* CONSTR_get_J_nnz_ptr(Constr* c) {
  if (c)
    return &(c->J_nnz);
  else
    return 0;
}

int* CONSTR_get_Jbar_nnz_ptr(Constr* c) {
  if (c)
    return &(c->Jbar_nnz);
  else
    return 0;
}

int* CONSTR_get_H_nnz(Constr* c) {
  if (c)
    return c->H_nnz;
  else
    return NULL;
}

int CONSTR_get_H_nnz_size(Constr* c) {
  if (c)
    return c->H_nnz_size;
  else
    return 0;
}

int CONSTR_get_A_row(Constr* c) {
  if (c)
    return c->A_row;
  else
    return 0;
}

int* CONSTR_get_A_row_ptr(Constr* c) {
  if (c)
    return &(c->A_row);
  else
    return NULL;
}

int CONSTR_get_G_row(Constr* c) {
  if (c)
    return c->G_row;
  else
    return 0;
}

int* CONSTR_get_G_row_ptr(Constr* c) {
  if (c)
    return &(c->G_row);
  else
    return NULL;
}

int CONSTR_get_J_row(Constr* c) {
  if (c)
    return c->J_row;
  else
    return 0;
}


int* CONSTR_get_J_row_ptr(Constr* c) {
  if (c)
    return &(c->J_row);
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

Constr* CONSTR_list_add(Constr* clist, Constr* nc) {
  LIST_add(Constr,clist,nc,next);
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

void CONSTR_list_count_step(Constr* clist, Branch* br, int t) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_count_step(cc,br,t);
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

void CONSTR_list_analyze_step(Constr* clist, Branch* br, int t) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_analyze_step(cc,br,t);
}

void CONSTR_list_eval_step(Constr* clist, Branch* br, int t, Vec* values) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_eval_step(cc,br,t,values);
}

void CONSTR_list_store_sens_step(Constr* clist, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  Constr* cc;
  Vec* vA;
  Vec* vf;
  Vec* vGu;
  Vec* vGl;
  int size_sA = 0;
  int size_sf = 0;
  int size_sG = 0;
  int offset_sA = 0;
  int offset_sf = 0;
  int offset_sG = 0;
  REAL* sA_data = VEC_get_data(sA);
  REAL* sf_data = VEC_get_data(sf);
  REAL* sGu_data = VEC_get_data(sGu);
  REAL* sGl_data = VEC_get_data(sGl);
  
  // Sizes
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    size_sA += MAT_get_size1(CONSTR_get_A(cc));
    size_sf += VEC_get_size(CONSTR_get_f(cc));
    size_sG += MAT_get_size1(CONSTR_get_G(cc));
  }
  
  // Map
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
  
    // Ax = b
    if (offset_sA + MAT_get_size1(CONSTR_get_A(cc)) <= VEC_get_size(sA))
      vA = VEC_new_from_array(&(sA_data[offset_sA]),MAT_get_size1(CONSTR_get_A(cc)));
    else
      vA = NULL;

    // f(x) = 0
    if (offset_sf + VEC_get_size(CONSTR_get_f(cc)) <= VEC_get_size(sf))
      vf = VEC_new_from_array(&(sf_data[offset_sf]),VEC_get_size(CONSTR_get_f(cc)));
    else
      vf = NULL;

    // Gx <= u
    if (offset_sG + MAT_get_size1(CONSTR_get_G(cc)) <= VEC_get_size(sGu))
      vGu = VEC_new_from_array(&(sGu_data[offset_sG]),MAT_get_size1(CONSTR_get_G(cc)));
    else
      vGu = NULL;

    // l <= Gx
    if (offset_sG + MAT_get_size1(CONSTR_get_G(cc)) <= VEC_get_size(sGl))
      vGl = VEC_new_from_array(&(sGl_data[offset_sG]),MAT_get_size1(CONSTR_get_G(cc)));
    else
      vGl = NULL;

    CONSTR_store_sens_step(cc,br,t,vA,vf,vGu,vGl);

    offset_sA += MAT_get_size1(CONSTR_get_A(cc));
    offset_sf += VEC_get_size(CONSTR_get_f(cc));
    offset_sG += MAT_get_size1(CONSTR_get_G(cc));
  }
}

Constr* CONSTR_new(int type, Net* net) {

  Constr* c = (Constr*)malloc(sizeof(Constr));

  // Error
  c->error_flag = FALSE;
  strcpy(c->error_string,"");

  // Network
  c->net = net;

  // Vars
  c->num_extra_vars = 0;

  // Fields
  c->type = type;
  strcpy(c->type_str,"");
  c->f = NULL;
  c->J = NULL;
  c->Jbar = NULL;
  c->H_array = NULL;  
  c->H_array_size = 0;
  c->H_combined = NULL;
  c->A = NULL;
  c->b = NULL;
  c->G = NULL;
  c->Gbar = NULL;
  c->l = NULL;
  c->u = NULL;
  c->A_nnz = 0;
  c->J_nnz = 0;
  c->Jbar_nnz = 0;
  c->G_nnz = 0;
  c->Gbar_nnz = 0;
  c->H_nnz = NULL;
  c->H_nnz_size = 0;
  c->A_row = 0;
  c->J_row = 0;
  c->G_row = 0;
  c->data = NULL;
  c->next = NULL;

  // Bus counted flags
  c->bus_counted_size = 0;
  c->bus_counted = NULL;

  // Constraint functions
  switch (type) {

    // AC power flow
  case CONSTR_TYPE_PF: 
    c->func_init = &CONSTR_PF_init;
    c->func_count_step = &CONSTR_PF_count_step;
    c->func_allocate = &CONSTR_PF_allocate;
    c->func_clear = &CONSTR_PF_clear;
    c->func_analyze_step = &CONSTR_PF_analyze_step;
    c->func_eval_step = &CONSTR_PF_eval_step;
    c->func_store_sens_step = &CONSTR_PF_store_sens_step;
    c->func_free = &CONSTR_PF_free;
    break;
  
    // DC power flow
  case CONSTR_TYPE_DCPF:
    c->func_init = &CONSTR_DCPF_init;
    c->func_count_step = &CONSTR_DCPF_count_step;
    c->func_allocate = &CONSTR_DCPF_allocate;
    c->func_clear = &CONSTR_DCPF_clear;
    c->func_analyze_step = &CONSTR_DCPF_analyze_step;
    c->func_eval_step = &CONSTR_DCPF_eval_step;
    c->func_store_sens_step = &CONSTR_DCPF_store_sens_step;
    c->func_free = &CONSTR_DCPF_free;
    break;
  
    // linearized power flow
  case CONSTR_TYPE_LINPF:
    c->func_init = &CONSTR_LINPF_init;
    c->func_count_step = &CONSTR_LINPF_count_step;
    c->func_allocate = &CONSTR_LINPF_allocate;
    c->func_clear = &CONSTR_LINPF_clear;
    c->func_analyze_step = &CONSTR_LINPF_analyze_step;
    c->func_eval_step = &CONSTR_LINPF_eval_step;
    c->func_store_sens_step = &CONSTR_LINPF_store_sens_step;
    c->func_free = &CONSTR_LINPF_free;
    break;
  
    // generator participation (active power)
  case CONSTR_TYPE_PAR_GEN_P:
    c->func_init = &CONSTR_PAR_GEN_P_init;
    c->func_count_step = &CONSTR_PAR_GEN_P_count_step;
    c->func_allocate = &CONSTR_PAR_GEN_P_allocate;
    c->func_clear = &CONSTR_PAR_GEN_P_clear;
    c->func_analyze_step = &CONSTR_PAR_GEN_P_analyze_step;
    c->func_eval_step = &CONSTR_PAR_GEN_P_eval_step;
    c->func_store_sens_step = &CONSTR_PAR_GEN_P_store_sens_step;
    c->func_free = &CONSTR_PAR_GEN_P_free;
    break;

    // generator participation (reactive power)
  case CONSTR_TYPE_PAR_GEN_Q:
    c->func_init = &CONSTR_PAR_GEN_Q_init;
    c->func_count_step = &CONSTR_PAR_GEN_Q_count_step;
    c->func_allocate = &CONSTR_PAR_GEN_Q_allocate;
    c->func_clear = &CONSTR_PAR_GEN_Q_clear;
    c->func_analyze_step = &CONSTR_PAR_GEN_Q_analyze_step;
    c->func_eval_step = &CONSTR_PAR_GEN_Q_eval_step;
    c->func_store_sens_step = &CONSTR_PAR_GEN_Q_store_sens_step;
    c->func_free = &CONSTR_PAR_GEN_Q_free;
    break;

    // variable fixing
  case CONSTR_TYPE_FIX:
    c->func_init = &CONSTR_FIX_init;
    c->func_count_step = &CONSTR_FIX_count_step;
    c->func_allocate = &CONSTR_FIX_allocate;
    c->func_clear = &CONSTR_FIX_clear;
    c->func_analyze_step = &CONSTR_FIX_analyze_step;
    c->func_eval_step = &CONSTR_FIX_eval_step;
    c->func_store_sens_step = &CONSTR_FIX_store_sens_step;
    c->func_free = &CONSTR_FIX_free;
    break;
    
    // voltage regulation by generator
  case CONSTR_TYPE_REG_GEN:
    c->func_init = &CONSTR_REG_GEN_init;
    c->func_count_step = &CONSTR_REG_GEN_count_step;
    c->func_allocate = &CONSTR_REG_GEN_allocate;
    c->func_clear = &CONSTR_REG_GEN_clear;
    c->func_analyze_step = &CONSTR_REG_GEN_analyze_step;
    c->func_eval_step = &CONSTR_REG_GEN_eval_step;
    c->func_store_sens_step = &CONSTR_REG_GEN_store_sens_step;
    c->func_free = &CONSTR_REG_GEN_free;
    break;
    
    // voltage regulation by transformer
  case CONSTR_TYPE_REG_TRAN:
    c->func_init = &CONSTR_REG_TRAN_init;
    c->func_count_step = &CONSTR_REG_TRAN_count_step;
    c->func_allocate = &CONSTR_REG_TRAN_allocate;
    c->func_clear = &CONSTR_REG_TRAN_clear;
    c->func_analyze_step = &CONSTR_REG_TRAN_analyze_step;
    c->func_eval_step = &CONSTR_REG_TRAN_eval_step;
    c->func_store_sens_step = &CONSTR_REG_TRAN_store_sens_step;
    c->func_free = &CONSTR_REG_TRAN_free;
    break;

    // voltage regulation by shunt device
  case CONSTR_TYPE_REG_SHUNT:
    c->func_init = &CONSTR_REG_SHUNT_init;
    c->func_count_step = &CONSTR_REG_SHUNT_count_step;
    c->func_allocate = &CONSTR_REG_SHUNT_allocate;
    c->func_clear = &CONSTR_REG_SHUNT_clear;
    c->func_analyze_step = &CONSTR_REG_SHUNT_analyze_step;
    c->func_eval_step = &CONSTR_REG_SHUNT_eval_step;
    c->func_store_sens_step = &CONSTR_REG_SHUNT_store_sens_step;
    c->func_free = &CONSTR_REG_SHUNT_free;
    break;
  
    // variable bounds (nonlin eq)
  case CONSTR_TYPE_BOUND:
    c->func_init = &CONSTR_BOUND_init;
    c->func_count_step = &CONSTR_BOUND_count_step;
    c->func_allocate = &CONSTR_BOUND_allocate;
    c->func_clear = &CONSTR_BOUND_clear;
    c->func_analyze_step = &CONSTR_BOUND_analyze_step;
    c->func_eval_step = &CONSTR_BOUND_eval_step;
    c->func_store_sens_step = &CONSTR_BOUND_store_sens_step;
    c->func_free = &CONSTR_BOUND_free;
    break;

    // DC branch flow limits
  case CONSTR_TYPE_DC_FLOW_LIM:
    c->func_init = &CONSTR_DC_FLOW_LIM_init;
    c->func_count_step = &CONSTR_DC_FLOW_LIM_count_step;
    c->func_allocate = &CONSTR_DC_FLOW_LIM_allocate;
    c->func_clear = &CONSTR_DC_FLOW_LIM_clear;
    c->func_analyze_step = &CONSTR_DC_FLOW_LIM_analyze_step;
    c->func_eval_step = &CONSTR_DC_FLOW_LIM_eval_step;
    c->func_store_sens_step = &CONSTR_DC_FLOW_LIM_store_sens_step;
    c->func_free = &CONSTR_DC_FLOW_LIM_free;
    break;

    // AC branch flow limits
  case CONSTR_TYPE_AC_FLOW_LIM:
    c->func_init = &CONSTR_AC_FLOW_LIM_init;
    c->func_count_step = &CONSTR_AC_FLOW_LIM_count_step;
    c->func_allocate = &CONSTR_AC_FLOW_LIM_allocate;
    c->func_clear = &CONSTR_AC_FLOW_LIM_clear;
    c->func_analyze_step = &CONSTR_AC_FLOW_LIM_analyze_step;
    c->func_eval_step = &CONSTR_AC_FLOW_LIM_eval_step;
    c->func_store_sens_step = &CONSTR_AC_FLOW_LIM_store_sens_step;
    c->func_free = &CONSTR_AC_FLOW_LIM_free;
    break;

    // variable bounds (lin ineq)
  case CONSTR_TYPE_LBOUND:
    c->func_init = &CONSTR_LBOUND_init;
    c->func_count_step = &CONSTR_LBOUND_count_step;
    c->func_allocate = &CONSTR_LBOUND_allocate;
    c->func_clear = &CONSTR_LBOUND_clear;
    c->func_analyze_step = &CONSTR_LBOUND_analyze_step;
    c->func_eval_step = &CONSTR_LBOUND_eval_step;
    c->func_store_sens_step = &CONSTR_LBOUND_store_sens_step;
    c->func_free = &CONSTR_LBOUND_free;
    break;

    // generator ramp constraints (lin ineq)
  case CONSTR_TYPE_GEN_RAMP:
    c->func_init = &CONSTR_GEN_RAMP_init;
    c->func_count_step = &CONSTR_GEN_RAMP_count_step;
    c->func_allocate = &CONSTR_GEN_RAMP_allocate;
    c->func_clear = &CONSTR_GEN_RAMP_clear;
    c->func_analyze_step = &CONSTR_GEN_RAMP_analyze_step;
    c->func_eval_step = &CONSTR_GEN_RAMP_eval_step;
    c->func_store_sens_step = &CONSTR_GEN_RAMP_store_sens_step;
    c->func_free = &CONSTR_GEN_RAMP_free;
    break;
    
    // unknown
  default: 
    c->func_init = NULL;
    c->func_count_step = NULL;
    c->func_allocate = NULL;
    c->func_clear = NULL;
    c->func_analyze_step = NULL;
    c->func_eval_step = NULL;
    c->func_store_sens_step = NULL;
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

void CONSTR_set_l(Constr* c, Vec* l) {
  if (c)
    c->l = l;
}

void CONSTR_set_u(Constr* c, Vec* u) {
  if (c)
    c->u = u;
}

void CONSTR_set_G(Constr* c, Mat* G) {
  if (c)
    c->G = G;
}

void CONSTR_set_Gbar(Constr* c, Mat* Gbar) {
  if (c)
    c->Gbar = Gbar;
}

void CONSTR_set_f(Constr* c, Vec* f) {
  if (c)
    c->f = f;
}

void CONSTR_set_J(Constr* c, Mat* J) {
  if (c)
    c->J = J;
}

void CONSTR_set_Jbar(Constr* c, Mat* Jbar) {
  if (c)
    c->Jbar = Jbar;
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

void CONSTR_set_A_nnz(Constr* c, int nnz) {
  if (c)
    c->A_nnz = nnz;
}

void CONSTR_set_G_nnz(Constr* c, int nnz) {
  if (c)
    c->G_nnz = nnz;
}

void CONSTR_set_Gbar_nnz(Constr* c, int nnz) {
  if (c)
    c->Gbar_nnz = nnz;
}

void CONSTR_set_J_nnz(Constr* c, int nnz) {
  if (c)
    c->J_nnz = nnz;
}

void CONSTR_set_Jbar_nnz(Constr* c, int nnz) {
  if (c)
    c->Jbar_nnz = nnz;
}

void CONSTR_set_H_nnz(Constr* c, int* nnz, int size) {
  if (c) {
    if (c->H_nnz)
      free(c->H_nnz);
    c->H_nnz = nnz;
    c->H_nnz_size = size;
  }
}

void CONSTR_set_A_row(Constr* c, int index) {
  if (c)
    c->A_row = index;
}

void CONSTR_set_G_row(Constr* c, int index) {
  if (c)
    c->G_row = index;
}

void CONSTR_set_J_row(Constr* c, int index) {
  if (c)
    c->J_row = index;
}

void CONSTR_set_bus_counted(Constr* c, char* counted, int size) {
  if (c) {
    if (c->bus_counted)
      free(c->bus_counted);
    c->bus_counted = counted;
    c->bus_counted_size = size;
  }
}

void CONSTR_set_data(Constr* c, void* data) {
  if (c)
    c->data = data;
}

void CONSTR_count(Constr* c) {
  int i;
  int t;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      CONSTR_count_step(c,NET_get_branch(net,i),t);
  }
}

void CONSTR_count_step(Constr* c, Branch* br, int t) {
  if (c && c->func_count_step && CONSTR_is_safe_to_count(c))
    (*(c->func_count_step))(c,br,t);
}

void CONSTR_allocate(Constr* c) {
  if (c && c->func_allocate && CONSTR_is_safe_to_count(c)) {
    CONSTR_del_matvec(c);
    (*(c->func_allocate))(c);
  }
}

void CONSTR_clear(Constr* c) {
  if (c && c->func_clear)
    (*(c->func_clear))(c);
}

void CONSTR_analyze(Constr* c) {
  int i;
  int t;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      CONSTR_analyze_step(c,NET_get_branch(net,i),t);
  }
}

void CONSTR_analyze_step(Constr* c, Branch* br, int t) {
  if (c && c->func_analyze_step && CONSTR_is_safe_to_analyze(c))
    (*(c->func_analyze_step))(c,br,t);
}

void CONSTR_eval(Constr* c, Vec* values) {
  int i;
  int t;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      CONSTR_eval_step(c,NET_get_branch(net,i),t,values);
  }
}

void CONSTR_eval_step(Constr* c, Branch* br, int t, Vec* values) {
  if (c && c->func_eval_step && CONSTR_is_safe_to_eval(c,values))
    (*(c->func_eval_step))(c,br,t,values);
}

void CONSTR_store_sens(Constr* c, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  int i;
  int t;
  Net* net = CONSTR_get_network(c);

  // No c
  if (!c)
    return;

  // Check sizes
  if ((VEC_get_size(sA) != MAT_get_size1(c->A)) ||
      (VEC_get_size(sf) != MAT_get_size1(c->J)) ||
      (VEC_get_size(sGu) != MAT_get_size1(c->G)) ||
      (VEC_get_size(sGl) != MAT_get_size1(c->G))) {
    sprintf(c->error_string,"invalid vector size");
    c->error_flag = TRUE;
    return;
  }

  // Clear
  CONSTR_clear(c);

  // Store sensitivities
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      CONSTR_store_sens_step(c,NET_get_branch(net,i),t,sA,sf,sGu,sGl);
  }
}

void CONSTR_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  if (c && c->func_store_sens_step && CONSTR_is_safe_to_count(c))
    (*(c->func_store_sens_step))(c,br,t,sA,sf,sGu,sGl);
}

BOOL CONSTR_is_safe_to_count(Constr* c) {
  Net* net = CONSTR_get_network(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net)*NET_get_num_periods(net))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to count");
    c->error_flag = TRUE;
    return FALSE;
  }  
}

BOOL CONSTR_is_safe_to_analyze(Constr* c) {
  Net* net = CONSTR_get_network(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
      MAT_get_size2(c->A) == NET_get_num_vars(net) &&
      MAT_get_size2(c->J) == NET_get_num_vars(net) &&
      MAT_get_size2(c->Jbar) == CONSTR_get_num_extra_vars(c))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to analyze");
    c->error_flag = TRUE;
    return FALSE;
  }
}

BOOL CONSTR_is_safe_to_eval(Constr* c, Vec* values) {
  Net* net = CONSTR_get_network(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
      MAT_get_size2(c->A) == NET_get_num_vars(net) &&
      MAT_get_size2(c->J) == NET_get_num_vars(net) &&
      MAT_get_size2(c->Jbar) == CONSTR_get_num_extra_vars(c) &&
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
  
  // No c
  if (!c)
    return;

  // Bus counted
  if (c->bus_counted)
    free(c->bus_counted);
  c->bus_counted_size = NET_get_num_buses(c->net)*NET_get_num_periods(c->net);
  ARRAY_zalloc(c->bus_counted,char,c->bus_counted_size);

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

