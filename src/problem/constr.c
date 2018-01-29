/** @file constr.c
 *  @brief This file defines the Constr data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr.h>

struct Constr {

  // Error
  BOOL error_flag;                       /**< @brief Error flag */
  char error_string[CONSTR_BUFFER_SIZE]; /**< @brief Error string */

  // Name
  char name[CONSTR_BUFFER_SIZE]; /**< @brief Name string */
  
  // Network
  Net* net;    /**< @brief Power network */
  
  // Nonlinear (f(x,y) = 0)
  Vec* f;           /**< @brief Vector of nonlinear constraint violations */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints */
  Mat* H_array;     /**< @brief Array of Hessian matrices of nonlinear constraints */
  int H_array_size; /**< @brief Size of Hessian array */
  Mat* H_combined;  /**< @brief Linear combination of Hessians of the nonlinear constraints */
  
  // Linear equality (A (x,y) = b)
  Mat* A;           /**< @brief Matrix of constraint normals of linear equality constraints */
  Vec* b;           /**< @brief Right-hand side vector of linear equality constraints */

  // Linear inequalities (l <= G (x,y) <= h)
  Mat* G;           /** @brief Matrix of constraint normals of linear inequality constraints wrt. variables */
  Vec* l;           /** @brief Lower bound for linear inequality contraints */
  Vec* u;           /** @brief Upper bound for linear inequality contraints */
  
  // Extra variables
  int num_extra_vars;   /** @brief Number of extra variables (set during "count") */
  Vec* l_extra_vars;    /** @brief Lower bounds for extra variables (set during "analyze") */
  Vec* u_extra_vars;    /** @brief Upper bounds for extra varaibles (set during "analyze") */
  Vec* init_extra_vars; /** @brief Extra variable initial values */
  
  // Counters and flags
  int A_nnz;             /**< @brief Counter for nonzeros of matrix A */
  int J_nnz;             /**< @brief Counter for nonzeros of matrix J */
  int G_nnz;             /**< @brief Counter for nonzeros of matrix G */
  int* H_nnz;            /**< @brief Array of counters of nonzeros of nonlinear constraint Hessians */
  int H_nnz_size;        /**< @brief Size of array of counter of Hessian nonzeros */
  int A_row;             /**< @brief Counter for linear equality constraints */
  int J_row;             /**< @brief Counter for nonlinear constraints */
  int G_row;             /**< @brief Counter for linear inequality constraints */
  char* bus_counted;     /**< @brief Flag for processing buses */
  int bus_counted_size;  /**< @brief Size of array of flags for processing buses */

  // Row info
  char* A_row_info; /**< @brief Array for info strings of rows of A (x,y) = b */
  char* J_row_info; /**< @brief Array for info strings of rows of f(x,y) = 0 */
  char* G_row_info; /**< @brief Array for info strings of rows of l <= G (x,y) <= u */
  
  // Type functions
  void (*func_init)(Constr* c);                                          /**< @brief Initialization function */
  void (*func_count_step)(Constr* c, Branch* br, int t);                 /**< @brief Function for counting nonzero entries */
  void (*func_allocate)(Constr* c);                                      /**< @brief Function for allocating required arrays */
  void (*func_clear)(Constr* c);                                         /**< @brief Function for clearing flags, counters, and function values */
  void (*func_analyze_step)(Constr* c, Branch* br, int t);               /**< @brief Function for analyzing sparsity pattern */
  void (*func_eval_step)(Constr* c, Branch* br, int t, Vec* v, Vec* ve); /**< @brief Function for evaluating constraint */
  void (*func_store_sens_step)(Constr* c, Branch* br, int t,
			       Vec* sA, Vec* sf, Vec* sGu, Vec* sGl);    /**< @brief Func. for storing sensitivities */
  void (*func_free)(Constr* c);                                          /**< @brief Function for de-allocating any data used */
  void (*func_set_parameter)(Constr* c, char* key, void* value);         /**< @brief Function for setting constraint parameter */

  // Type data
  void* data; /**< @brief Type-dependent constraint data structure */

  // List
  Constr* next; /**< @brief List of constraints */
};

void CONSTR_allocate_H_array(Constr* c, int size) {
  if (c) {
    if (c->H_array)
      MAT_array_del(c->H_array,c->H_array_size);
    c->H_array = MAT_array_new(size);
    c->H_array_size = size;
  }
}

void CONSTR_allocate_H_combined(Constr* c) {
  int i;
  int H_comb_nnz = 0;
  int num_vars = NET_get_num_vars(CONSTR_get_network(c));
  int num_extra_vars = CONSTR_get_num_extra_vars(c);
  if (c) {
    for (i = 0; i < c->H_array_size; i++)
      H_comb_nnz += MAT_get_nnz(MAT_array_get(c->H_array,i));
    if (c->H_combined)
      MAT_del(c->H_combined);
    CONSTR_set_H_combined(c,MAT_new(num_vars+num_extra_vars, // size1 (rows)
				    num_vars+num_extra_vars, // size2 (cols)
				    H_comb_nnz));
  }
}

void CONSTR_finalize_structure_of_Hessians(Constr* c) {

  // Local variables
  int* Hi_comb;
  int* Hj_comb;
  int H_nnz_comb;
  Mat* H_array;
  int* Hi;
  int* Hj;
  int temp;
  int k;
  int m;
  
  // Check
  if (!c)
    return;

  // Ensure lower triangular and save struct of H comb
  H_nnz_comb = 0;
  H_array = CONSTR_get_H_array(c);
  Hi_comb = MAT_get_row_array(CONSTR_get_H_combined(c));
  Hj_comb = MAT_get_col_array(CONSTR_get_H_combined(c));
  for (k = 0; k < CONSTR_get_H_array_size(c); k++) {
    Hi = MAT_get_row_array(MAT_array_get(H_array,k));
    Hj = MAT_get_col_array(MAT_array_get(H_array,k));
    for (m = 0; m < MAT_get_nnz(MAT_array_get(H_array,k)); m++) {
      if (Hi[m] < Hj[m]) {
	temp = Hi[m];
	Hi[m] = Hj[m];
	Hj[m] = temp;
      }
      Hi_comb[H_nnz_comb] = Hi[m];
      Hj_comb[H_nnz_comb] = Hj[m];
      H_nnz_comb++;
    }
  }
}

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
    MAT_del(c->G);
    VEC_del(c->l);
    VEC_del(c->u);
    VEC_del(c->l_extra_vars);
    VEC_del(c->u_extra_vars);
    VEC_del(c->init_extra_vars);
    MAT_array_del(c->H_array,c->H_array_size);
    MAT_del(c->H_combined);
    c->b = NULL;
    c->A = NULL;
    c->f = NULL;
    c->J = NULL;
    c->G = NULL;
    c->l = NULL;
    c->u = NULL;
    c->l_extra_vars = NULL;
    c->u_extra_vars = NULL;
    c->init_extra_vars = NULL;
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

    // Row infos
    if (c->A_row_info)
      free(c->A_row_info);
    if (c->J_row_info)
      free(c->J_row_info);
    if (c->G_row_info)
      free(c->G_row_info);

    // Data
    if (c->func_free)
      (*(c->func_free))(c);

    free(c);
  }
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

Vec* CONSTR_get_l_extra_vars(Constr* c) {
  if (c)
    return c->l_extra_vars;
  else
    return NULL;
}

Vec* CONSTR_get_u_extra_vars(Constr* c) {
  if (c)
    return c->u_extra_vars;
  else
    return NULL;
}

Vec* CONSTR_get_init_extra_vars(Constr* c) {
  if (c)
    return c->init_extra_vars;
  else
    return NULL;
}

Mat* CONSTR_get_G(Constr* c) {
  if (c)
    return c->G;
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

int* CONSTR_get_G_nnz_ptr(Constr* c) {
  if (c)
    return &(c->G_nnz);
  else
    return NULL;
}

int CONSTR_get_J_nnz(Constr* c) {
  if (c)
    return c->J_nnz;
  else
    return 0;
}

int* CONSTR_get_J_nnz_ptr(Constr* c) {
  if (c)
    return &(c->J_nnz);
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

Mat* CONSTR_get_var_projection(Constr* c) {

  // Local variables
  Mat* P;
  int i;

  // Check
  if (!c)
    return NULL;

  // Allocate
  P = MAT_new(NET_get_num_vars(c->net),
	      NET_get_num_vars(c->net)+c->num_extra_vars,
	      NET_get_num_vars(c->net));

  // Fill
  for (i = 0; i < MAT_get_nnz(P); i++) {
    MAT_set_i(P,i,i);
    MAT_set_j(P,i,i);
    MAT_set_d(P,i,1.);
  }

  // Return
  return P;
}

Mat* CONSTR_get_extra_var_projection(Constr* c) {

    // Local variables
  Mat* P;
  int i;

  // Check
  if (!c)
    return NULL;

  // Allocate
  P = MAT_new(c->num_extra_vars,
	      NET_get_num_vars(c->net)+c->num_extra_vars,
	      c->num_extra_vars);

  // Fill
  for (i = 0; i < MAT_get_nnz(P); i++) {
    MAT_set_i(P,i,i);
    MAT_set_j(P,i,NET_get_num_vars(c->net)+i);
    MAT_set_d(P,i,1.);
  }

  // Return
  return P;
}

char* CONSTR_get_A_row_info_string(Constr* c, int index) {
  if (c && c->A_row_info && 0 <= index && index < MAT_get_size1(c->A))
    return c->A_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char);
  else
    return "";
}

char* CONSTR_get_J_row_info_string(Constr* c, int index) {
  if (c && c->J_row_info && 0 <= index && index < MAT_get_size1(c->J))
    return c->J_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char);
  else
    return "";
}

char* CONSTR_get_G_row_info_string(Constr* c, int index) {
  if (c && c->G_row_info && 0 <= index && index < MAT_get_size1(c->G))
    return c->G_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char);
  else
    return "";
}

void CONSTR_list_finalize_structure_of_Hessians(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_finalize_structure_of_Hessians(cc);
}

void CONSTR_list_clear_error(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc))
    CONSTR_clear_error(cc);
}

BOOL CONSTR_list_has_error(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (CONSTR_has_error(cc))
      return TRUE;
  }
  return FALSE;
}

char* CONSTR_list_get_error_string(Constr* clist) {
  Constr* cc;
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (CONSTR_has_error(cc))
      return CONSTR_get_error_string(cc);
  }
  return "";
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
  int offset = 0;
  REAL* coeffd = VEC_get_data(coeff);
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (offset + VEC_get_size(CONSTR_get_f(cc)) <= VEC_get_size(coeff))
      v = VEC_new_from_array(&(coeffd[offset]),VEC_get_size(CONSTR_get_f(cc)));
    else
      v = NULL;       
    CONSTR_combine_H(cc,v,ensure_psd);
    offset += VEC_get_size(CONSTR_get_f(cc));
    if (v)
      free(v);
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

void CONSTR_list_eval_step(Constr* clist, Branch* br, int t, Vec* v, Vec* ve) {
  Constr* cc;
  Vec* ve_c;
  int offset = 0;
  REAL* ve_data = VEC_get_data(ve);
  for (cc = clist; cc != NULL; cc = CONSTR_get_next(cc)) {
    if (offset + CONSTR_get_num_extra_vars(cc) <= VEC_get_size(ve))
      ve_c = VEC_new_from_array(&(ve_data[offset]),CONSTR_get_num_extra_vars(cc));
    else
      ve_c = NULL;
    CONSTR_eval_step(cc,br,t,v,ve_c);
    offset += CONSTR_get_num_extra_vars(cc);
    if (ve_c)
      free(ve_c);
  }
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

    if (vA)
      free(vA);
    if (vf)
      free(vf);
    if (vGu)
      free(vGu);
    if (vGl)
      free(vGl);
  }
}

Constr* CONSTR_new(Net* net) {

  Constr* c = (Constr*)malloc(sizeof(Constr));

  // Error
  c->error_flag = FALSE;
  strcpy(c->error_string,"");

  // Name
  strcpy(c->name,"unknown");

  // Network
  c->net = net;

  // Vars
  c->num_extra_vars = 0;

  // Fields
  c->f = NULL;
  c->J = NULL;
  c->H_array = NULL;  
  c->H_array_size = 0;
  c->H_combined = NULL;
  c->A = NULL;
  c->b = NULL;
  c->G = NULL;
  c->l = NULL;
  c->u = NULL;
  c->l_extra_vars = NULL;
  c->u_extra_vars = NULL;
  c->init_extra_vars = NULL;
  c->A_nnz = 0;
  c->J_nnz = 0;
  c->G_nnz = 0;
  c->H_nnz = NULL;
  c->H_nnz_size = 0;
  c->A_row = 0;
  c->J_row = 0;
  c->G_row = 0;
  c->next = NULL;

  // Row infos
  c->A_row_info = NULL;
  c->J_row_info = NULL;
  c->G_row_info = NULL;

  // Bus counted flags
  c->bus_counted_size = 0;
  c->bus_counted = NULL;

  // Methods
  c->func_init = NULL;
  c->func_count_step = NULL;
  c->func_allocate = NULL;
  c->func_clear = NULL;
  c->func_analyze_step = NULL;
  c->func_eval_step = NULL;
  c->func_store_sens_step = NULL;
  c->func_free = NULL;
  c->func_set_parameter = NULL;
  
  // Data
  c->data = NULL;

  // Update network
  CONSTR_update_network(c);
  
  return c;
}

void CONSTR_set_parameter(Constr* c, char* key, void* value) {
  if (c) {
    if (c->func_set_parameter)
      (*(c->func_set_parameter))(c, key, value);
    else {
      sprintf(c->error_string,"constraint does not support setting parameters");
      c->error_flag = TRUE;
      return;
    }
  }
}

void CONSTR_set_name(Constr* c, char* name) {
  if (c)
    strcpy(c->name,name);
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

void CONSTR_set_l_extra_vars(Constr* c, Vec* l) {
  if (c)
    c->l_extra_vars = l;
}

void CONSTR_set_u_extra_vars(Constr* c, Vec* u) {
  if (c)
    c->u_extra_vars = u;
}

void CONSTR_set_init_extra_vars(Constr* c, Vec* init) {
  if (c)
    c->init_extra_vars = init;
}

void CONSTR_set_G(Constr* c, Mat* G) {
  if (c)
    c->G = G;
}

void CONSTR_set_f(Constr* c, Vec* f) {
  if (c) {
    if (c->f)
      VEC_del(c->f);
    c->f = f;
  }
}

void CONSTR_set_J(Constr* c, Mat* J) {
  if (c) {
    if (c->J)
      MAT_del(c->J);
    c->J = J;
  }
}

void CONSTR_set_H_array(Constr* c, Mat* array, int size) {
  if (c) {
    if (c->H_array)
      MAT_array_del(c->H_array,c->H_array_size);
    c->H_array = array;
    c->H_array_size = size;
  }  
}

void CONSTR_set_H_single(Constr* c, int i, Mat* m) {
  Mat* H;
  if (c && 0 <= i && i < c->H_array_size) {
    H = MAT_array_get(c->H_array,i);
    MAT_set_size1(H,MAT_get_size1(m));
    MAT_set_size2(H,MAT_get_size2(m));
    MAT_set_row_array(H,MAT_get_row_array(m));
    MAT_set_col_array(H,MAT_get_col_array(m));
    MAT_set_data_array(H,MAT_get_data_array(m));
    MAT_set_nnz(H,MAT_get_nnz(m));
    MAT_set_owns_rowcol(H,MAT_get_owns_rowcol(m));
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

void CONSTR_set_J_nnz(Constr* c, int nnz) {
  if (c)
    c->J_nnz = nnz;
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

void CONSTR_set_A_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time) {
  char info[CONSTR_INFO_BUFFER_SIZE];
  if (c && c->A_row_info && 0 <= index && index < MAT_get_size1(c->A)) {
    snprintf(info,
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s:%s:%d:%s:%d",
	     c->name,
	     obj,
	     obj_id,
	     constr_info,
	     time);
    snprintf(c->A_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s",
	     info);
  }
}

void CONSTR_set_J_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time) {
  char info[CONSTR_INFO_BUFFER_SIZE];
  if (c && c->J_row_info && 0 <= index && index < MAT_get_size1(c->J)) {
    snprintf(info,
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s:%s:%d:%s:%d",
	     c->name,
	     obj,
	     obj_id,
	     constr_info,
	     time);
    snprintf(c->J_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s",
	     info);
  }
}

void CONSTR_set_G_row_info_string(Constr* c, int index, char* obj, int obj_id, char* constr_info, int time) {
  char info[CONSTR_INFO_BUFFER_SIZE];
  if (c && c->G_row_info && 0 <= index && index < MAT_get_size1(c->G)) {
    snprintf(info,
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s:%s:%d:%s:%d",
	     c->name,
	     obj,
	     obj_id,
	     constr_info,
	     time);
    snprintf(c->G_row_info+index*CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     CONSTR_INFO_BUFFER_SIZE*sizeof(char),
	     "%s",
	     info);
  }
}

void CONSTR_init(Constr* c) {
  if (c && c->func_free)
    (*(c->func_free))(c);
  if (c && c->func_init)
    (*(c->func_init))(c);
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
    CONSTR_allocate_H_combined(c);
    c->A_row_info = (char*)malloc(sizeof(char)*CONSTR_INFO_BUFFER_SIZE*MAT_get_size1(c->A));
    c->J_row_info = (char*)malloc(sizeof(char)*CONSTR_INFO_BUFFER_SIZE*MAT_get_size1(c->J));
    c->G_row_info = (char*)malloc(sizeof(char)*CONSTR_INFO_BUFFER_SIZE*MAT_get_size1(c->G));
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
  CONSTR_finalize_structure_of_Hessians(c);
}

void CONSTR_analyze_step(Constr* c, Branch* br, int t) {
  if (c && c->func_analyze_step && CONSTR_is_safe_to_analyze(c))
    (*(c->func_analyze_step))(c,br,t);
}

void CONSTR_eval(Constr* c, Vec* v, Vec* ve) {
  int i;
  int t;
  Net* net = CONSTR_get_network(c);
  CONSTR_clear(c);
  for (t = 0; t < NET_get_num_periods(net); t++) {
    for (i = 0; i < NET_get_num_branches(net); i++)
      CONSTR_eval_step(c,NET_get_branch(net,i),t,v,ve);
  }
}

void CONSTR_eval_step(Constr* c, Branch* br, int t, Vec* v, Vec* ve) {
  if (c && c->func_eval_step && CONSTR_is_safe_to_eval(c,v,ve))
    (*(c->func_eval_step))(c,br,t,v,ve);
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
  NET_clear_sensitivities(c->net);

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
  int num_vars = NET_get_num_vars(net)+CONSTR_get_num_extra_vars(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
      MAT_get_size2(CONSTR_get_A(c)) == num_vars &&
      MAT_get_size2(CONSTR_get_G(c)) == num_vars &&
      MAT_get_size2(CONSTR_get_J(c)) == num_vars &&
      VEC_get_size(CONSTR_get_l_extra_vars(c)) == CONSTR_get_num_extra_vars(c) &&
      VEC_get_size(CONSTR_get_u_extra_vars(c)) == CONSTR_get_num_extra_vars(c) &&
      VEC_get_size(CONSTR_get_init_extra_vars(c)) == CONSTR_get_num_extra_vars(c))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to analyze");
    c->error_flag = TRUE;
    return FALSE;
  }
}

BOOL CONSTR_is_safe_to_eval(Constr* c, Vec* v, Vec* ve) {
  Net* net = CONSTR_get_network(c);
  int num_vars = NET_get_num_vars(net)+CONSTR_get_num_extra_vars(c);
  if (CONSTR_get_bus_counted_size(c) == NET_get_num_buses(net)*NET_get_num_periods(net) &&
      MAT_get_size2(CONSTR_get_A(c)) == num_vars &&
      MAT_get_size2(CONSTR_get_G(c)) == num_vars &&
      MAT_get_size2(CONSTR_get_J(c)) == num_vars &&
      VEC_get_size(v) == NET_get_num_vars(net) &&
      (VEC_get_size(ve) == CONSTR_get_num_extra_vars(c) || 
       VEC_get_size(ve) == 0))
    return TRUE;
  else {
    sprintf(c->error_string,"constraint is not safe to eval");
    c->error_flag = TRUE;
    return FALSE;
  }
}

BOOL CONSTR_has_error(Constr* c) {
  if (c)
    return c->error_flag;
  else
    return FALSE;
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
  if (c)
    return c->error_string;
  else
    return NULL;
}

char* CONSTR_get_name(Constr* c) {
  if (c)
    return c->name;
  else
    return NULL;
}

Net* CONSTR_get_network(Constr* c) {
  if (c)
    return c->net;
  else
    return NULL;
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

  // Init
  CONSTR_init(c);
}

void CONSTR_set_func_init(Constr* c, void (*func)(Constr* c)) {
  if (c)
    c->func_init = func;
}

void CONSTR_set_func_count_step(Constr* c, void (*func)(Constr* c, Branch* br, int t)) {
  if (c)
    c->func_count_step = func;
}

void CONSTR_set_func_allocate(Constr* c, void (*func)(Constr* c)) {
  if (c)
    c->func_allocate = func;
}

void CONSTR_set_func_clear(Constr* c, void (*func)(Constr* c)) {
  if (c)
    c->func_clear = func;
}

void CONSTR_set_func_analyze_step(Constr* c, void (*func)(Constr* c, Branch* br, int t)) {
  if (c)
    c->func_analyze_step = func;
}

void CONSTR_set_func_eval_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* v, Vec* ve)) {
  if (c)
    c->func_eval_step = func;
}

void CONSTR_set_func_store_sens_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl)) {
  if (c)
    c->func_store_sens_step = func;
}

void CONSTR_set_func_free(Constr* c, void (*func)(Constr* c)) {
  if (c)
    c->func_free = func;
}

void CONSTR_set_func_set_parameter(Constr* c, void (*func)(Constr* c, char* key, void* value)) {
  if (c)
    c->func_set_parameter = func;
}
