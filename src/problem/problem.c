/** @file problem.c
 *  @brief This file defines the Prob data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/problem.h>

struct Prob {

  // Error
  BOOL error_flag;                     /**< @brief Error flags */
  char error_string[PROB_BUFFER_SIZE]; /**< @brief Error string */

  // Constraints
  Constr* constr; /**< @brief List of constraints */

  // Functions
  Func* func;     /**< @brief List of objective functions */

  // Heurisics
  Heur* heur;     /**< @brief List of heuristics */
  
  // Network
  Net* net;       /**< @brief Power flow network */

  // Objective function
  REAL phi;       /**< @brief Combined objective value */
  Vec* gphi;      /**< @brief Gradient of combined objective function */
  Mat* Hphi;      /**< @brief Hessian of combined objective function */

  // Linear equality constraints
  Vec* b;         /**< @brief Combined right-hand side of linear constraints */
  Mat* A;         /**< @brief Combined linear constraint matrix */
  
  // Nonlinear equality constraints
  Vec* f;           /**< @brief Nonlinear constraint function values */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints */
  Mat* H_combined;  /**< @brief Combined Hessians of nonlinear constraints */
};

void PROB_add_constr(Prob* p, int type) {
  if (p) {
    if (!PROB_find_constr(p,type))
      p->constr = CONSTR_list_add(p->constr,CONSTR_new(type,p->net));
  }
}

void PROB_add_func(Prob* p, int type, REAL weight) {
  if (p)
    p->func = FUNC_list_add(p->func,FUNC_new(type,weight,p->net));
}

void PROB_add_heur(Prob* p, int type) {
  if (p)
    p->heur = HEUR_list_add(p->heur,
			    HEUR_new(type,p->net));
}

void PROB_analyze(Prob* p) {

  // Local variables
  Branch* br;
  Constr* c;
  Func* f;
  int Arow;
  int Annz;
  int Jrow;
  int Jnnz;
  int Hphinnz;
  int Hcombnnz;
  int num_vars;
  int i;
  
  if (!p)
    return;

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  
  // Count
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_count_branch(p->constr,br);
    FUNC_list_count_branch(p->func,br);
  }
    
  // Allocate
  CONSTR_list_allocate(p->constr);
  FUNC_list_allocate(p->func);

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);

  // Analyze
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_analyze_branch(p->constr,br);
    FUNC_list_analyze_branch(p->func,br);
  }

  // Allocate problem matvec
  Annz = 0;
  Arow = 0;
  Jrow = 0;
  Jnnz = 0;
  Hphinnz = 0;
  Hcombnnz = 0;
  num_vars = NET_get_num_vars(p->net);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    Arow += MAT_get_size1(CONSTR_get_A(c));
    Annz += MAT_get_nnz(CONSTR_get_A(c));
    Jrow += MAT_get_size1(CONSTR_get_J(c));
    Jnnz += MAT_get_nnz(CONSTR_get_J(c));
    Hcombnnz += MAT_get_nnz(CONSTR_get_H_combined(c));
  }
  for (f = p->func; f != NULL; f = FUNC_get_next(f))
    Hphinnz += MAT_get_nnz(FUNC_get_Hphi(f));
  p->phi = 0;
  p->gphi = VEC_new(num_vars);
  p->Hphi = MAT_new(num_vars,num_vars,Hphinnz);
  p->b = VEC_new(Arow);
  p->A = MAT_new(Arow,num_vars,Annz);
  p->f = VEC_new(Jrow);
  p->J = MAT_new(Jrow,num_vars,Jnnz);
  p->H_combined = MAT_new(num_vars,num_vars,Hcombnnz);

  // Update 
  PROB_update_lin(p);
  PROB_update_nonlin_struc(p);
}

void PROB_apply_heuristics(Prob* p, Vec* point) {

  // Local variables
  Branch* br;
  Constr* c;
  int i;
  
  if (!p)
    return;

  // Clear
  HEUR_list_clear(p->heur,p->net);

  // Apply
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    HEUR_list_apply_to_branch(p->heur,p->constr,p->net,br,point);
  }

  // Udpate A and b
  PROB_update_lin(p);
}

void PROB_eval(Prob* p, Vec* point) {

  // Local variables
  Branch* br;
  int i;
  
  if (!p)
    return;

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  NET_clear_properties(p->net);

  // Eval
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_eval_branch(p->constr,br,point);
    FUNC_list_eval_branch(p->func,br,point);
    NET_update_properties_branch(p->net,br,point);
  }

  // Update 
  PROB_update_nonlin_data(p);
}

void PROB_store_sens(Prob* p, Vec* sens) {

  // Local variables
  Branch* br;
  int i;
  
  if (!p || !sens)
    return;

  // Check size
  if (VEC_get_size(sens) != VEC_get_size(p->f)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }

  // Clear
  CONSTR_list_clear(p->constr);

  // Store sens
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_store_sens_branch(p->constr,br,sens);
  }
}

void PROB_del(Prob* p) {
  if (p) {
    PROB_clear(p);
    free(p);
    p = NULL;
  }
}

void PROB_clear(Prob* p) {
  if (p) {
    
    // Free data
    CONSTR_list_del(p->constr);
    FUNC_list_del(p->func);
    HEUR_list_del(p->heur);
    VEC_del(p->b);
    MAT_del(p->A);
    VEC_del(p->f);
    MAT_del(p->J);
    MAT_del(p->H_combined);
    VEC_del(p->gphi);
    MAT_del(p->Hphi);    

    // Re-initialize
    PROB_init(p);
  }
}

void PROB_combine_H(Prob* p, Vec* coeff, BOOL ensure_psd) {
  
  // Local variables
  Constr* c;
  int Hcombnnz;
  REAL* Hcomb;
  REAL* Hcomb_constr;
  int i;
  
  if (!p || !coeff)
    return;

  // Check size
  if (VEC_get_size(coeff) != VEC_get_size(p->f)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }
  
  // Combine
  CONSTR_list_combine_H(p->constr,coeff,ensure_psd);

  // Combine and update
  Hcombnnz = 0;
  Hcomb = MAT_get_data_array(p->H_combined);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    Hcomb_constr = MAT_get_data_array(CONSTR_get_H_combined(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_H_combined(c)); i++) {
      Hcomb[Hcombnnz] = Hcomb_constr[i];
      Hcombnnz++;
    }
  }
}

Constr* PROB_find_constr(Prob* p, int constr_type) {
  Constr* c;
  if (p) {
    for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
      if (CONSTR_get_type(c) == constr_type)
	return c;
    }
    return NULL;
  }
  else
    return NULL;
}

Constr* PROB_get_constr(Prob* p) {
  if (p)
    return p->constr;
  else
    return NULL;
}

char* PROB_get_error_string(Prob* p) {
  if (!p)
    return NULL;
  else
    return p->error_string;
}

Func* PROB_get_func(Prob* p) {
  if (p)
    return p->func;
  else
    return NULL;
}

Heur* PROB_get_heur(Prob* p) {
  if (p)
    return p->heur;
  else
    return NULL;
}

Vec* PROB_get_init_point(Prob* p) {
  if (p)
    return NET_get_var_values(p->net,CURRENT);
  else
    return NULL;
}

Vec* PROB_get_upper_limits(Prob* p) {
  if (p)
    return NET_get_var_values(p->net,UPPER_LIMITS);
  else
    return NULL;
}

Vec* PROB_get_lower_limits(Prob* p) {
  if (p)
    return NET_get_var_values(p->net,LOWER_LIMITS);
  else
    return NULL;
}

Net* PROB_get_network(Prob* p) {
  if (p)
    return p->net;
  else
    return NULL;
}

REAL PROB_get_phi(Prob* p) {
  if (p)
    return p->phi;
  else
    return 0;
}

Vec* PROB_get_gphi(Prob* p) {
  if (p)
    return p->gphi;
  else
    return NULL;
}

Mat* PROB_get_Hphi(Prob* p) {
  if (p)
    return p->Hphi;
  else
    return NULL;
}

Mat* PROB_get_A(Prob* p) {
  if (p)
    return p->A;
  else 
    return NULL;
}

Vec* PROB_get_b(Prob* p) {
  if (p)
    return p->b;
  else 
    return NULL;
}

Mat* PROB_get_J(Prob* p) {
  if (p)
    return p->J;
  else 
    return NULL;
}

Vec* PROB_get_f(Prob* p) {
  if (p)
    return p->f;
  else 
    return NULL;
}

Mat* PROB_get_H_combined(Prob* p) {
  if (p)
    return p->H_combined;
  else
    return NULL;
}

BOOL PROB_has_error(Prob* p) {
  if (!p)
    return FALSE;
  else
    return p->error_flag;
}

void PROB_init(Prob* p) {
  if (p) {

    // Error
    p->error_flag = FALSE;
    strcpy(p->error_string,"");
    
    p->constr = NULL;
    p->func = NULL;
    p->heur = NULL;
    p->net = NULL;
    p->phi = 0;
    p->gphi = NULL;
    p->Hphi = NULL;
    p->b = NULL;
    p->A = NULL;
    p->f = NULL;
    p->J = NULL;
    p->H_combined = NULL;
  }
}

Prob* PROB_new(void) {
  Prob* p = (Prob*)malloc(sizeof(Prob));
  PROB_init(p);
  return p;
}

void PROB_set_network(Prob* p, Net* net) {
  if (p)
    p->net = net;
}

void PROB_show(Prob* p) {
  Func* f;
  if (p) {
    printf("\nProblem\n");
    printf("functions  : %d\n",FUNC_list_len(p->func));
    for (f = p->func; f != NULL; f = FUNC_get_next(f)) {
      printf("\ttype %d \t value %.4e\n",FUNC_get_type(f),FUNC_get_phi(f));
    }
    printf("constraints: %d\n",CONSTR_list_len(p->constr));  
    printf("heuristics : %d\n",HEUR_list_len(p->heur));
  }
}

void PROB_update_nonlin_struc(Prob* p) {
  
  // Local variables
  Func* f;
  Constr* c;
  int* Hphi_i;
  int* Hphi_j;
  int* Hphi_i_func;
  int* Hphi_j_func;
  int Hphinnz;
  int* Ji;
  int* Jj;
  int* Ji_constr;
  int* Jj_constr;
  int Jrow;
  int Jnnz;
  int* Hcomb_i;
  int* Hcomb_j;
  int* Hcomb_i_constr;
  int* Hcomb_j_constr;
  int Hcombnnz;
  int i;

  // Hphi
  Hphinnz = 0;
  Hphi_i = MAT_get_row_array(p->Hphi);
  Hphi_j = MAT_get_col_array(p->Hphi);
  for (f = p->func; f != NULL; f = FUNC_get_next(f)) {
    Hphi_i_func = MAT_get_row_array(FUNC_get_Hphi(f));
    Hphi_j_func = MAT_get_col_array(FUNC_get_Hphi(f));
    for (i = 0; i < MAT_get_nnz(FUNC_get_Hphi(f)); i++) {
      Hphi_i[Hphinnz] = Hphi_i_func[i];
      Hphi_j[Hphinnz] = Hphi_j_func[i];
      Hphinnz++;
    }     
  }

  // J and Hcomb
  Jnnz = 0;
  Jrow = 0;
  Ji = MAT_get_row_array(p->J);
  Jj = MAT_get_col_array(p->J);
  Hcombnnz = 0;
  Hcomb_i = MAT_get_row_array(p->H_combined);
  Hcomb_j = MAT_get_col_array(p->H_combined);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // J
    Ji_constr = MAT_get_row_array(CONSTR_get_J(c));
    Jj_constr = MAT_get_col_array(CONSTR_get_J(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_J(c)); i++) {
      Ji[Jnnz] = Ji_constr[i]+Jrow;
      Jj[Jnnz] = Jj_constr[i];
      Jnnz++;
    }
    Jrow += MAT_get_size1(CONSTR_get_J(c));

    // H comb
    Hcomb_i_constr = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hcomb_j_constr = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_H_combined(c)); i++) {
      Hcomb_i[Hcombnnz] = Hcomb_i_constr[i];
      Hcomb_j[Hcombnnz] = Hcomb_j_constr[i];
      Hcombnnz++;
    }
  }
}

void PROB_update_nonlin_data(Prob* p) {
  
  // Local variables
  Func* func;
  Constr* c;
  REAL* gphi;
  REAL* gphi_func;
  REAL* Hphi;
  REAL* Hphi_func;
  REAL weight;
  int Hphinnz;
  REAL* f;
  REAL* f_constr;
  REAL* J;
  REAL* J_constr;
  int Jnnz;
  int Jrow;
  int i;
  
  // phi and derivatives
  p->phi = 0;
  Hphinnz = 0;
  VEC_set_zero(p->gphi);
  gphi = VEC_get_data(p->gphi);
  Hphi = MAT_get_data_array(p->Hphi);
  for (func = p->func; func != NULL; func = FUNC_get_next(func)) {

    // Weight
    weight = FUNC_get_weight(func);

    // phi
    p->phi += weight*FUNC_get_phi(func);

    //gphi
    gphi_func = VEC_get_data(FUNC_get_gphi(func));
    for (i = 0; i < NET_get_num_vars(p->net); i++)
      gphi[i] += weight*gphi_func[i];

    // Hphi
    Hphi_func = MAT_get_data_array(FUNC_get_Hphi(func));
    for (i = 0; i < MAT_get_nnz(FUNC_get_Hphi(func)); i++) {
      Hphi[Hphinnz] = weight*Hphi_func[i];
      Hphinnz++;
    }     
  }

  // f and derivatives
  Jnnz = 0;
  Jrow = 0;
  f = VEC_get_data(p->f);
  J = MAT_get_data_array(p->J);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // J and f of constraint
    f_constr = VEC_get_data(CONSTR_get_f(c));
    J_constr = MAT_get_data_array(CONSTR_get_J(c));

    // Update J
    for (i = 0; i < MAT_get_nnz(CONSTR_get_J(c)); i++) {
      J[Jnnz] = J_constr[i];
      Jnnz++;
    }

    // Update f
    for (i = 0; i < MAT_get_size1(CONSTR_get_J(c)); i++) {
      f[Jrow] = f_constr[i];
      Jrow++;
    }
  }
}

void PROB_update_lin(Prob* p) {
  
  // Local variables
  Constr* c;
  REAL* b;
  REAL* b_constr;
  int* Ai;
  int* Aj;
  REAL* Ad;
  int* Ai_constr;
  int* Aj_constr;
  REAL* Ad_constr;
  int Annz;
  int Arow;
  int i;
  
  if (!p)
    return;

  // Update
  Annz = 0;
  Arow = 0;
  b = VEC_get_data(p->b);
  Ai = MAT_get_row_array(p->A);
  Aj = MAT_get_col_array(p->A);
  Ad = MAT_get_data_array(p->A);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // A and b of constraint
    b_constr = VEC_get_data(CONSTR_get_b(c));
    Ai_constr = MAT_get_row_array(CONSTR_get_A(c));
    Aj_constr = MAT_get_col_array(CONSTR_get_A(c));
    Ad_constr = MAT_get_data_array(CONSTR_get_A(c));

    // Update A
    for (i = 0; i < MAT_get_nnz(CONSTR_get_A(c)); i++) {
      Ai[Annz] = Ai_constr[i]+Arow;
      Aj[Annz] = Aj_constr[i];
      Ad[Annz] = Ad_constr[i];
      Annz++;
    }

    // Update b
    for (i = 0; i < MAT_get_size1(CONSTR_get_A(c)); i++) {
      b[Arow] = b_constr[i];
      Arow++;
    }
  }
}

