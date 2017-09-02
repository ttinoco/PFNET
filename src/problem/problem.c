/** @file problem.c
 *  @brief This file defines the Prob data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/problem.h>

struct Prob {

  // Error
  BOOL error_flag;                     /**< @brief Error flags */
  char error_string[PROB_BUFFER_SIZE]; /**< @brief Error string */

  // Output
  char output_string[PROB_BUFFER_SIZE]; /**< @brief Output string */

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

  // Linear equality constraints (Ax = 0)
  Vec* b;         /**< @brief Combined right-hand side of linear constraints */
  Mat* A;         /**< @brief Combined linear constraint matrix */
  
  // Nonlinear equality constraints (f(x) = 0)
  Vec* f;           /**< @brief Nonlinear constraint function values */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints */
  Mat* H_combined;  /**< @brief Combined Hessians of nonlinear constraints */

  // Linear inequality constraints (l <= Gx <= u)
  Mat* G;           /** @brief Matrix of constraint normals of linear inequality constraints */
  Vec* l;           /** @brief Lower bound for linear inequality contraints */
  Vec* u;           /** @brief Upper bound for linear inequality contraints */

  // Extra variables
  int num_extra_vars;          /** @brief Number of extra variables */
};

void PROB_add_constr(Prob* p, Constr* c) {
  if (p) {
    if (PROB_get_network(p) != CONSTR_get_network(c)) {
      sprintf(p->error_string,"constraint is associated with different network");
      p->error_flag = TRUE;
      return;
    }
    if (!PROB_find_constr(p,CONSTR_get_name(c)))
      p->constr = CONSTR_list_add(p->constr,c);
  }
}

void PROB_add_func(Prob* p, Func* f) {
  if (p) {
    if (PROB_get_network(p) != FUNC_get_network(f)) {
      sprintf(p->error_string,"function is associated with different network");
      p->error_flag = TRUE;
      return;
    }
    p->func = FUNC_list_add(p->func,f);
  }
}

void PROB_add_heur(Prob* p, int type) {
  if (p)
    p->heur = HEUR_list_add(p->heur,HEUR_new(type,p->net));
}

void PROB_analyze(Prob* p) {

  // Local variables
  Branch* br;
  Constr* c;
  Func* f;
  int Arow;
  int Annz;
  int Grow;
  int Gnnz;
  int Jrow;
  int Jnnz;
  int Hphinnz;
  int Hcombnnz;
  int num_vars;
  int num_extra_vars;
  int k;
  int t;
  
  // No p
  if (!p)
    return;

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  
  // Count
  for (t = 0; t < NET_get_num_periods(p->net); t++) {
    for (k = 0; k < NET_get_num_branches(p->net); k++) {
      
      br = NET_get_branch(p->net,k);
      
      // Constraints
      CONSTR_list_count_step(p->constr,br,t);
      if (CONSTR_list_has_error(p->constr)) {
	strcpy(p->error_string,CONSTR_list_get_error_string(p->constr));
	p->error_flag = TRUE;
	return;
      }
      
      // Functions
      FUNC_list_count_step(p->func,br,t);
      if (FUNC_list_has_error(p->func)) {
	strcpy(p->error_string,FUNC_list_get_error_string(p->func));
	p->error_flag = TRUE;
	return;
      }
    }
  }

  // Extra vars
  num_extra_vars = 0;
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c))
    num_extra_vars += CONSTR_get_num_extra_vars(c);
  p->num_extra_vars = num_extra_vars;

  // Allocate
  CONSTR_list_allocate(p->constr);
  FUNC_list_allocate(p->func);

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);

  // Analyze
  for (t = 0; t < NET_get_num_periods(p->net); t++) {
    for (k = 0; k < NET_get_num_branches(p->net); k++) {
      
      br = NET_get_branch(p->net,k);
      
      // Constraints
      CONSTR_list_analyze_step(p->constr,br,t);
      if (CONSTR_list_has_error(p->constr)) {
	strcpy(p->error_string,CONSTR_list_get_error_string(p->constr));
	p->error_flag = TRUE;
	return;
      }

      // Functions
      FUNC_list_analyze_step(p->func,br,t);
      if (FUNC_list_has_error(p->func)) {
	strcpy(p->error_string,FUNC_list_get_error_string(p->func));
	p->error_flag = TRUE;
	return;
      }
    }
  }
  CONSTR_list_finalize_structure_of_Hessians(p->constr);
  FUNC_list_finalize_structure_of_Hessian(p->func);

  // Delete matvec
  PROB_del_matvec(p);

  // Allocate matvec
  Arow = 0;
  Annz = 0;
  Grow = 0;
  Gnnz = 0;
  Jrow = 0;
  Jnnz = 0;
  Hphinnz = 0;
  Hcombnnz = 0;
  num_vars = NET_get_num_vars(p->net);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    Arow += MAT_get_size1(CONSTR_get_A(c));
    Annz += MAT_get_nnz(CONSTR_get_A(c));

    Grow += MAT_get_size1(CONSTR_get_G(c));
    Gnnz += MAT_get_nnz(CONSTR_get_G(c)); 

    Jrow += MAT_get_size1(CONSTR_get_J(c));
    Jnnz += MAT_get_nnz(CONSTR_get_J(c));
    Hcombnnz += MAT_get_nnz(CONSTR_get_H_combined(c));
  }
  for (f = p->func; f != NULL; f = FUNC_get_next(f))
    Hphinnz += MAT_get_nnz(FUNC_get_Hphi(f));

  p->phi = 0;
  p->gphi = VEC_new(num_vars+num_extra_vars);
  p->Hphi = MAT_new(num_vars+num_extra_vars,num_vars+num_extra_vars,Hphinnz);

  p->b = VEC_new(Arow);
  p->A = MAT_new(Arow,num_vars+num_extra_vars,Annz);

  p->l = VEC_new(Grow);
  p->u = VEC_new(Grow);
  p->G = MAT_new(Grow,num_vars+num_extra_vars,Gnnz);
  
  p->f = VEC_new(Jrow);
  p->J = MAT_new(Jrow,num_vars+num_extra_vars,Jnnz);
  p->H_combined = MAT_new(num_vars+num_extra_vars,num_vars+num_extra_vars,Hcombnnz);

  // Update
  PROB_update_lin(p);
  PROB_update_nonlin_struc(p);
}

void PROB_apply_heuristics(Prob* p, Vec* point) {

  // Local variables
  Branch* br;
  int i;
  int t;
  
  // No p
  if (!p)
    return;

  // Clear
  HEUR_list_clear(p->heur,p->net);

  // Apply
  for (t = 0; t < NET_get_num_periods(p->net); t++) {
    for (i = 0; i < NET_get_num_branches(p->net); i++) {
      br = NET_get_branch(p->net,i);
      HEUR_list_apply_step(p->heur,p->constr,p->net,br,t,point);
    }
  }
  
  // Udpate A and b
  PROB_update_lin(p);
}

void PROB_clear_error(Prob* p) {
  if (p) {

    // Problem
    p->error_flag = FALSE;
    strcpy(p->error_string,"");

    // Constraints
    CONSTR_list_clear_error(p->constr);

    // Functions
    FUNC_list_clear_error(p->func);
    
    // Network
    NET_clear_error(p->net);
  }
}

void PROB_eval(Prob* p, Vec* point) {

  // Local variables
  REAL* point_data;
  Branch* br;
  int num_vars;
  Vec* x;
  Vec* y;
  int k;
  int t;
  
  // No p
  if (!p)
    return;

  // Check sizes
  if (PROB_get_num_primal_variables(p) != VEC_get_size(point)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }

  // Extract x (network) and y (extra)
  point_data = VEC_get_data(point);
  num_vars = NET_get_num_vars(p->net);
  x = VEC_new_from_array(&(point_data[0]),num_vars);
  y = VEC_new_from_array(&(point_data[num_vars]),VEC_get_size(point)-num_vars);
 
  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  NET_clear_properties(p->net);

  // Eval
  for (t = 0; t < NET_get_num_periods(p->net); t++) {
    for (k = 0; k < NET_get_num_branches(p->net); k++) {
      
      br = NET_get_branch(p->net,k);
      
      // Constraints
      CONSTR_list_eval_step(p->constr,br,t,x,y);
      if (CONSTR_list_has_error(p->constr)) {
      	strcpy(p->error_string,CONSTR_list_get_error_string(p->constr));
	p->error_flag = TRUE;
	return;
      }
      
      // Functions
      FUNC_list_eval_step(p->func,br,t,x);
      if (FUNC_list_has_error(p->func)) {
      	strcpy(p->error_string,FUNC_list_get_error_string(p->func));
	p->error_flag = TRUE;
	return;
      }
      
      // Network
      NET_update_properties_step(p->net,br,t,x);
      if (NET_has_error(p->net)) {
	strcpy(p->error_string,NET_get_error_string(p->net));
	p->error_flag = TRUE;
	return;
      }
    }
  }

  // Clear
  free(x);
  free(y);

  // Update
  PROB_update_nonlin_data(p,point);
}

void PROB_store_sens(Prob* p, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Branch* br;
  int i;
  int t;
  
  // No p
  if (!p)
    return;

  // Check sizes
  if ((VEC_get_size(sA) != MAT_get_size1(p->A)) ||
      (VEC_get_size(sf) != MAT_get_size1(p->J)) ||
      (VEC_get_size(sGu) != MAT_get_size1(p->G)) ||
      (VEC_get_size(sGl) != MAT_get_size1(p->G))) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }

  // Clear
  CONSTR_list_clear(p->constr);

  // Store sens
  for (t = 0; t < NET_get_num_periods(p->net); t++) {
    for (i = 0; i < NET_get_num_branches(p->net); i++) {

      br = NET_get_branch(p->net,i);
      
      // Constraints
      CONSTR_list_store_sens_step(p->constr,br,t,sA,sf,sGu,sGl);
      if (CONSTR_list_has_error(p->constr)) {
	strcpy(p->error_string,CONSTR_list_get_error_string(p->constr));
	p->error_flag = TRUE;
	return;
      }
    }
  }
}

void PROB_del(Prob* p) {
  if (p) {    
    PROB_clear(p);
    free(p);
  }
}

void PROB_del_matvec(Prob* p) {
  if (p) {

    VEC_del(p->b);
    MAT_del(p->A);
    p->b = NULL;
    p->A = NULL;
    
    VEC_del(p->u);
    VEC_del(p->l);
    MAT_del(p->G);
    p->u = NULL;
    p->l = NULL;
    p->G = NULL;
 
    VEC_del(p->f);
    MAT_del(p->J);
    MAT_del(p->H_combined);
    p->f = NULL;
    p->J = NULL;
    p->H_combined = NULL;

    VEC_del(p->gphi);
    MAT_del(p->Hphi);
    p->gphi = NULL;
    p->Hphi = NULL;
  }
}

void PROB_clear(Prob* p) {
  if (p) {
    
    // Free constr, functions and heuristics
    CONSTR_list_del(p->constr);
    FUNC_list_del(p->func);
    HEUR_list_del(p->heur);
    
    // Free matvec
    PROB_del_matvec(p);

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
  int k;
  
  // Check inputs
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
    for (k = 0; k < MAT_get_nnz(CONSTR_get_H_combined(c)); k++) {
      Hcomb[Hcombnnz] = Hcomb_constr[k];
      Hcombnnz++;
    }
  }
}

Constr* PROB_find_constr(Prob* p, char* name) {
  Constr* cc;
  if (p) {
    for (cc = p->constr; cc != NULL; cc = CONSTR_get_next(cc)) {
      if (strcmp(name,CONSTR_get_name(cc)) == 0)
	return cc;
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
  if (p)
    return p->error_string;
  else
    return NULL;
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
  Vec* out;
  Vec* x;
  Vec* init_extra;
  Constr* c;
  int i;
  int j;
  if (!p)
    return NULL;
  out = VEC_new(PROB_get_num_primal_variables(p));
  x = NET_get_var_values(p->net,CURRENT);
  for (i = 0; i < VEC_get_size(x); i++)
    VEC_set(out,i,VEC_get(x,i));
  i = VEC_get_size(x);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    init_extra = CONSTR_get_init_extra_vars(c);
    for (j = 0; j < VEC_get_size(init_extra); j++) {
      VEC_set(out,i,VEC_get(init_extra,j));
      i++;
    }
  }
  free(x);  
  return out;
}

Vec* PROB_get_upper_limits(Prob* p) {
  Vec* out;
  Vec* x;
  Vec* u_extra;
  Constr* c;
  int i;
  int j;
  if (!p)
    return NULL;
  out = VEC_new(PROB_get_num_primal_variables(p));
  x = NET_get_var_values(p->net,UPPER_LIMITS);
  for (i = 0; i < VEC_get_size(x); i++)
    VEC_set(out,i,VEC_get(x,i));
  i = VEC_get_size(x);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    u_extra = CONSTR_get_u_extra_vars(c);
    for (j = 0; j < VEC_get_size(u_extra); j++) {
      VEC_set(out,i,VEC_get(u_extra,j));
      i++;
    }
  }
  free(x);
  return out;
}

Vec* PROB_get_lower_limits(Prob* p) {
  Vec* out;
  Vec* x;
  Vec* l_extra;
  Constr* c;
  int i;
  int j;
  if (!p)
    return NULL;
  out = VEC_new(PROB_get_num_primal_variables(p));
  x = NET_get_var_values(p->net,LOWER_LIMITS);
  for (i = 0; i < VEC_get_size(x); i++)
    VEC_set(out,i,VEC_get(x,i));
  i = VEC_get_size(x);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    l_extra = CONSTR_get_l_extra_vars(c);
    for (j = 0; j < VEC_get_size(l_extra); j++) {
      VEC_set(out,i,VEC_get(l_extra,j));
      i++;
    }
  }
  free(x);
  return out;
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

Mat* PROB_get_G(Prob* p) {
  if (p)
    return p->G;
  else 
    return NULL;
}

Vec* PROB_get_l(Prob* p) {
  if (p)
    return p->l;
  else 
    return NULL;
}

Vec* PROB_get_u(Prob* p) {
  if (p)
    return p->u;
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

int PROB_get_num_extra_vars(Prob* p) {
  if (p)
    return p->num_extra_vars;
  else
    return 0;
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
   
    // Output
    strcpy(p->output_string,"");
 
    p->constr = NULL;
    p->func = NULL;
    p->heur = NULL;
    
    p->phi = 0;
    p->gphi = NULL;
    p->Hphi = NULL;
    
    p->b = NULL;
    p->A = NULL;
    
    p->l = NULL;
    p->u = NULL;
    p->G = NULL;
    
    p->f = NULL;
    p->J = NULL;
    p->H_combined = NULL;

    p->num_extra_vars = 0;
  }
}

Prob* PROB_new(Net* net) {
  Prob* p = (Prob*)malloc(sizeof(Prob));
  p->net = net;
  PROB_init(p);
  return p;
}

char* PROB_get_show_str(Prob* p) {

  Func* f;
  Constr* c;
  char* out;

  if (!p)
    return NULL;

  out = p->output_string;
  strcpy(out,"");

  sprintf(out+strlen(out),"\nProblem\n");
  sprintf(out+strlen(out),"functions: %d\n",FUNC_list_len(p->func));
  for (f = p->func; f != NULL; f = FUNC_get_next(f))
    sprintf(out+strlen(out),"   %s\n",FUNC_get_name(f));
  sprintf(out+strlen(out),"constraints: %d\n",CONSTR_list_len(p->constr));  
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c))
    sprintf(out+strlen(out),"   %s\n",CONSTR_get_name(c));

  return out;
}

void PROB_show(Prob* p) {
  
  printf("%s",PROB_get_show_str(p));
}

void PROB_update_nonlin_struc(Prob* p) {
  /* This function fills in problem Jacobians and Hessians
     structure with constraint structure */

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

  int num_vars;
  int offset;
  int k;

  // Hphi
  Hphinnz = 0;
  Hphi_i = MAT_get_row_array(p->Hphi);
  Hphi_j = MAT_get_col_array(p->Hphi);
  for (f = p->func; f != NULL; f = FUNC_get_next(f)) {
    Hphi_i_func = MAT_get_row_array(FUNC_get_Hphi(f));
    Hphi_j_func = MAT_get_col_array(FUNC_get_Hphi(f));
    for (k = 0; k < MAT_get_nnz(FUNC_get_Hphi(f)); k++) {
      Hphi_i[Hphinnz] = Hphi_i_func[k];
      Hphi_j[Hphinnz] = Hphi_j_func[k];
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
  num_vars = NET_get_num_vars(p->net);
  offset = num_vars;
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // J
    Ji_constr = MAT_get_row_array(CONSTR_get_J(c));
    Jj_constr = MAT_get_col_array(CONSTR_get_J(c));
    for (k = 0; k < MAT_get_nnz(CONSTR_get_J(c)); k++) {
      Ji[Jnnz] = Jrow+Ji_constr[k];
      if (Jj_constr[k] < num_vars)
	Jj[Jnnz] = Jj_constr[k];                 // x var
      else
	Jj[Jnnz] = offset+Jj_constr[k]-num_vars; // y var
      Jnnz++;
    }
    Jrow += MAT_get_size1(CONSTR_get_J(c));

    // H comb
    Hcomb_i_constr = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hcomb_j_constr = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (k = 0; k < MAT_get_nnz(CONSTR_get_H_combined(c)); k++) {
      if (Hcomb_i_constr[k] < num_vars)
	Hcomb_i[Hcombnnz] = Hcomb_i_constr[k];                 // x var
      else
	Hcomb_i[Hcombnnz] = offset+Hcomb_i_constr[k]-num_vars; // y var
      if (Hcomb_j_constr[k] < num_vars)
	Hcomb_j[Hcombnnz] = Hcomb_j_constr[k];                 // x var
      else
	Hcomb_j[Hcombnnz] = offset+Hcomb_j_constr[k]-num_vars; // y var
      Hcombnnz++;
    }

    // Update offset
    offset += CONSTR_get_num_extra_vars(c);
  }
}

void PROB_update_nonlin_data(Prob* p, Vec* point) {
  /* This function fills in problem Jacobians and Hessians
     data with constraint data */
  
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
  int k;

  // Check sizes
  if (PROB_get_num_primal_variables(p) != VEC_get_size(point)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }
  
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
    for (k = 0; k < VEC_get_size(FUNC_get_gphi(func)); k++)
      gphi[k] += weight*gphi_func[k];

    // Hphi
    Hphi_func = MAT_get_data_array(FUNC_get_Hphi(func));
    for (k = 0; k < MAT_get_nnz(FUNC_get_Hphi(func)); k++) {
      Hphi[Hphinnz] = weight*Hphi_func[k];
      Hphinnz++;
    }     
  }

  // f and derivatives
  Jnnz = 0;
  Jrow = 0;
  f = VEC_get_data(p->f);
  J = MAT_get_data_array(p->J);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    
    // J, f of constraint
    f_constr = VEC_get_data(CONSTR_get_f(c));
    J_constr = MAT_get_data_array(CONSTR_get_J(c));

    // Update f
    for (k = 0; k < VEC_get_size(CONSTR_get_f(c)); k++) {
      f[Jrow+k] = f_constr[k];
    }

    // Update J 
    for (k = 0; k < MAT_get_nnz(CONSTR_get_J(c)); k++) {
      J[Jnnz] = J_constr[k];
      Jnnz++;
    }

    // Update row
    Jrow += MAT_get_size1(CONSTR_get_J(c));
  }
}

void PROB_update_lin(Prob* p) {
  /* This function updates problem A,b,G,l,u with 
     constraint A,b,G,l,u. */
  
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
  
  REAL* l;
  REAL* l_constr;
  REAL* u;
  REAL* u_constr;
  int* Gi;
  int* Gj;
  REAL* Gd;
  int* Gi_constr;
  int* Gj_constr;
  REAL* Gd_constr;
  int Gnnz;
  int Grow;

  int num_vars;
  int offset;
  int k;
  
  // Check
  if (!p)
    return;

  // Init problem data
  Annz = 0;
  Arow = 0;
  b = VEC_get_data(p->b);
  Ai = MAT_get_row_array(p->A);
  Aj = MAT_get_col_array(p->A);
  Ad = MAT_get_data_array(p->A);

  Gnnz = 0;
  Grow = 0;
  l = VEC_get_data(p->l);
  u = VEC_get_data(p->u);
  Gi = MAT_get_row_array(p->G);
  Gj = MAT_get_col_array(p->G);
  Gd = MAT_get_data_array(p->G);
  
  // Process constraints
  num_vars = NET_get_num_vars(p->net);
  offset = num_vars;
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // A and b of constraint
    b_constr = VEC_get_data(CONSTR_get_b(c));
    Ai_constr = MAT_get_row_array(CONSTR_get_A(c));
    Aj_constr = MAT_get_col_array(CONSTR_get_A(c));
    Ad_constr = MAT_get_data_array(CONSTR_get_A(c));

    // G, l, u of constraint
    l_constr = VEC_get_data(CONSTR_get_l(c));
    u_constr = VEC_get_data(CONSTR_get_u(c));
    Gi_constr = MAT_get_row_array(CONSTR_get_G(c));
    Gj_constr = MAT_get_col_array(CONSTR_get_G(c));
    Gd_constr = MAT_get_data_array(CONSTR_get_G(c));

    // Update A
    for (k = 0; k < MAT_get_nnz(CONSTR_get_A(c)); k++) {
      Ai[Annz] = Arow+Ai_constr[k];
      if (Aj_constr[k] < num_vars)
	Aj[Annz] = Aj_constr[k];                 // x var
      else
	Aj[Annz] = offset+Aj_constr[k]-num_vars; // y var
      Ad[Annz] = Ad_constr[k];
      Annz++;
    }

    // Update b
    for (k = 0; k < MAT_get_size1(CONSTR_get_A(c)); k++) {
      b[Arow] = b_constr[k];
      Arow++;
    }
    
    // Update G
    for (k = 0; k < MAT_get_nnz(CONSTR_get_G(c)); k++) {
      Gi[Gnnz] = Grow+Gi_constr[k];
      if (Gj_constr[k] < num_vars)
	Gj[Gnnz] = Gj_constr[k];                 // x var
      else
	Gj[Gnnz] = offset+Gj_constr[k]-num_vars; // y var
      Gd[Gnnz] = Gd_constr[k];
      Gnnz++;
    }

    // Update l,u
    for (k = 0; k < MAT_get_size1(CONSTR_get_G(c)); k++) {
      l[Grow] = l_constr[k];
      u[Grow] = u_constr[k];
      Grow++;
    }

    // Update offset
    offset += CONSTR_get_num_extra_vars(c);
  }
}

int PROB_get_num_primal_variables(Prob* p) {
  if (!p)
    return 0;
  if (p->gphi)
    return VEC_get_size(p->gphi);
  if (p->Hphi)
    return MAT_get_size1(p->Hphi);
  if (p->A)
    return MAT_get_size2(p->A);
  if (p->J)
    return MAT_get_size2(p->J);
  return NET_get_num_vars(p->net)+p->num_extra_vars;
}

int PROB_get_num_linear_equality_constraints(Prob* p) {
  if (!p)
    return 0;
  if (p->A)
    return MAT_get_size1(p->A);
  return 0;
}

int PROB_get_num_nonlinear_equality_constraints(Prob* p) {
  if (!p)
    return 0;
  if (p->f)
    return VEC_get_size(p->f);
  return 0;
}

