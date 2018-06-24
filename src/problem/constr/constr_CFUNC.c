/** @file constr_CFUNC.c
 *  @brief This file defines the data structure and routines associated with the constraint of type CFUNC.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_CFUNC.h>

struct Constr_CFUNC_Data {
  Func* func; 
  REAL rhs;
  char op[CONSTR_CFUNC_BUFFER_SIZE];
};

Constr* CONSTR_CFUNC_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_CFUNC_init);
  CONSTR_set_func_count_step(c, &CONSTR_CFUNC_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_CFUNC_analyze_step);
  CONSTR_set_func_allocate(c, &CONSTR_CFUNC_allocate);
  CONSTR_set_func_clear(c, &CONSTR_CFUNC_clear);
  CONSTR_set_func_eval_step(c, &CONSTR_CFUNC_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_CFUNC_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_CFUNC_free);
  CONSTR_set_func_set_parameter(c, &CONSTR_CFUNC_set_parameter);
  CONSTR_init(c);
  return c;
}

void CONSTR_CFUNC_init(Constr* c) {
  
  // Local variaables
  Constr_CFUNC_Data* data;
  
  data = (Constr_CFUNC_Data*)malloc(sizeof(Constr_CFUNC_Data));
  data->func = NULL;
  data->rhs = 0;
  strcpy(data->op,"=");
  CONSTR_set_H_nnz(c,(int*)calloc(1,sizeof(int)),1);
  CONSTR_set_name(c,"constrained function");
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_CFUNC_clear(Constr* c) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  
  // Check
  if (!data)
    return;

  // Data
  FUNC_clear(data->func);
}

void CONSTR_CFUNC_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  int num_vars;
  int num_extra_vars;
  int* Hnnz;
  
  // Check
  if (!data)
    return;

  // Function count step
  FUNC_count_step(data->func,br,t);

  // Post-processing
  if ((t == BRANCH_get_num_periods(br)-1) && (BRANCH_get_index(br) == NET_get_num_branches(net)-1)) {

    // Extra vars
    if (strcmp(data->op,">=") == 0 || strcmp(data->op,"<=") == 0)
      CONSTR_set_num_extra_vars(c, 1);
    else
      CONSTR_set_num_extra_vars(c, 0);

    num_vars = NET_get_num_vars(CONSTR_get_network(c));
    num_extra_vars = CONSTR_get_num_extra_vars(c);

    // A
    // nothing

    // G 
    CONSTR_set_G_row(c,num_extra_vars);
    CONSTR_set_G_nnz(c,num_extra_vars);

    // J
    CONSTR_set_J_row(c,1);
    CONSTR_set_J_nnz(c,num_vars+num_extra_vars); // dense

    // H
    H_nnz = CONSTR_get_H_nnz(c);
    H_nnz[0] = MAT_get_nnz(FUNC_get_Hphi(data->func));
  }
}

void CONSTR_CFUNC_allocate(Constr* c) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  
  // Check
  if (!data)
    return;

  // Function
  FUNC_allocate(data->func);
}

void CONSTR_CFUNC_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  Net* net = CONSTR_get_network(c);
  int num_vars = NET_get_num_vars(net);
  Mat* J;
  Mat* H;
  Mat* Hphi;
  Mat* G;
  Vec* l;
  Vec* u;
  int i;
  
  // Check
  if (!data)
    return;

  // Count step
  FUNC_analyze_step(data->func,br,t);
  
  // Post-processing
  if ((t == BRANCH_get_num_periods(br)-1) && (BRANCH_get_index(br) == NET_get_num_branches(net)-1)) {

    G = CONSTR_get_G(c);
    l = CONSTR_get_l(c);
    u = CONSTR_get_u(c);
    J = CONSTR_get_J(c);
    H = CONSTR_get_H_single(c,0);
    Hphi = FUNC_get_Hphi(data->func);
    
    // Inequality >=
    if (strcmp(data->op,">=") == 0) {

      // Extra vars
      VEC_set(CONSTR_get_l_extra_vars(c),0,0);
      VEC_set(CONSTR_get_u_extra_vars(c),0,CONSTR_CFUNC_EXTRA_VAR_INF);
      VEC_set(CONSTR_get_init_extra_vars(c),0,0);
      
      // G l u
      MAT_set_i(G,0,0);
      MAT_set_j(G,0,num_vars);
      MAT_set_d(G,0,1.);
      VEC_set(l,0,0);
      VEC_set(u,0,CONSTR_CFUNC_EXTRA_VAR_INF);

      // J extra vars
      MAT_set_i(J,num_vars,0);
      MAT_set_j(J,num_vars,num_vars);
    } 

    // Inequality <=
    else if (strcmp(data->op,"<=") == 0) {

      // Extra vars
      VEC_set(CONSTR_get_l_extra_vars(c),0,-CONSTR_CFUNC_EXTRA_VAR_INF);
      VEC_set(CONSTR_get_u_extra_vars(c),0,0);
      VEC_set(CONSTR_get_init_extra_vars(c),0,0);

      // G l u
      MAT_set_i(G,0,0);
      MAT_set_j(G,0,num_vars);
      MAT_set_d(G,0,1.);
      VEC_set(l,0,-CONSTR_CFUNC_EXTRA_VAR_INF);
      VEC_set(u,0,0);

      // J extra vars
      MAT_set_i(J,num_vars,0);
      MAT_set_j(J,num_vars,num_vars);
    }

    // J normal vars
    for (i = 0; i < num_vars; i++) {
      MAT_set_i(J,i,0);
      MAT_set_j(J,i,i);
    }

    // H
    for (i = 0; i < MAT_get_nnz(Hphi); i++) {
      MAT_set_i(H,i,MAT_get_i(Hphi,i));
      MAT_set_j(H,i,MAT_get_j(Hphi,i));
    }
  }
}

void CONSTR_CFUNC_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  
  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  Net* net = CONSTR_get_network(c);
  int num_vars = NET_get_num_vars(net);
  REAL extra_var = 0;
  REAL* f;
  REAL* J;
  REAL* H;
  int Hnnz;
  REAL* Hphi;
  REAL* gphi;
  int i;
  
  // Check
  if (!data)
    return;

  // Count step
  FUNC_eval_step(data->func,br,t,values);
  
  // Post-processing
  if ((t == BRANCH_get_num_periods(br)-1) && (BRANCH_get_index(br) == NET_get_num_branches(net)-1)) {

    f = VEC_get_data(CONSTR_get_f(c));
    J = MAT_get_data_array(CONSTR_get_J(c));
    H = MAT_get_data_array(CONSTR_get_H_single(c,0));
    Hnnz = MAT_get_nnz(CONSTR_get_H_single(c,0));
    gphi = VEC_get_data(FUNC_get_gphi(data->func));
    Hphi = MAT_get_data_array(FUNC_get_Hphi(data->func));
    
    // Inequality >= or <=
    if ((strcmp(data->op,">=") == 0) || (strcmp(data->op,"<=") == 0)) {

      // Extra var
      if (VEC_get_size(values_extra) > 0)
	extra_var = VEC_get(values_extra,0);
      else
	extra_var = 0;

      // J extra vars
      J[num_vars] = -1;
    }

    // f
    f[0] = FUNC_get_phi(data->func) - data->rhs - extra_var;
    
    // J normal vars
    for (i = 0; i < num_vars; i++)
      J[i] = gphi[i];

    // H
    for (i = 0; i < Hnnz; i++)
      H[i] = Hphi[i];
  }
}

void CONSTR_CFUNC_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_CFUNC_set_parameter(Constr* c, char* key, void* value) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  
  // Check
  if (!data)
    return;

  // Set 
  if (strcmp(key,"func") == 0) // function
    data->func = (Func*)value;
	
  else if (strcmp(key,"rhs") == 0) // right-hand-side
    data->rhs = *((REAL*)value);

  else if (strcmp(key,"op") == 0) // operator
    strncpy(data->op,(char*)value,CONSTR_CFUNC_BUFFER_SIZE);
  
  else // unknown
    CONSTR_set_error(c,"invalid parameter");  
}

void CONSTR_CFUNC_free(Constr* c) {

  // Local variables
  Constr_CFUNC_Data* data = (Constr_CFUNC_Data*)CONSTR_get_data(c);
  
  // Free
  if (data) {
    FUNC_del(data->func);
    free(data);
  }

  // Set data
  CONSTR_set_data(c,NULL);
}
