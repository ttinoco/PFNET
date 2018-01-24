/** @file constr_AC_LIN_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_LIN_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_AC_LIN_FLOW_LIM.h>

#if !HAVE_LINE_FLOW

// Dummies
void LF_set_branch_parameters(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
	double g, double b, double b_sh, double t_ratio, double t_shift, double I_max, LF_Branch *branch){

}

LF_Results* LF_construct(LF_Branch* branch, int flow_side, LF_Options *options) {
  return NULL;
}

void LF_free_results(LF_Results *results) {

}

double* LF_get_A_matrix(LF_Results *results) {
  return NULL;
}

double* LF_get_b_vector(LF_Results *results) {
  return NULL;
}

LF_ResultFlag LF_get_flag(LF_Results *results) {
  return error_other;
}

int LF_get_number_constraints(LF_Results *results) {
  return 0;
}

double LF_get_error(LF_Results *results) {
  return 0.;
}

char* LF_get_message(LF_Results *results) {
  return "LINE_FLOW library not available";
}

#endif

// Constraint data
struct Constr_AC_LIN_FLOW_LIM_Data {
  LF_Results** results; // array of pointers to LF results (1 per branch and time)
  int size;
};

// Constraint routines

Constr* CONSTR_AC_LIN_FLOW_LIM_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_AC_LIN_FLOW_LIM_init);
  CONSTR_set_func_count_step(c,&CONSTR_AC_LIN_FLOW_LIM_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_AC_LIN_FLOW_LIM_allocate);
  CONSTR_set_func_clear(c,&CONSTR_AC_LIN_FLOW_LIM_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_AC_LIN_FLOW_LIM_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_AC_LIN_FLOW_LIM_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_AC_LIN_FLOW_LIM_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_AC_LIN_FLOW_LIM_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_AC_LIN_FLOW_LIM_init(Constr* c) {
    
  // Local variables
  int i;
  int num;
  Net* net;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Init
  net = CONSTR_get_network(c);
  num = NET_get_num_branches(net)*NET_get_num_periods(net); // max number of LINE_FLOW results
  data = (Constr_AC_LIN_FLOW_LIM_Data*)malloc(sizeof(Constr_AC_LIN_FLOW_LIM_Data));
  data->results = (LF_Results**)malloc(sizeof(LF_Results*)*num);
  data->size = num;
  for (i = 0; i < num; i++)
    data->results[i] = NULL;
  CONSTR_set_name(c,"linearized AC branch flow limits");
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_AC_LIN_FLOW_LIM_clear(Constr* c) {
  
  // Counters
  CONSTR_set_G_nnz(c,0); // number of nonzeros of G matrix (from l <= Gx <= u)
  CONSTR_set_G_row(c,0); // number of rows of G matrix (from l <= Gx <= u)
}

void CONSTR_AC_LIN_FLOW_LIM_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  int k;
  int* G_nnz;
  int* G_row;
  Net* net;
  Bus* bus[2];
  REAL V_min[2];
  REAL V_max[2];
  int offset;
  int num_constr;
  LF_Results* results;
  LF_Branch branch;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Constr data
  net = CONSTR_get_network(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  data = (Constr_AC_LIN_FLOW_LIM_Data*)CONSTR_get_data(c);

  // Offset
  offset = NET_get_num_branches(net)*t;
  
  // Check pointer
  if (!G_nnz || !G_row || !data) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint has undefined pointers");  
    return;
  }

  // Outage (skip)
  if (BRANCH_is_on_outage(br))
    return;

  // Zero limits (skip)
  if (BRANCH_get_ratingA(br) == 0.)
    return;
  
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);

  // Check tap ratio and phase shift
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) ||
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint does not support tap ratios or phase shifts as variables");
    return;
  }

  // Check voltage angles
  if (!BUS_has_flags(bus[0],FLAG_VARS,BUS_VAR_VANG) &&
      !BUS_has_flags(bus[1],FLAG_VARS,BUS_VAR_VANG)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires at least one voltage angle as variable accross branch");
    return;
  }
  
  // Voltage magnitude limits
  for (k = 0; k < 2; k++) {
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) {
      if (!BUS_has_flags(bus[k],FLAG_BOUNDED,BUS_VAR_VMAG)) {
	CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires variable voltage magnitudes to be bounded");
	return;
      }
      V_min[k] = BUS_get_v_min_norm(bus[k]);
      V_max[k] = BUS_get_v_max_norm(bus[k]);
    }
    else {
      V_min[k] = BUS_get_v_mag(bus[k],t);
      V_max[k] = BUS_get_v_mag(bus[k],t);
    }
  }

  // Set parameters of branch structure of Line_Flow library
  LF_set_branch_parameters(V_min[0], 
			   V_max[0], 
			   V_min[1], 
			   V_max[1],
                           BRANCH_get_g(br), 
			   BRANCH_get_b(br),
			   BRANCH_get_b_k(br)+BRANCH_get_b_m(br), // total line charging
			   1./BRANCH_get_ratio(br,t),
			   BRANCH_get_phase(br,t)*180./PI, // degrees
			   BRANCH_get_ratingA(br),
                           &branch);

  // Construct linear constraints Ikm <= Imax
  results = LF_construct(&branch, 3, NULL);

  // Save results
  if (data->results[BRANCH_get_index(br)+offset])
    LF_free_results(data->results[BRANCH_get_index(br)+offset]); // delete existing
  data->results[BRANCH_get_index(br)+offset] = results; // save new
 
  // Check success
  if (LF_get_flag(results) != success && LF_get_flag(results) != non_binding) {
    CONSTR_set_error(c,LF_get_message(results));
    return;
  }
   
  // Count rows and nnz
  if (LF_get_flag(results) == success) {
    num_constr = LF_get_number_constraints(results);
    for (k = 0; k < 2; k++) {
      if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) // v_mag variable
	(*G_nnz) += num_constr;
      if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) // v_ang variable
	(*G_nnz) += num_constr;
    }
    (*G_row) += num_constr;
  }
}

void CONSTR_AC_LIN_FLOW_LIM_allocate(Constr* c) {

  // Local variables
  Net* net;
  int num_vars;
  int G_nnz;
  int G_row;

  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);
  G_nnz = CONSTR_get_G_nnz(c);
  G_row = CONSTR_get_G_row(c);

  // J f (from f(x) = 0)
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // A b (from Ax = b)
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // G l u (from l <= Gx <= u)
  CONSTR_set_G(c,MAT_new(G_row,    // size1 (rows)
			 num_vars, // size2 (cols)
			 G_nnz));  // nnz
  CONSTR_set_l(c,VEC_new(G_row));
  CONSTR_set_u(c,VEC_new(G_row));
}

void CONSTR_AC_LIN_FLOW_LIM_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Mat* G;
  Vec* l;
  Vec* u;
  REAL* A;
  REAL* b;
  int j;
  int k;
  int* G_nnz;
  int* G_row;
  Net* net;
  Bus* bus[2];
  int offset;
  int num_constr;
  LF_Results* results;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Constr data
  net = CONSTR_get_network(c);
  G = CONSTR_get_G(c);
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  data = (Constr_AC_LIN_FLOW_LIM_Data*)CONSTR_get_data(c);

  // Offset
  offset = NET_get_num_branches(net)*t;
  
  // Check pointer
  if (!G_nnz || !G_row || !data) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint has undefined pointers");  
    return;
  }

  // Outage (skip)
  if (BRANCH_is_on_outage(br))
    return;

  // Zero limits (skip)
  if (BRANCH_get_ratingA(br) == 0.)
    return;
  
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);

  // Check tap ratio and phase shift
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) ||
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint does not support tap ratios or phase shifts as variables");
    return;
  }

  // Check voltage angles
  if (!BUS_has_flags(bus[0],FLAG_VARS,BUS_VAR_VANG) &&
      !BUS_has_flags(bus[1],FLAG_VARS,BUS_VAR_VANG)) {
    CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires at least one voltage angle as variable accross branch");
    return;
  }
  
  // Check voltage magnitudes
  for (k = 0; k < 2; k++) {
    if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) {
      if (!BUS_has_flags(bus[k],FLAG_BOUNDED,BUS_VAR_VMAG)) {
	CONSTR_set_error(c,"AC_LIN_FLOW_LIM constraint requires variable voltage magnitudes to be bounded");
	return;
      }
    }
  }


  // Get results
  results = data->results[BRANCH_get_index(br)+offset];
  
  // Construct l <= Gx <= u
  if (LF_get_flag(results) == success) {

    num_constr = LF_get_number_constraints(results);
    A = LF_get_A_matrix(results);
    b = LF_get_b_vector(results);
    
    for (j = 0; j < num_constr; j++) {
  
      // l and u
      VEC_set(l,*G_row,-CONSTR_AC_LIN_FLOW_LIM_INF); // lower limit (minus infinity)
      VEC_set(u,*G_row,b[j]);                        // upper limit

      for (k = 0; k < 2; k++) {
  
	if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) { // v_mag variable

	  // G
	  MAT_set_i(G,*G_nnz,*G_row);
	  MAT_set_j(G,*G_nnz,BUS_get_index_v_mag(bus[k],t)); // vk time t
	  MAT_set_d(G,*G_nnz,A[3*j+k]);                      // k = 0 gives Vi, k = 1 gives Vj
	  (*G_nnz)++;
	}
	  
	if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) { // v_ang variable
	    
	  // G
	  MAT_set_i(G,*G_nnz,*G_row);
	  MAT_set_j(G,*G_nnz,BUS_get_index_v_ang(bus[k],t)); // wk time t
	  if (k == 0) 
	    MAT_set_d(G,*G_nnz,A[3*j+2]);                    // theta_i
	  else
	    MAT_set_d(G,*G_nnz,-A[3*j+2]);                   // theta_j
	  (*G_nnz)++;
	}
	else {
	  if (k == 0)
	    VEC_add_to_entry(u,*G_row,-A[3*j+2]*BUS_get_v_ang(bus[k],t));
	  else
	    VEC_add_to_entry(u,*G_row,A[3*j+2]*BUS_get_v_ang(bus[k],t));
	}
      }
      (*G_row)++;
    }
  }
}

void CONSTR_AC_LIN_FLOW_LIM_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing
}

void CONSTR_AC_LIN_FLOW_LIM_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing yet
}

void CONSTR_AC_LIN_FLOW_LIM_free(Constr* c) {

  // Local variables
  int i;
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Get data
  data = (Constr_AC_LIN_FLOW_LIM_Data*)CONSTR_get_data(c);

  // Free
  if (data) {
    for (i = 0; i < data->size; i++) {
      if (data->results[i])
	LF_free_results(data->results[i]);
    }
    free(data->results);
    free(data);
  }
   
  // Set data
  CONSTR_set_data(c,NULL);
}
