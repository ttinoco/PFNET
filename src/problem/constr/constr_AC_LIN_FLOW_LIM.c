/** @file constr_AC_LIN_FLOW_LIM.c
 *  @brief This file defines the data structure and routines associated with the constraint of type AC_LIN_FLOW_LIM.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_AC_LIN_FLOW_LIM.h>

#if !HAVE_LINE_FLOW

// Dummies
LINE_FLOW_Results* LINE_FLOW_construct(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
				       double g, double b, double B_sh, double Kt_real, double Kt_shift, 
				       double I_max_user, int flow_side,
				       LINE_FLOW_Params* params) {
  return NULL;
}

void LINE_FLOW_free_results(LINE_FLOW_Results* results) {

}

double* LINE_FLOW_get_A_matrix(LINE_FLOW_Results* results) {
  return NULL;
}

double* LINE_FLOW_get_b_vector(LINE_FLOW_Results* results) {
  return NULL;
}

LINE_FLOW_Flag LINE_FLOW_get_flag(LINE_FLOW_Results* results) {
  return error_other;
}

int LINE_FLOW_get_number_constraints(LINE_FLOW_Results* results) {
  return 0;
}

double LINE_FLOW_get_error(LINE_FLOW_Results* results) {
  return 0.;
}

char* LINE_FLOW_get_message(LINE_FLOW_Results* results) {
  return "LINE_FLOW library not available";
}

#endif

// Constraint data
struct Constr_AC_LIN_FLOW_LIM_Data {
  LINE_FLOW_Results** results; // array of pointers to LINE_FLOW results (2 per branch and time)
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
  num = 2*NET_get_num_branches(net)*NET_get_num_periods(net); // max number of LINE_FLOW results
  data = (Constr_AC_LIN_FLOW_LIM_Data*)malloc(sizeof(Constr_AC_LIN_FLOW_LIM_Data));
  data->results = (LINE_FLOW_Results**)malloc(sizeof(LINE_FLOW_Results*)*num);
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
  int i;
  int k;
  int* G_nnz;
  int* G_row;
  Net* net;
  Bus* bus[2];
  REAL V_min[2];
  REAL V_max[2];
  int offset;
  int num_constr;
  LINE_FLOW_Results* results[2];
  Constr_AC_LIN_FLOW_LIM_Data* data;

  // Constr data
  net = CONSTR_get_network(c);
  G_nnz = CONSTR_get_G_nnz_ptr(c);
  G_row = CONSTR_get_G_row_ptr(c);
  data = (Constr_AC_LIN_FLOW_LIM_Data*)CONSTR_get_data(c);

  // Offset
  offset = 2*NET_get_num_branches(net)*t;
  
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

  // Construct constraints Ikm <= Imax
  results[0] = LINE_FLOW_construct(V_min[0], 
				   V_max[0], 
				   V_min[1], 
				   V_max[1],
				   BRANCH_get_g(br), 
				   BRANCH_get_b(br),
				   BRANCH_get_b_k(br)+BRANCH_get_b_m(br), // total line charging
				   1./BRANCH_get_ratio(br,t),
				   BRANCH_get_phase(br,t)*180./PI, // degrees
				   BRANCH_get_ratingA(br),
				   1,     // "k" side
				   NULL);

  // Construct constraints Imk <= Imax
  results[1] = LINE_FLOW_construct(V_min[0], 
				   V_max[0], 
				   V_min[1], 
				   V_max[1],
				   BRANCH_get_g(br), 
				   BRANCH_get_b(br),
				   BRANCH_get_b_k(br)+BRANCH_get_b_m(br), // total line charging
				   1./BRANCH_get_ratio(br,t),
				   BRANCH_get_phase(br,t)*180./PI, // degrees
				   BRANCH_get_ratingA(br),
				   2,     // "m" side
				   NULL);

  // Save results
  for (i = 0; i < 2; i++) {
    if (data->results[2*BRANCH_get_index(br)+i+offset])
      LINE_FLOW_free_results(data->results[2*BRANCH_get_index(br)+i+offset]); // delete existing
    data->results[2*BRANCH_get_index(br)+i+offset] = results[i];              // save new
  }
 
  // Check success
  for (i = 0; i < 2; i++) {
    if (LINE_FLOW_get_flag(results[i]) != success && LINE_FLOW_get_flag(results[i]) != non_binding) {
      CONSTR_set_error(c,LINE_FLOW_get_message(results[i]));
      return;
    }
  }
 
  // Count rows and nnz
  for (i = 0; i < 2; i++) {
    if (LINE_FLOW_get_flag(results[i]) == success) {
      num_constr = LINE_FLOW_get_number_constraints(results[i]);
      for (k = 0; k < 2; k++) {
	if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) // v_mag variable
	  (*G_nnz) += num_constr;
	if (BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VANG)) // v_ang variable
	  (*G_nnz) += num_constr;
      }
      (*G_row) += num_constr;
    }
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
  int i;
  int j;
  int k;
  int* G_nnz;
  int* G_row;
  Net* net;
  Bus* bus[2];
  int offset;
  int num_constr;
  LINE_FLOW_Results* results[2];
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
  offset = 2*NET_get_num_branches(net)*t;
  
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
  results[0] = data->results[2*BRANCH_get_index(br)+0+offset];
  results[1] = data->results[2*BRANCH_get_index(br)+1+offset];
  
  // Construct l <= Gx <= u
  for (i = 0; i < 2; i++) {

    if (LINE_FLOW_get_flag(results[i]) == success) {

      num_constr = LINE_FLOW_get_number_constraints(results[i]);
      A = LINE_FLOW_get_A_matrix(results[i]);
      b = LINE_FLOW_get_b_vector(results[i]);
    
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
	LINE_FLOW_free_results(data->results[i]);
    }
    free(data->results);
    free(data);
  }
   
  // Set data
  CONSTR_set_data(c,NULL);
}
