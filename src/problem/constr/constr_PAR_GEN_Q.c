/** @file constr_PAR_GEN_Q.c
 *  @brief This file defines the data structure and routines associated with the constraint of type PAR_GEN_Q.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_PAR_GEN_Q.h>

struct Constr_PAR_GEN_Q_Data {

  int type;
};

Constr* CONSTR_PAR_GEN_Q_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_PAR_GEN_Q_init);
  CONSTR_set_func_count_step(c, &CONSTR_PAR_GEN_Q_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_PAR_GEN_Q_allocate);
  CONSTR_set_func_clear(c, &CONSTR_PAR_GEN_Q_clear);
  CONSTR_set_func_analyze_step(c, &CONSTR_PAR_GEN_Q_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_PAR_GEN_Q_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_PAR_GEN_Q_store_sens_step);
  CONSTR_set_func_set_parameter(c, &CONSTR_PAR_GEN_Q_set_parameter);
  CONSTR_set_func_free(c, &CONSTR_PAR_GEN_Q_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_PAR_GEN_Q_init(Constr* c) {
  
  // Local variables
  Constr_PAR_GEN_Q_Data* data;
  
  // Data
  data = (Constr_PAR_GEN_Q_Data*)malloc(sizeof(Constr_PAR_GEN_Q_Data));
  data->type = CONSTR_PAR_GEN_Q_TYPE_RANGE;

  // Init
  CONSTR_set_name(c,"generator reactive power participation");
  CONSTR_set_data(c,data);
}

void CONSTR_PAR_GEN_Q_clear(Constr* c) {

  // Counters
  CONSTR_set_A_nnz(c,0);
  CONSTR_set_A_row(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_PAR_GEN_Q_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen1;
  Gen* gen2;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  Constr_PAR_GEN_Q_Data* data;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_PAR_GEN_Q_Data*)CONSTR_get_data(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus_counted || !data)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Reactive power of regulating generators
      if (BUS_is_regulated_by_gen(bus)) {


	// TYPE: RANGE or FRACTION
	//************************
	if (data->type == CONSTR_PAR_GEN_Q_TYPE_RANGE ||
	    data->type == CONSTR_PAR_GEN_Q_TYPE_FRACTION) {

	  // Reference
	  for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
	    if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q))
	      break;
	  }

	  // Constraint for each pair
	  for (gen2 = GEN_get_reg_next(gen1); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
	    if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_Q)) {
	      (*A_nnz)++;
	      (*A_nnz)++;
	      (*A_row)++;
	    }
	  }
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_Q_allocate(Constr* c) {

  // Local variables
  int num_constr;
  int num_vars;
  int A_nnz;

  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_constr = CONSTR_get_A_row(c);
  A_nnz = CONSTR_get_A_nnz(c);

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // G u l
  CONSTR_set_G(c,MAT_new(0,num_vars,0));
  CONSTR_set_u(c,VEC_new(0));
  CONSTR_set_l(c,VEC_new(0));

  // b
  CONSTR_set_b(c,VEC_new(num_constr));

  // A
  CONSTR_set_A(c,MAT_new(num_constr, // size1 (rows)
			 num_vars,   // size2 (rows)
			 A_nnz)); // nnz
}

void CONSTR_PAR_GEN_Q_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen1;
  Gen* gen2;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  Constr_PAR_GEN_Q_Data* data;
  Vec* b;
  Mat* A;
  int i;
  int T;
  
  REAL Qmin1;
  REAL Qmin2;
  REAL dQ1;
  REAL dQ2;

  REAL alpha1;
  REAL alpha2;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  data = (Constr_PAR_GEN_Q_Data*)CONSTR_get_data(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus_counted || !data)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Reactive power of regulating generators
      if (BUS_is_regulated_by_gen(bus)) {
	
	// TYPE: RANGE
	//************
	if (data->type == CONSTR_PAR_GEN_Q_TYPE_RANGE) {

	  // Get reference
	  for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
	    if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q))
	      break;
	  }

	  Qmin1 = GEN_get_Q_min(gen1);
	  dQ1 = GEN_get_Q_max(gen1)-Qmin1;
	  if (dQ1 < CONSTR_PAR_GEN_Q_PARAM)
	    dQ1 = CONSTR_PAR_GEN_Q_PARAM;
	  
	  // Add constraint for each pair
	  for (gen2 = GEN_get_reg_next(gen1); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
	    
	    if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_Q)) {
	    
	      Qmin2 = GEN_get_Q_min(gen2);
	      dQ2 = GEN_get_Q_max(gen2)-Qmin2;
	      if (dQ2 < CONSTR_PAR_GEN_Q_PARAM)
		dQ2 = CONSTR_PAR_GEN_Q_PARAM;
	      
	      VEC_set(b,*A_row,Qmin1/dQ1-Qmin2/dQ2);
	      
	      MAT_set_i(A,*A_nnz,*A_row);
	      MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen1,t));
	      MAT_set_d(A,*A_nnz,1./dQ1);
	      (*A_nnz)++;
	      
	      MAT_set_i(A,*A_nnz,*A_row);
	      MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen2,t));
	      MAT_set_d(A,*A_nnz,-1./dQ2);
	      (*A_nnz)++;
	      
	      (*A_row)++;
	    }
	  }
	}
	
	// TYPE: FRACTION
	//***************
	else if (data->type == CONSTR_PAR_GEN_Q_TYPE_FRACTION) {

	  // Get reference
	  for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
	    if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q))
	      break;
	  }

	  alpha1 = GEN_get_Q_par(gen1);
	  if (alpha1 < CONSTR_PAR_GEN_Q_PARAM)
	    alpha1 = CONSTR_PAR_GEN_Q_PARAM;
	  
	  // Add constraint for each pair
	  for (gen2 = GEN_get_reg_next(gen1); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
	    
	    if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_Q)) {
	    
	      alpha2 = GEN_get_Q_par(gen2);
	      if (alpha2 < CONSTR_PAR_GEN_Q_PARAM)
		alpha2 = CONSTR_PAR_GEN_Q_PARAM;
	      
	      VEC_set(b,*A_row,0.);
	      
	      MAT_set_i(A,*A_nnz,*A_row);
	      MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen1,t));
	      MAT_set_d(A,*A_nnz,alpha2);
	      (*A_nnz)++;
	      
	      MAT_set_i(A,*A_nnz,*A_row);
	      MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen2,t));
	      MAT_set_d(A,*A_nnz,-alpha1);
	      (*A_nnz)++;
	      
	      (*A_row)++;
	    }
	  }	  
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_PAR_GEN_Q_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_PAR_GEN_Q_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_PAR_GEN_Q_set_parameter(Constr* c, char* key, void* value) {

  // Local variables
  Constr_PAR_GEN_Q_Data* data = (Constr_PAR_GEN_Q_Data*)CONSTR_get_data(c);

  // Check
  if (!data)
    return;

  // Set
  if (strcmp(key,"type") == 0) {
    data->type = *((int*)value);
  }

  else // unknown
    CONSTR_set_error(c,"invalid parameter");
}

void CONSTR_PAR_GEN_Q_free(Constr* c) {

  // Local variables
  Constr_PAR_GEN_Q_Data* data;

  // Get data
  data = (Constr_PAR_GEN_Q_Data*)CONSTR_get_data(c);

  // Free
  if (data)
    free(data);

  // Set data
  CONSTR_set_data(c,NULL);  
}
