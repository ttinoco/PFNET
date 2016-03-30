/** @file constr_BOUND.c
 *  @brief This file defines the data structure and routines associated with the constraint of type BOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_BOUND.h>

void CONSTR_BOUND_init(Constr* c) {

  // Init
  CONSTR_set_Hcounter(c,NULL,0);
  CONSTR_set_data(c,NULL);
}

void CONSTR_BOUND_clear(Constr* c) {
  
  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));
  
  // Counters
  CONSTR_set_Jcounter(c,0);
  CONSTR_set_branch_counter(c,0);
  
  // Flags
  CONSTR_clear_bus_counted(c);  
}

void CONSTR_BOUND_count_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int* Jcounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  
  // Constr data
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Jcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
 
  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*Jcounter)++; // upper bound
    (*Jcounter)++; // lower bound
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*Jcounter)++; // upper bound
    (*Jcounter)++; // lower bound
  }
  
  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	(*Jcounter)++; // upper bound
	(*Jcounter)++; // lower bound
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*Jcounter)++; // upper bound
	(*Jcounter)++; // lower bound
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
      }
      
      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Shunt suscepatnace (b)
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
      }      
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_BOUND_allocate(Constr* c) {

  // Local variables
  int Jcounter;
  Mat* H_array;
  Mat* H;
  int num_vars;
  int i;
  
  Jcounter = CONSTR_get_Jcounter(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));
  
  // f
  CONSTR_set_f(c,VEC_new(Jcounter));

  // J
  CONSTR_set_J(c,MAT_new(Jcounter,   // size1 (rows)
			 num_vars,   // size2 (cols)
			 Jcounter)); // nnz
  
  // H
  H_array = MAT_array_new(Jcounter);
  CONSTR_set_H_array(c,H_array,Jcounter);
  for (i = 0; i < Jcounter; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,1);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(1,sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(1,sizeof(int)));
    MAT_set_data_array(H,(REAL*)malloc(1*sizeof(REAL)));
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,   // size1 (rows)
				  num_vars,   // size2 (cols)
				  Jcounter)); // nnz
}

void CONSTR_BOUND_analyze_branch(Constr* c, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  Mat* J;
  Mat* H_array;
  Mat* H;
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int* Jcounter;
  char* bus_counted;
  int bus_index[2];
  int index_var;
  int k;

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  CONSTR_inc_branch_counter(c);

  // Check pointers
  if (!Jcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Branch
  //*******

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {

    index_var = BRANCH_get_index_ratio(br);
    
    // J
    MAT_set_i(J,*Jcounter,*Jcounter);
    MAT_set_j(J,*Jcounter,index_var);

    MAT_set_i(J,*Jcounter+1,*Jcounter+1);
    MAT_set_j(J,*Jcounter+1,index_var);

    // H
    H = MAT_array_get(H_array,*Jcounter);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    H = MAT_array_get(H_array,*Jcounter+1);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    (*Jcounter)++;
    (*Jcounter)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {

    index_var = BRANCH_get_index_phase(br);

    // J
    MAT_set_i(J,*Jcounter,*Jcounter);
    MAT_set_j(J,*Jcounter,index_var);

    MAT_set_i(J,*Jcounter+1,*Jcounter+1);
    MAT_set_j(J,*Jcounter+1,index_var);

    // H
    H = MAT_array_get(H_array,*Jcounter);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    H = MAT_array_get(H_array,*Jcounter+1);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);
    
    (*Jcounter)++;
    (*Jcounter)++;
  }
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

	index_var = BUS_get_index_v_mag(bus);

	// J
	MAT_set_i(J,*Jcounter,*Jcounter);
	MAT_set_j(J,*Jcounter,index_var);

	MAT_set_i(J,*Jcounter+1,*Jcounter+1);
	MAT_set_j(J,*Jcounter+1,index_var);
	
	// H
	H = MAT_array_get(H_array,*Jcounter);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	H = MAT_array_get(H_array,*Jcounter+1);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);
	
	(*Jcounter)++;
	(*Jcounter)++;
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {

	index_var = BUS_get_index_v_ang(bus);

	// J
	MAT_set_i(J,*Jcounter,*Jcounter);
	MAT_set_j(J,*Jcounter,index_var);

	MAT_set_i(J,*Jcounter+1,*Jcounter+1);
	MAT_set_j(J,*Jcounter+1,index_var);
	
	// H
	H = MAT_array_get(H_array,*Jcounter);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	H = MAT_array_get(H_array,*Jcounter+1);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);
	
	(*Jcounter)++;
	(*Jcounter)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {

	  index_var = GEN_get_index_P(gen);
	  
	  // J
	  MAT_set_i(J,*Jcounter,*Jcounter);
	  MAT_set_j(J,*Jcounter,index_var);

	  MAT_set_i(J,*Jcounter+1,*Jcounter+1);
	  MAT_set_j(J,*Jcounter+1,index_var);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {

	  index_var = GEN_get_index_Q(gen);
	  
	  // J
	  MAT_set_i(J,*Jcounter,*Jcounter);
	  MAT_set_j(J,*Jcounter,index_var);

	  MAT_set_i(J,*Jcounter+1,*Jcounter+1);
	  MAT_set_j(J,*Jcounter+1,index_var);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
		
	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  
	  index_var = SHUNT_get_index_b(shunt);
	  
	  // J
	  MAT_set_i(J,*Jcounter,*Jcounter);
	  MAT_set_j(J,*Jcounter,index_var);

	  MAT_set_i(J,*Jcounter+1,*Jcounter+1);
	  MAT_set_j(J,*Jcounter+1,index_var);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
      }
    }	  
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }

  // Done
  if (CONSTR_get_branch_counter(c) == NET_get_num_branches(CONSTR_get_network(c))) {
    
    // Ensure lower triangular and save struct of H comb
    Hi_comb = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hj_comb = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (k = 0; k < CONSTR_get_H_array_size(c); k++) {
      Hi = MAT_get_row_array(MAT_array_get(H_array,k));
      Hj = MAT_get_col_array(MAT_array_get(H_array,k));
      Hi_comb[k] = Hi[0];
      Hj_comb[k] = Hj[0];
    }
  }
}

void CONSTR_BOUND_eval_branch(Constr* c, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  Mat* H_array;
  REAL* f;
  REAL* J;
  Mat* H;
  int* Jcounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  REAL u;
  REAL umin;
  REAL umax;
  REAL du;
  REAL a1;
  REAL a2;
  REAL b;
  REAL eps;
  REAL sqrterm1;
  REAL sqrterm2;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  CONSTR_inc_branch_counter(c);

  // Check pointers
  if (!f || !J || !Jcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Param
  eps = CONSTR_BOUND_PARAM;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Branch
  //*******

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    
    u = VEC_get(var_values,BRANCH_get_index_ratio(br));
    umax = BRANCH_get_ratio_max(br);
    umin = BRANCH_get_ratio_min(br);
    du = (umax-umin > eps) ? umax-umin : eps;
    
    a1 = umax-u;
    a2 = u-umin;
    b = eps*eps/du;
    sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
    sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

    // f
    f[*Jcounter]   = a1 + b - sqrterm1; // upper
    f[*Jcounter+1] = a2 + b - sqrterm2; // lower
    
    // J
    J[*Jcounter]   = -(1-a1/sqrterm1);
    J[*Jcounter+1] = (1-a2/sqrterm2);

    // H
    H = MAT_array_get(H_array,*Jcounter);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

    H = MAT_array_get(H_array,*Jcounter+1);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

    (*Jcounter)++;
    (*Jcounter)++;
  }
  
  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    
    u = VEC_get(var_values,BRANCH_get_index_phase(br));
    umax = BRANCH_get_phase_max(br);
    umin = BRANCH_get_phase_min(br);
    du = (umax-umin > eps) ? umax-umin : eps;
    
    a1 = umax-u;
    a2 = u-umin;
    b = eps*eps/du;
    sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
    sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

    // f
    f[*Jcounter]   = a1 + b - sqrterm1; // upper
    f[*Jcounter+1] = a2 + b - sqrterm2; // lower
    
    // J
    J[*Jcounter]   = -(1-a1/sqrterm1);
    J[*Jcounter+1] = (1-a2/sqrterm2);

    // H
    H = MAT_array_get(H_array,*Jcounter);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

    H = MAT_array_get(H_array,*Jcounter+1);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

    (*Jcounter)++;
    (*Jcounter)++;    
  }
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
	u = VEC_get(var_values,BUS_get_index_v_mag(bus));
	umax = BUS_get_v_max(bus);
	umin = BUS_get_v_min(bus);
	du = (umax-umin > eps) ? umax-umin : eps;
	
	a1 = umax-u;
	a2 = u-umin;
	b = eps*eps/du;
	sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	sqrterm2 = sqrt(a2*a2+b*b+eps*eps);
	
	// f
	f[*Jcounter]   = a1 + b - sqrterm1; // upper
	f[*Jcounter+1] = a2 + b - sqrterm2; // lower
	
	// J
	J[*Jcounter]   = -(1-a1/sqrterm1);
	J[*Jcounter+1] = (1-a2/sqrterm2);
	
	// H
	H = MAT_array_get(H_array,*Jcounter);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));
	
	H = MAT_array_get(H_array,*Jcounter+1);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));
	
	(*Jcounter)++;
	(*Jcounter)++;	
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	
	u = VEC_get(var_values,BUS_get_index_v_ang(bus));
	umax = PI;
	umin = -PI;
	du = (umax-umin > eps) ? umax-umin : eps;
	
	a1 = umax-u;
	a2 = u-umin;
	b = eps*eps/du;
	sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	sqrterm2 = sqrt(a2*a2+b*b+eps*eps);
	
	// f
	f[*Jcounter]   = a1 + b - sqrterm1; // upper
	f[*Jcounter+1] = a2 + b - sqrterm2; // lower
	
	// J
	J[*Jcounter]   = -(1-a1/sqrterm1);
	J[*Jcounter+1] = (1-a2/sqrterm2);
	
	// H
	H = MAT_array_get(H_array,*Jcounter);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));
	
	H = MAT_array_get(H_array,*Jcounter+1);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));
	
	(*Jcounter)++;
	(*Jcounter)++;	
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  
	  u = VEC_get(var_values,GEN_get_index_P(gen));
	  umax = GEN_get_P_max(gen);
	  umin = GEN_get_P_min(gen);
	  du = (umax-umin > eps) ? umax-umin : eps;
	  
	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);
	  
	  // f
	  f[*Jcounter]   = a1 + b - sqrterm1; // upper
	  f[*Jcounter+1] = a2 + b - sqrterm2; // lower
	  
	  // J
	  J[*Jcounter]   = -(1-a1/sqrterm1);
	  J[*Jcounter+1] = (1-a2/sqrterm2);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));
	  
	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  
	  u = VEC_get(var_values,GEN_get_index_Q(gen));
	  umax = GEN_get_Q_max(gen);
	  umin = GEN_get_Q_min(gen);
	  du = (umax-umin > eps) ? umax-umin : eps;
	  
	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);
	  
	  // f
	  f[*Jcounter]   = a1 + b - sqrterm1; // upper
	  f[*Jcounter+1] = a2 + b - sqrterm2; // lower
	  
	  // J
	  J[*Jcounter]   = -(1-a1/sqrterm1);
	  J[*Jcounter+1] = (1-a2/sqrterm2);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));
	  
	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  
	  u = VEC_get(var_values,SHUNT_get_index_b(shunt));
	  umax = SHUNT_get_b_max(shunt);
	  umin = SHUNT_get_b_min(shunt);
	  du = (umax-umin > eps) ? umax-umin : eps;
	  
	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);
	  
	  // f
	  f[*Jcounter]   = a1 + b - sqrterm1; // upper
	  f[*Jcounter+1] = a2 + b - sqrterm2; // lower
	  
	  // J
	  J[*Jcounter]   = -(1-a1/sqrterm1);
	  J[*Jcounter+1] = (1-a2/sqrterm2);
	  
	  // H
	  H = MAT_array_get(H_array,*Jcounter);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));
	  
	  H = MAT_array_get(H_array,*Jcounter+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));
	  
	  (*Jcounter)++;
	  (*Jcounter)++;
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_BOUND_store_sens_branch(Constr* c, Branch* br, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int* Jcounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  
  // Constr data
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Jcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
 
  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*Jcounter)++; // upper bound
    (*Jcounter)++; // lower bound
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*Jcounter)++; // upper bound
    (*Jcounter)++; // lower bound
  }
  
  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	BUS_set_sens_v_mag_u_bound(bus,VEC_get(sf,*Jcounter));
	(*Jcounter)++; // upper bound
	BUS_set_sens_v_mag_l_bound(bus,VEC_get(sf,*Jcounter));
	(*Jcounter)++; // lower bound
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*Jcounter)++; // upper bound
	(*Jcounter)++; // lower bound
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
      }
      
      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*Jcounter)++; // upper bound
	  (*Jcounter)++; // lower bound
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }  
}

void CONSTR_BOUND_free(Constr* c) {
  // Nothing
}
