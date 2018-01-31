/** @file constr_REG_GEN_EX.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_GEN_EX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_GEN_EX.h>

Constr* CONSTR_REG_GEN_EX_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_REG_GEN_EX_init);
  CONSTR_set_func_count_step(c,&CONSTR_REG_GEN_EX_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_REG_GEN_EX_allocate);
  CONSTR_set_func_clear(c,&CONSTR_REG_GEN_EX_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_REG_GEN_EX_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_REG_GEN_EX_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_REG_GEN_EX_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_REG_GEN_EX_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_REG_GEN_EX_init(Constr* c) {

  // Local variables
  Net* net;
  int num_Jconstr;

  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = (NET_get_num_reg_gens(net)-NET_get_num_slack_gens(net))*NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_name(c,"extended voltage regulation by generators");
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_GEN_EX_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_J_nnz(c,0);
  CONSTR_set_J_row(c,0);
  CONSTR_clear_H_nnz(c);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_REG_GEN_EX_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !J_row || !H_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    // Bus not counted yet
    if (!bus_counted[bus_index_t[k]]) {

      // Bus is regulated and not slack
      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) {
	
	// Regulating generators
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	    
	    // J
	    (*J_nnz)++; // dComp/dv
	    
	    // H
	    H_nnz[*J_row]++; // v and v
	  }

	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var
		    
	    // J
	    (*J_nnz)++; // dComp/dQ
	    
	    // H
	    H_nnz[*J_row]++; // Q and Q

	    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	      // H
	      H_nnz[*J_row]++; // Q and v
	    }
	  }
	  
	  // Count row
	  (*J_row)++; // dComp
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_EX_allocate(Constr* c) {

  // Local variables
  int J_nnz;
  int J_row;
  int* H_nnz;
  Mat* H;
  int num_vars;
  int i;

  J_nnz = CONSTR_get_J_nnz(c);
  J_row = CONSTR_get_J_row(c);
  H_nnz = CONSTR_get_H_nnz(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // G u l
  CONSTR_set_G(c,MAT_new(0,num_vars,0));
  CONSTR_set_u(c,VEC_new(0));
  CONSTR_set_l(c,VEC_new(0));
  
  // b
  CONSTR_set_b(c,VEC_new(0));

  // A
  CONSTR_set_A(c,MAT_new(0,num_vars,0));

  // f
  CONSTR_set_f(c,VEC_new(J_row));

  // J
  CONSTR_set_J(c,MAT_new(J_row,    // size1 (rows)
			 num_vars, // size2 (cols)
			 J_nnz));  // nnz

  // H
  CONSTR_allocate_H_array(c,J_row);
  for (i = 0; i < J_row; i++) {
    H = CONSTR_get_H_single(c,i);
    MAT_set_nnz(H,H_nnz[i]);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)calloc(H_nnz[i],sizeof(REAL)));
  }
}

void CONSTR_REG_GEN_EX_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Mat* J;
  Mat* H_array;
  Mat* H;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods and vars
  T = BRANCH_get_num_periods(br);

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !H_array || !J_row || !H_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    // Bus not counted yet
    if (!bus_counted[bus_index_t[k]]) {

      // Bus regulated and not slack
      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) {

	// Regulating generators
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  
	  // Hessians
	  H = MAT_array_get(H_array,*J_row);
	
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	  	  
	    // J
	    MAT_set_i(J,*J_nnz,*J_row);
	    MAT_set_j(J,*J_nnz,BUS_get_index_v_mag(bus,t));
	    (*J_nnz)++; // dComp/dv

	    // H
	    MAT_set_i(H,H_nnz[*J_row],BUS_get_index_v_mag(bus,t));
	    MAT_set_j(H,H_nnz[*J_row],BUS_get_index_v_mag(bus,t));
	    H_nnz[*J_row]++; // v and v
	  }
	    
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var
	    
	    // J
	    MAT_set_i(J,*J_nnz,*J_row);
	    MAT_set_j(J,*J_nnz,GEN_get_index_Q(rg,t));
	    (*J_nnz)++; // dComp/dQ
	    
	    // H
	    MAT_set_i(H,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	    MAT_set_j(H,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	    H_nnz[*J_row]++; // Q and Q
	    
	    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	      
	      // H
	      MAT_set_i(H,H_nnz[*J_row],BUS_get_index_v_mag(bus,t));
	      MAT_set_j(H,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	      H_nnz[*J_row]++; // Q and v
	    }
	  }
	  
	  // Count row
	  (*J_row)++; // dComp
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_EX_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* H;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL eta = CONSTR_REG_GEN_EX_PARAM;
  REAL v;
  REAL vs;
  REAL a;
  REAL b;
  REAL Q;
  REAL Qmin;
  REAL Qmax;
  REAL sqrt_term;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!f || !J || !J_nnz || !J_row || !H_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    // Bus not counted yet
    if (!bus_counted[bus_index_t[k]]) {

      // Bus regulated and not slack
      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) {

	// Regulating generators
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {

	  // Hessian
	  H = MAT_get_data_array(MAT_array_get(H_array,*J_row));

	  // v
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	    v = VEC_get(values,BUS_get_index_v_mag(bus,t));
	  else
	    v = BUS_get_v_mag(bus,t);
	  vs = BUS_get_v_set(bus,t);
	
	  // Q values
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q))
	    Q = VEC_get(values,GEN_get_index_Q(rg,t)); // p.u.
	  else
	    Q = GEN_get_Q(rg,t);    // p.u.
	  Qmax = GEN_get_Q_max(rg); // p.u.
	  Qmin = GEN_get_Q_min(rg); // p.u.

	  // a b
	  a = (Qmax-Q)*(Q-Qmin);
	  b = (v-vs)*(v-vs);

	  // Sqrt term
	  sqrt_term = sqrt( a*a + b*b + 2*eta );

	  // f
	  f[*J_row] = a + b - sqrt_term;   // Comp

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	    // J
	    J[*J_nnz] = (1. - b/sqrt_term)*2*(v-vs);
	    (*J_nnz)++; // dComp/dv

	    // H
	    H[H_nnz[*J_row]] = -((a*a+2*eta)/pow(sqrt_term,3.))*pow(2*(v-vs),2) + 2*(1.-b/sqrt_term);
	    H_nnz[*J_row]++; // v and v
	  }
	
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var

	    // J
	    J[*J_nnz] = (1. - a/sqrt_term)*(Qmax+Qmin-2.*Q);
	    (*J_nnz)++; // dcomp/dQ

	    // H
	    H[H_nnz[*J_row]] = -((b*b+2*eta)/pow(sqrt_term,3.))*pow(Qmax+Qmin-2*Q,2) - 2*(1.-a/sqrt_term);
	    H_nnz[*J_row]++; // Q and Q

	    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	      // H
	      H[H_nnz[*J_row]] = (a*b/pow(sqrt_term,3.))*(Qmax+Qmin-2*Q)*2*(v-vs);
	      H_nnz[*J_row]++; // Q and v
	    }
	  }

	  // Count row
	  (*J_row)++; // dComp
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_EX_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_REG_GEN_EX_free(Constr* c) {
  // Nothing
}
