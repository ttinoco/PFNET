/** @file constr_REG_GEN.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_GEN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_GEN.h>

Constr* CONSTR_REG_GEN_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_REG_GEN_init);
  CONSTR_set_func_count_step(c,&CONSTR_REG_GEN_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_REG_GEN_allocate);
  CONSTR_set_func_clear(c,&CONSTR_REG_GEN_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_REG_GEN_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_REG_GEN_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_REG_GEN_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_REG_GEN_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_REG_GEN_init(Constr* c) {

  // Local variables
  Net* net;
  int num_Jconstr;

  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 2*(NET_get_num_buses_reg_by_gen(net)-NET_get_num_slack_buses(net))*NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_name(c,"voltage regulation by generators");
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_GEN_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_A_nnz(c,0);
  CONSTR_set_J_nnz(c,0);
  CONSTR_set_A_row(c,0);
  CONSTR_set_J_row(c,0);
  CONSTR_clear_H_nnz(c);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_REG_GEN_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !J_row || !H_nnz || !bus_counted)
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

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) { // regulator and not slack
	
	// Linear
	//*******
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	  // A
	  (*A_nnz)++; // v
	}
	  
	// A
	(*A_nnz)++; // y
	(*A_nnz)++; // z

	// Count
	(*A_row)++;

	// Nonlinear
	//**********

	// J
	(*J_nnz)++; // dCompY/dy
	(*J_nnz)++; // dCompZ/dz
	
	// H
	H_nnz[*J_row]++;     // y and y (CompY)
	H_nnz[*J_row+1]++;   // z and z (CompZ)
	
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var
	    
	    // J
	    (*J_nnz)++; // dCompY/dQ
	    (*J_nnz)++; // dCompZ/dQ
	    
	    // H
	    H_nnz[*J_row]++;   // Q and Q (CompY)
	    H_nnz[*J_row]++;   // y and Q (CompY)
	    
	    H_nnz[*J_row+1]++; // Q and Q (CompZ)
	    H_nnz[*J_row+1]++; // z and Q (CompZ)

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) { // Q1 var
		H_nnz[*J_row]++;   // Q and Q1 (CompY)
		H_nnz[*J_row+1]++; // Q and Q1 (CompZ)
	      }
	    }	    
	  }
	}

	// Count
	(*J_row)++; // dCompY
	(*J_row)++; // dCompZ

	// Num extra vars
	CONSTR_set_num_extra_vars(c,*J_row);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_allocate(Constr* c) {

  // Local variables
  int A_nnz;
  int J_nnz;
  int A_row;
  int J_row;
  int* H_nnz;
  Mat* H;
  int num_vars;
  int num_extra_vars;
  int i;

  A_nnz = CONSTR_get_A_nnz(c);
  J_nnz = CONSTR_get_J_nnz(c);
  A_row = CONSTR_get_A_row(c);
  J_row = CONSTR_get_J_row(c);
  H_nnz = CONSTR_get_H_nnz(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_extra_vars = CONSTR_get_num_extra_vars(c);

  // Extra vars
  CONSTR_set_l_extra_vars(c,VEC_new(num_extra_vars));
  CONSTR_set_u_extra_vars(c,VEC_new(num_extra_vars));
  CONSTR_set_init_extra_vars(c,VEC_new(num_extra_vars));

  // G u l
  CONSTR_set_G(c,MAT_new(0,num_vars+num_extra_vars,0));
  CONSTR_set_u(c,VEC_new(0));
  CONSTR_set_l(c,VEC_new(0));
  
  // b
  CONSTR_set_b(c,VEC_new(A_row));

  // A
  CONSTR_set_A(c,MAT_new(A_row,                   // size1 (rows)
			 num_vars+num_extra_vars, // size2 (cols)
			 A_nnz));                 // nnz

  // f
  CONSTR_set_f(c,VEC_new(J_row));

  // J
  CONSTR_set_J(c,MAT_new(J_row,                   // size1 (rows)
			 num_vars+num_extra_vars, // size2 (cols)
			 J_nnz));                 // nnz

  // H
  CONSTR_allocate_H_array(c,J_row);
  for (i = 0; i < J_row; i++) {
    H = CONSTR_get_H_single(c,i);
    MAT_set_nnz(H,H_nnz[i]);
    MAT_set_size1(H,num_vars+num_extra_vars);
    MAT_set_size2(H,num_vars+num_extra_vars);
    MAT_set_row_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)calloc(H_nnz[i],sizeof(REAL)));
  }
}

void CONSTR_REG_GEN_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hy;
  Mat* Hz;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int index_y;
  int index_z;
  int k;
  int T;
  int num_vars;

  // Number of periods and vars
  T = BRANCH_get_num_periods(br);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // Constr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !H_array || !J_row || !H_nnz || !bus_counted)
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

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) { // regulator and not slack

	// Hessians
	Hy = MAT_array_get(H_array,*J_row);
	Hz = MAT_array_get(H_array,*J_row+1);

	// Indices
	index_y = num_vars+(*J_row);
	index_z = num_vars+(*J_row+1);

	// Linear
	//*******

	// b
	VEC_set(b,*A_row,BUS_get_v_set(bus,t));
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	  
	  // A
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++; // v
	}
	else
	  VEC_add_to_entry(b,*A_row,-BUS_get_v_mag(bus,t));

	// A
	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,index_y);
	MAT_set_d(A,*A_nnz,-1.);
	(*A_nnz)++; // y

	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,index_z);
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++; // z
	
	// Count
	(*A_row)++;

	// Nonlinear
	//**********

	// J
	MAT_set_i(J,*J_nnz,*J_row);
	MAT_set_j(J,*J_nnz,index_y);
	(*J_nnz)++; // dCompY/dy

	MAT_set_i(J,*J_nnz,*J_row+1);
	MAT_set_j(J,*J_nnz,index_z);
	(*J_nnz)++; // dCompZ/dz

	// H
	MAT_set_i(Hy,H_nnz[*J_row],index_y);
	MAT_set_j(Hy,H_nnz[*J_row],index_y);
	H_nnz[*J_row]++; // y and y (CompY)

	MAT_set_i(Hz,H_nnz[*J_row+1],index_z);
	MAT_set_j(Hz,H_nnz[*J_row+1],index_z);
	H_nnz[*J_row+1]++; // z and z (CompZ)

	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var

	    // J
	    MAT_set_i(J,*J_nnz,*J_row);
	    MAT_set_j(J,*J_nnz,GEN_get_index_Q(rg,t));
	    (*J_nnz)++; // dCompY/dQ
	    
	    MAT_set_i(J,*J_nnz,*J_row+1);
	    MAT_set_j(J,*J_nnz,GEN_get_index_Q(rg,t));
	    (*J_nnz)++; // dCompZ/dQ
	    
	    // H
	    MAT_set_i(Hy,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	    MAT_set_j(Hy,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	    H_nnz[*J_row]++; // Q and Q (CompY)
	    
	    MAT_set_i(Hy,H_nnz[*J_row],index_y);
	    MAT_set_j(Hy,H_nnz[*J_row],GEN_get_index_Q(rg,t));
	    H_nnz[*J_row]++; // y and Q (CompY)

	    MAT_set_i(Hz,H_nnz[*J_row+1],GEN_get_index_Q(rg,t));
	    MAT_set_j(Hz,H_nnz[*J_row+1],GEN_get_index_Q(rg,t));
	    H_nnz[*J_row+1]++; // Q and Q (CompZ)

	    MAT_set_i(Hz,H_nnz[*J_row+1],index_z);
	    MAT_set_j(Hz,H_nnz[*J_row+1],GEN_get_index_Q(rg,t));
	    H_nnz[*J_row+1]++; // z and Q (CompZ)

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {

		MAT_set_i(Hy,H_nnz[*J_row],GEN_get_index_Q(rg,t));
		MAT_set_j(Hy,H_nnz[*J_row],GEN_get_index_Q(rg1,t));
		H_nnz[*J_row]++; // Q and Q1 (CompY)
		
		MAT_set_i(Hz,H_nnz[*J_row+1],GEN_get_index_Q(rg,t));
		MAT_set_j(Hz,H_nnz[*J_row+1],GEN_get_index_Q(rg1,t));
		H_nnz[*J_row+1]++; // Q and Q1 (CompZ)
	      }
	    }
	  }
	}
	 
	// Extra var limits
	VEC_set(CONSTR_get_l_extra_vars(c),*J_row,-CONSTR_REG_GEN_MAX_YZ);   // y
	VEC_set(CONSTR_get_l_extra_vars(c),*J_row+1,-CONSTR_REG_GEN_MAX_YZ); // z

	VEC_set(CONSTR_get_u_extra_vars(c),*J_row,CONSTR_REG_GEN_MAX_YZ);   // y
	VEC_set(CONSTR_get_u_extra_vars(c),*J_row+1,CONSTR_REG_GEN_MAX_YZ); // z
 
	// Count
	(*J_row)++; // dCompY
	(*J_row)++; // dCompZ
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* Hy;
  REAL* Hz;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL y;
  REAL z;
  REAL Qsum;
  REAL Qmin;
  REAL Qmax;
  REAL Qy;
  REAL Qz;
  REAL sqrt_termY;
  REAL sqrt_termZ;
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

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) { // regulator and not slack

	// Hessians
	Hy = MAT_get_data_array(MAT_array_get(H_array,*J_row));
	Hz = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));

	// Extra vars
	if (VEC_get_size(values_extra) > 0) {
	  y = VEC_get(values_extra,*J_row);
	  z = VEC_get(values_extra,*J_row+1);
	}
	else {
	  y = 0.;
	  z = 0.;
	}
	
	// Q value
	Qsum = 0;
	Qmax = 0;
	Qmin = 0;
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q))
	    Qsum += VEC_get(values,GEN_get_index_Q(rg,t)); // p.u.
	  else
	    Qsum += GEN_get_Q(rg,t); // p.u.
	  Qmax += GEN_get_Q_max(rg); // p.u.
	  Qmin += GEN_get_Q_min(rg); // p.u.
	}
	Qy = (Qsum-Qmin);
	Qz = (Qmax-Qsum);

	// Terms
	sqrt_termY = sqrt( Qy*Qy + y*y + 2*CONSTR_REG_GEN_PARAM );
	sqrt_termZ = sqrt( Qz*Qz + z*z + 2*CONSTR_REG_GEN_PARAM );

	// f
	f[*J_row] = Qy + y - sqrt_termY;   // CompY
	f[*J_row+1] = Qz + z - sqrt_termZ; // CompZ

	// J
	J[*J_nnz] = 1. - y/sqrt_termY;
	(*J_nnz)++; // dCompY/dy

	J[*J_nnz] = 1. - z/sqrt_termZ;
	(*J_nnz)++; // dCompZ/dz

	// H
	Hy[H_nnz[*J_row]] = -(Qy*Qy+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	H_nnz[*J_row]++; // y and y (CompY)

	Hz[H_nnz[*J_row+1]] = -(Qz*Qz+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	H_nnz[*J_row+1]++; // z and z (CompZ)
	
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Q var

	    // J
	    J[*J_nnz] = 1. - Qy/sqrt_termY;
	    (*J_nnz)++; // dcompY/dQ

	    J[*J_nnz] = -1. + Qz/sqrt_termZ;
	    (*J_nnz)++; // dcompZ/dQ

	    // H
	    Hy[H_nnz[*J_row]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	    H_nnz[*J_row]++; // Q and Q (CompY)
	    
	    Hy[H_nnz[*J_row]] = Qy*y/pow(sqrt_termY,3.);
	    H_nnz[*J_row]++; // y and Q (CompZ)

	    Hz[H_nnz[*J_row+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	    H_nnz[*J_row+1]++; // Q and Q (CompZ)

	    Hz[H_nnz[*J_row+1]] = -Qz*z/pow(sqrt_termZ,3.);
	    H_nnz[*J_row+1]++; // z and Q (CompZ)

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) { // Q1 var

		Hy[H_nnz[*J_row]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
		H_nnz[*J_row]++; // Q and Q1 (CompY)

		Hz[H_nnz[*J_row+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
		H_nnz[*J_row+1]++; // Q and Q1 (CompZ)
	      }
	    }
	  }
	}

	// Count
	(*J_row)++; // dCompY
	(*J_row)++; // dCompZ
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int* J_row;
  char* bus_counted;
  int bus_index_t[2];
  REAL lamCompY;
  REAL lamCompZ;
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J_row = CONSTR_get_J_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_row || !bus_counted)
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

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      if (BUS_is_regulated_by_gen(bus) && !BUS_is_slack(bus)) { // regulator and not slack

	lamCompY = VEC_get(sf,*J_row);
	(*J_row)++; // dCompY
	lamCompZ = VEC_get(sf,*J_row);
	(*J_row)++; // dCompZ

	if (fabs(lamCompY) > fabs(lamCompZ))
	  BUS_set_sens_v_reg_by_gen(bus,lamCompY,t);
	else
	  BUS_set_sens_v_reg_by_gen(bus,lamCompZ,t);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_free(Constr* c) {
  // Nothing
}
