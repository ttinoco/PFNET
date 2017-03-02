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

void CONSTR_REG_GEN_init(Constr* c) {

  // Local variables
  Net* net;
  int num_Jconstr;

  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 2*(NET_get_num_buses_reg_by_gen(net)-NET_get_num_slack_buses(net))*NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
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
  CONSTR_set_Aconstr_index(c,0);
  CONSTR_set_Jconstr_index(c,0);
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
  int* Aconstr_index;
  int* Jconstr_index;
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
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !Aconstr_index ||
      !Jconstr_index || !H_nnz || !bus_counted)
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

      if (BUS_is_regulated_by_gen(bus) && // reg gen
	  !BUS_is_slack(bus)) {           // not slack

	// Linear
	//*******
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && // v var
	    BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	  // A
	  (*A_nnz)++; // v
	  (*A_nnz)++; // y
	  (*A_nnz)++; // z

	  // Inc A cosntr index
	  (*Aconstr_index)++;
	}

	// Nonlinear
	//**********
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	  // J
	  (*J_nnz)++; // dCompY/dy
	  (*J_nnz)++; // dCompZ/dz

	  // H
	  H_nnz[*Jconstr_index]++;     // y and y
	  H_nnz[*Jconstr_index+1]++;   // z and z
	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {
	      H_nnz[*Jconstr_index]++;   // y and Q
	      H_nnz[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}

	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	    // J
	    (*J_nnz)++; // dcompY/dQ
	    (*J_nnz)++; // dcompZ/dQ

	    // H
	    H_nnz[*Jconstr_index]++;   // Q and Q
	    H_nnz[*Jconstr_index+1]++; // Q and Q
	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {
		H_nnz[*Jconstr_index]++;   // Q and Q1
		H_nnz[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
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
  int Aconstr_index;
  int Jconstr_index;
  int* H_nnz;
  Mat* H_array;
  Mat* H;
  int H_comb_nnz;
  int num_vars;
  int i;

  A_nnz = CONSTR_get_A_nnz(c);
  J_nnz = CONSTR_get_J_nnz(c);
  Aconstr_index = CONSTR_get_Aconstr_index(c);
  Jconstr_index = CONSTR_get_Jconstr_index(c);
  H_nnz = CONSTR_get_H_nnz(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // b
  CONSTR_set_b(c,VEC_new(Aconstr_index));

  // A
  CONSTR_set_A(c,MAT_new(Aconstr_index, // size1 (rows)
			 num_vars,      // size2 (cols)
			 A_nnz));    // nnz

  // f
  CONSTR_set_f(c,VEC_new(Jconstr_index));

  // J
  CONSTR_set_J(c,MAT_new(Jconstr_index, // size1 (rows)
			 num_vars,      // size2 (cols)
			 J_nnz));    // nnz

  // H
  H_comb_nnz = 0;
  H_array = MAT_array_new(Jconstr_index);
  CONSTR_set_H_array(c,H_array,Jconstr_index);
  for (i = 0; i < Jconstr_index; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,H_nnz[i]);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)calloc(H_nnz[i],sizeof(REAL)));
    H_comb_nnz += H_nnz[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
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
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int* A_nnz;
  int* J_nnz;
  int* Aconstr_index;
  int* Jconstr_index;
  int* H_nnz;
  int H_nnz_comb;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int m;
  int temp;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !Aconstr_index ||
      !Jconstr_index || !H_nnz || !bus_counted)
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

      if (BUS_is_regulated_by_gen(bus) &&  // reg by gen
	  !BUS_is_slack(bus)) {            // not slack

	// Hessians
	Hy = MAT_array_get(H_array,*Jconstr_index);
	Hz = MAT_array_get(H_array,*Jconstr_index+1);

	// Linear
	//*******
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && // v var
	    BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	  // b
	  VEC_set(b,*Aconstr_index,BUS_get_v_set(bus,t));

	  // A
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++; // v

	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BUS_get_index_y(bus,t));
	  MAT_set_d(A,*A_nnz,-1.);
	  (*A_nnz)++; // y

	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BUS_get_index_z(bus,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++; // z

	  // Inc A constr index
	  (*Aconstr_index)++;
	}

	// Nonlinear
	//**********
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {

	  // J
	  MAT_set_i(J,*J_nnz,*Jconstr_index);
	  MAT_set_j(J,*J_nnz,BUS_get_index_y(bus,t));
	  (*J_nnz)++; // dCompY/dy

	  MAT_set_i(J,*J_nnz,*Jconstr_index+1);
	  MAT_set_j(J,*J_nnz,BUS_get_index_z(bus,t));
	  (*J_nnz)++; // dCompZ/dz

	  // H
	  MAT_set_i(Hy,H_nnz[*Jconstr_index],BUS_get_index_y(bus,t));
	  MAT_set_j(Hy,H_nnz[*Jconstr_index],BUS_get_index_y(bus,t));
	  H_nnz[*Jconstr_index]++;     // y and y

	  MAT_set_i(Hz,H_nnz[*Jconstr_index+1],BUS_get_index_z(bus,t));
	  MAT_set_j(Hz,H_nnz[*Jconstr_index+1],BUS_get_index_z(bus,t));
	  H_nnz[*Jconstr_index+1]++;   // z and z

	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {

	      MAT_set_i(Hy,H_nnz[*Jconstr_index],BUS_get_index_y(bus,t));
	      MAT_set_j(Hy,H_nnz[*Jconstr_index],GEN_get_index_Q(rg,t));
	      H_nnz[*Jconstr_index]++;   // y and Q

	      MAT_set_i(Hz,H_nnz[*Jconstr_index+1],BUS_get_index_z(bus,t));
	      MAT_set_j(Hz,H_nnz[*Jconstr_index+1],GEN_get_index_Q(rg,t));
	      H_nnz[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}

	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	    // J
	    MAT_set_i(J,*J_nnz,*Jconstr_index);
	    MAT_set_j(J,*J_nnz,GEN_get_index_Q(rg,t));
	    (*J_nnz)++; // dcompY/dQ

	    MAT_set_i(J,*J_nnz,*Jconstr_index+1);
	    MAT_set_j(J,*J_nnz,GEN_get_index_Q(rg,t));
	    (*J_nnz)++; // dcompZ/dQ

	    // H
	    MAT_set_i(Hy,H_nnz[*Jconstr_index],GEN_get_index_Q(rg,t));
	    MAT_set_j(Hy,H_nnz[*Jconstr_index],GEN_get_index_Q(rg,t));
	    H_nnz[*Jconstr_index]++;   // Q and Q

	    MAT_set_i(Hz,H_nnz[*Jconstr_index+1],GEN_get_index_Q(rg,t));
	    MAT_set_j(Hz,H_nnz[*Jconstr_index+1],GEN_get_index_Q(rg,t));
	    H_nnz[*Jconstr_index+1]++; // Q and Q

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {

		MAT_set_i(Hy,H_nnz[*Jconstr_index],GEN_get_index_Q(rg,t));
		MAT_set_j(Hy,H_nnz[*Jconstr_index],GEN_get_index_Q(rg1,t));
		H_nnz[*Jconstr_index]++;   // Q and Q1

		MAT_set_i(Hz,H_nnz[*Jconstr_index+1],GEN_get_index_Q(rg,t));
		MAT_set_j(Hz,H_nnz[*Jconstr_index+1],GEN_get_index_Q(rg1,t));
		H_nnz[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }

  // Done
  if ((t == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {

    // Ensure lower triangular and save struct of H comb
    H_nnz_comb = 0;
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
}

void CONSTR_REG_GEN_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {

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
  int* Jconstr_index;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL v;
  REAL v_set;
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
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!f || !J || !J_nnz || !Jconstr_index ||
      !H_nnz || !bus_counted)
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

      if (BUS_is_regulated_by_gen(bus) &&  // reg by gen
	  !BUS_is_slack(bus)) {            // not slack

	// Hessians
	Hy = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index));
	Hz = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+1));

	// Q value
	Qsum = 0;
	Qmax = 0;
	Qmin = 0;
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q))
	    Qsum += VEC_get(var_values,GEN_get_index_Q(rg,t)); // p.u.
	  else
	    Qsum += GEN_get_Q(rg,t);      // p.u.
	  Qmax += GEN_get_Q_max(rg); // p.u.
	  Qmin += GEN_get_Q_min(rg); // p.u.
	}
	Qy = (Qsum-Qmin);
	Qz = (Qmax-Qsum);

	// yz values
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	  y = VEC_get(var_values,BUS_get_index_y(bus,t)); // p.u.
	  z = VEC_get(var_values,BUS_get_index_z(bus,t));	// p.u.
	}
	else {
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	    v = VEC_get(var_values,BUS_get_index_v_mag(bus,t)); // p.u.
	  else
	    v = BUS_get_v_mag(bus,t);   // p.u.
	  v_set = BUS_get_v_set(bus,t); // p.u.
	  if (v > v_set) {
	    y = v-v_set;
	    z = 0;
	  }
	  else {
	    y = 0;
	    z = v_set-v;
	  }
	}

	// Terms
	sqrt_termY = sqrt( Qy*Qy + y*y + 2*CONSTR_REG_GEN_PARAM );
	sqrt_termZ = sqrt( Qz*Qz + z*z + 2*CONSTR_REG_GEN_PARAM );

	// f
	f[*Jconstr_index] = Qy + y - sqrt_termY;   // CompY
	f[*Jconstr_index+1] = Qz + z - sqrt_termZ; // CompZ

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {

	  // J
	  J[*J_nnz] = 1. - y/sqrt_termY;
	  (*J_nnz)++; // dCompY/dy

	  J[*J_nnz] = 1. - z/sqrt_termZ;
	  (*J_nnz)++; // dCompZ/dz

	  // H
	  Hy[H_nnz[*Jconstr_index]] = -(Qy*Qy+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	  H_nnz[*Jconstr_index]++;     // y and y

	  Hz[H_nnz[*Jconstr_index+1]] = -(Qz*Qz+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	  H_nnz[*Jconstr_index+1]++;   // z and z
	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {

	      Hy[H_nnz[*Jconstr_index]] = Qy*y/pow(sqrt_termY,3.);
	      H_nnz[*Jconstr_index]++;   // y and Q

	      Hz[H_nnz[*Jconstr_index+1]] = -Qz*z/pow(sqrt_termZ,3.);
	      H_nnz[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}

	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var

	    // J
	    J[*J_nnz] = 1. - Qy/sqrt_termY;
	    (*J_nnz)++; // dcompY/dQ

	    J[*J_nnz] = -1. + Qz/sqrt_termZ;
	    (*J_nnz)++; // dcompZ/dQ

	    // H
	    Hy[H_nnz[*Jconstr_index]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	    H_nnz[*Jconstr_index]++;   // Q and Q

	    Hz[H_nnz[*Jconstr_index+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	    H_nnz[*Jconstr_index+1]++; // Q and Q

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {

		Hy[H_nnz[*Jconstr_index]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
		H_nnz[*Jconstr_index]++;   // Q and Q1

		Hz[H_nnz[*Jconstr_index+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
		H_nnz[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
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
  int* Jconstr_index;
  char* bus_counted;
  int bus_index_t[2];
  REAL lamCompY;
  REAL lamCompZ;
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Jconstr_index || !bus_counted)
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

      if (BUS_is_regulated_by_gen(bus) && // reg gen
	  !BUS_is_slack(bus)) {           // not slack

	lamCompY = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // dCompY
	lamCompZ = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // dCompZ

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
