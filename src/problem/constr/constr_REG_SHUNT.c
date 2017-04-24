/** @file constr_REG_SHUNT.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_SHUNT.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_SHUNT.h>

Constr* CONSTR_REG_SHUNT_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c,&CONSTR_REG_SHUNT_init);
  CONSTR_set_func_count_step(c,&CONSTR_REG_SHUNT_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_REG_SHUNT_allocate);
  CONSTR_set_func_clear(c,&CONSTR_REG_SHUNT_clear);
  CONSTR_set_func_analyze_step(c,&CONSTR_REG_SHUNT_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_REG_SHUNT_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_REG_SHUNT_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_REG_SHUNT_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_REG_SHUNT_init(Constr* c) {

  // Local variables
  Net* net;
  int num_Jconstr;

  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 4*NET_get_num_switched_shunts(CONSTR_get_network(c))*NET_get_num_periods(net);
  CONSTR_set_H_nnz(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_name(c,"voltage regulation by shunts");
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_SHUNT_clear(Constr* c) {

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

void CONSTR_REG_SHUNT_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
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

      // Shunts
      //*******
      
      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {

	// Linear
	//*******
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // A
	  (*A_nnz)++; // b
	}
	
	(*A_nnz)++; // y
	(*A_nnz)++; // z
	
	(*A_row)++;

	// Nonlinear constraints 1 (vmax,vmin)
	//************************************

	// J
	(*J_nnz)++; // dCompVmin/dy
	(*J_nnz)++; // dCompVmax/dz
	
	// H
	H_nnz[*J_row]++;   // y and y (vmin)
	H_nnz[*J_row+1]++; // z and z (vmax)
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	  H_nnz[*J_row]++;   // y and v (vmin)
	  H_nnz[*J_row+1]++; // z and v (vmax)
	}
	H_nnz[*J_row]++;   // y and vl (vmin)
	H_nnz[*J_row+1]++; // z and vh (vmax)

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	  // J
	  (*J_nnz)++; // dCompVmin/dv
	  (*J_nnz)++; // dCompVmax/dv

	  // H
	  H_nnz[*J_row]++;   // v and v (vmin)
	  H_nnz[*J_row+1]++; // v and v (vmax)
	  H_nnz[*J_row]++;   // v and vl (vmin)
	  H_nnz[*J_row+1]++; // v and vh (vmax)
	}

	// J
	(*J_nnz)++; // dCompVmin/dvl
	(*J_nnz)++; // dCompVmax/dvh

	// H
	H_nnz[*J_row]++;   // vl and vl (vmin)
	H_nnz[*J_row+1]++; // vh and vh (vmax)

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  (*J_nnz)++; // dCompBmax/db
	  (*J_nnz)++; // dCompBmin/db

	  // H
	  H_nnz[*J_row+2]++; // b and b (bmax)
	  H_nnz[*J_row+3]++; // b and b (bmin)
	  H_nnz[*J_row+2]++; // b and vl (bmax)
	  H_nnz[*J_row+3]++; // b and vh (bmin)
	}

	// J
	(*J_nnz)++; // dCompBmax/dvl
	(*J_nnz)++; // dCompBmin/dvh

	// H
	H_nnz[*J_row+2]++; // vl and vl (bmax)
	H_nnz[*J_row+3]++; // vh and vh (bmin)

	// Count
	(*J_row)++; // CompVmin
	(*J_row)++; // CompVmax
	(*J_row)++; // CompBmax
	(*J_row)++; // CompBmin

	// Num extra vars
	CONSTR_set_num_extra_vars(c,*J_row);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_SHUNT_allocate(Constr* c) {

  // Local variables
  int A_nnz;
  int J_nnz;
  int A_row;
  int J_row;
  int* H_nnz;
  Mat* H_array;
  Mat* H;
  int H_comb_nnz;
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
  H_comb_nnz = 0;
  H_array = MAT_array_new(J_row);
  CONSTR_set_H_array(c,H_array,J_row);
  for (i = 0; i < J_row; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,H_nnz[i]);
    MAT_set_size1(H,num_vars+num_extra_vars);
    MAT_set_size2(H,num_vars+num_extra_vars);
    MAT_set_row_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)malloc(H_nnz[i]*sizeof(REAL)));
    H_comb_nnz += H_nnz[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars+num_extra_vars, // size1 (rows)
				  num_vars+num_extra_vars, // size2 (cols)
				  H_comb_nnz));            // nnz
}

void CONSTR_REG_SHUNT_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hvmin;
  Mat* Hvmax;
  Mat* Hbmin;
  Mat* Hbmax;
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;
  int H_nnz_comb;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int m;
  int temp;
  int index_v;
  int index_vl;
  int index_vh;
  int index_b;
  int index_y;
  int index_z;
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

  // Branch
  //*******

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Shunts
      //*******

      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {

	// Hessians (NOTE ORDER!!!)
	Hvmin = MAT_array_get(H_array,*J_row);
	Hvmax = MAT_array_get(H_array,*J_row+1);
	Hbmax = MAT_array_get(H_array,*J_row+2);
	Hbmin = MAT_array_get(H_array,*J_row+3);

	// Indices
	index_v = BUS_get_index_v_mag(bus,t);
	index_vl = num_vars+(*J_row+2);
	index_vh = num_vars+(*J_row+3);
	index_b = SHUNT_get_index_b(shunt,t);
	index_y = num_vars+(*J_row);
	index_z = num_vars+(*J_row+1);

	// Linear (b = b_0 + y - z)
	//*************************
	
	// b
	VEC_set(b,*A_row,SHUNT_get_b(shunt,t)); // current susceptance value

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // A
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,index_b);
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++; // b
	}
	else
	  VEC_add_to_entry(b,*A_row,-SHUNT_get_b(shunt,t));

	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,index_y);
	MAT_set_d(A,*A_nnz,-1.);
	(*A_nnz)++; // y
	
	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,index_z);
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++; // z

	(*A_row)++;

	// Nonlinear constraints 1 (vmin,vmax)
	//************************************

	// J
	MAT_set_i(J,*J_nnz,*J_row);
	MAT_set_j(J,*J_nnz,index_y);
	(*J_nnz)++; // dcompVmin/dy

	MAT_set_i(J,*J_nnz,*J_row+1);
	MAT_set_j(J,*J_nnz,index_z);
	(*J_nnz)++; // dcompVmax/dz
	
	// H
	MAT_set_i(Hvmin,H_nnz[*J_row],index_y);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_y);
	H_nnz[*J_row]++;   // y and y (vmin)
	
	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_z);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_z);
	H_nnz[*J_row+1]++; // z and z (vmax)
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

	  MAT_set_i(Hvmin,H_nnz[*J_row],index_y);
	  MAT_set_j(Hvmin,H_nnz[*J_row],index_v);
	  H_nnz[*J_row]++;   // y and v (vmin)
	  
	  MAT_set_i(Hvmax,H_nnz[*J_row+1],index_z);
	  MAT_set_j(Hvmax,H_nnz[*J_row+1],index_v);
	  H_nnz[*J_row+1]++; // z and v (vmax)
	}

	MAT_set_i(Hvmin,H_nnz[*J_row],index_y);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
	H_nnz[*J_row]++;   // y and vl (vmin)
	
	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_z);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
	H_nnz[*J_row+1]++; // z and vh (vmax)

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	  // J
	  MAT_set_i(J,*J_nnz,*J_row);
	  MAT_set_j(J,*J_nnz,index_v);
	  (*J_nnz)++; // dcompVmin/dv

	  MAT_set_i(J,*J_nnz,*J_row+1);
	  MAT_set_j(J,*J_nnz,index_v);
	  (*J_nnz)++; // dcompVmax/dv

	  // H
	  MAT_set_i(Hvmin,H_nnz[*J_row],index_v);
	  MAT_set_j(Hvmin,H_nnz[*J_row],index_v);
	  H_nnz[*J_row]++;   // v and v (vmin)

	  MAT_set_i(Hvmax,H_nnz[*J_row+1],index_v);
	  MAT_set_j(Hvmax,H_nnz[*J_row+1],index_v);
	  H_nnz[*J_row+1]++; // v and v (vmax)

	  MAT_set_i(Hvmin,H_nnz[*J_row],index_v);
	  MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
	  H_nnz[*J_row]++;   // v and vl (vmin)
	  
	  MAT_set_i(Hvmax,H_nnz[*J_row+1],index_v);
	  MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
	  H_nnz[*J_row+1]++; // v and vh (vmax)
	}

	// J
	MAT_set_i(J,*J_nnz,*J_row);
	MAT_set_j(J,*J_nnz,index_vl);
	(*J_nnz)++; // dcompVmin/dvl
	
	MAT_set_i(J,*J_nnz,*J_row+1);
	MAT_set_j(J,*J_nnz,index_vh);
	(*J_nnz)++; // dcompVmax/dvh
	
	// H
	MAT_set_i(Hvmin,H_nnz[*J_row],index_vl);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
	H_nnz[*J_row]++;   // vl and vl (vmin)
	
	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_vh);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
	H_nnz[*J_row+1]++; // vh and vh (vmax)

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  MAT_set_i(J,*J_nnz,*J_row+2);
	  MAT_set_j(J,*J_nnz,index_b);
	  (*J_nnz)++; // dcompBmax/db

	  MAT_set_i(J,*J_nnz,*J_row+3);
	  MAT_set_j(J,*J_nnz,index_b);
	  (*J_nnz)++; // dcompBmin/db

	  // H
	  MAT_set_i(Hbmax,H_nnz[*J_row+2],index_b);
	  MAT_set_j(Hbmax,H_nnz[*J_row+2],index_b);
	  H_nnz[*J_row+2]++; // b and b (bmax)

	  MAT_set_i(Hbmin,H_nnz[*J_row+3],index_b);
	  MAT_set_j(Hbmin,H_nnz[*J_row+3],index_b);
	  H_nnz[*J_row+3]++; // b and b (bmin)

	  MAT_set_i(Hbmax,H_nnz[*J_row+2],index_b);
	  MAT_set_j(Hbmax,H_nnz[*J_row+2],index_vl);
	  H_nnz[*J_row+2]++; // b and vl (bmax)
	  
	  MAT_set_i(Hbmin,H_nnz[*J_row+3],index_b);
	  MAT_set_j(Hbmin,H_nnz[*J_row+3],index_vh);
	  H_nnz[*J_row+3]++; // b and vh (bmin)
	}

	// J
	MAT_set_i(J,*J_nnz,*J_row+2);
	MAT_set_j(J,*J_nnz,index_vl);
	(*J_nnz)++; // dcompBmax/dvl
	
	MAT_set_i(J,*J_nnz,*J_row+3);
	MAT_set_j(J,*J_nnz,index_vh);
	(*J_nnz)++; // dcompBmin/dvh
	
	// H
	MAT_set_i(Hbmax,H_nnz[*J_row+2],index_vl);
	MAT_set_j(Hbmax,H_nnz[*J_row+2],index_vl);
	H_nnz[*J_row+2]++; // vl and vl (bmax)
	
	MAT_set_i(Hbmin,H_nnz[*J_row+3],index_vh);
	MAT_set_j(Hbmin,H_nnz[*J_row+3],index_vh);
	H_nnz[*J_row+3]++; // vh and vh (bmin)

	// Count
	(*J_row)++; // compVmin
	(*J_row)++; // compVmax
	(*J_row)++; // compBmax
	(*J_row)++; // compBmin
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

void CONSTR_REG_SHUNT_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* Hvmin;
  REAL* Hvmax;
  REAL* Hbmin;
  REAL* Hbmax;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL v;
  REAL vl;
  REAL vh;
  REAL vmin;
  REAL vmax;
  REAL b;
  REAL bmax;
  REAL bmin;
  REAL y;
  REAL z;
  REAL sqrtermVmin;
  REAL sqrtermVmax;
  REAL sqrtermBmax;
  REAL sqrtermBmin;
  REAL norm = CONSTR_REG_SHUNT_NORM;
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

      // Shunts
      //*******

      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {

	// v values
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	  v = VEC_get(values,BUS_get_index_v_mag(bus,t));
	else
	  v = BUS_get_v_mag(bus,t);
	vmax = BUS_get_v_max(bus);
	vmin = BUS_get_v_min(bus);
	if (VEC_get_size(values_extra) > 0) {
	  vl = VEC_get(values_extra,*J_row+2);
	  vh = VEC_get(values_extra,*J_row+3);
	}
	else {
	  vl = 0;
	  vh = 0;
	}

	// b values
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  b = VEC_get(values,SHUNT_get_index_b(shunt,t));
	}
	else
	  b = SHUNT_get_b(shunt,t);
	bmax = SHUNT_get_b_max(shunt);
	bmin = SHUNT_get_b_min(shunt);
	if (VEC_get_size(values_extra) > 0) {
	  y = VEC_get(values_extra,*J_row);
	  z = VEC_get(values_extra,*J_row+1);
	}
	else {
	  y = 0.;
	  z = 0;
	}

	// Terms
	sqrtermVmin = sqrt( (v+vl-vmin)*(v+vl-vmin) + y*y + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermVmax = sqrt( (vmax-v+vh)*(vmax-v+vh) + z*z + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermBmax = sqrt( (bmax-b)*(bmax-b) + vl*vl + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermBmin = sqrt( (b-bmin)*(b-bmin) + vh*vh + 2*CONSTR_REG_SHUNT_PARAM );

	// Hessians (NOTE ORDER!!!)
	Hvmin = MAT_get_data_array(MAT_array_get(H_array,*J_row));
	Hvmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));
	Hbmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+2));
	Hbmin = MAT_get_data_array(MAT_array_get(H_array,*J_row+3));

	// f
	f[*J_row] = ((v+vl-vmin) + y - sqrtermVmin)*norm;   // vmin
	f[*J_row+1] = ((vmax-v+vh) + z - sqrtermVmax)*norm; // vmax
	f[*J_row+2] = ((bmax-b) + vl - sqrtermBmax)*norm;   // bmax
	f[*J_row+3] = ((b-bmin) + vh - sqrtermBmin)*norm;   // bmin

	// Nonlinear constraints 1 (vmin,vmax)
	//************************************

	// J
	J[*J_nnz] = (1.-y/sqrtermVmin)*norm;
	(*J_nnz)++; // dcompVmin/dy
	
	J[*J_nnz] = (1.-z/sqrtermVmax)*norm;
	(*J_nnz)++; // dcompVmax/dz

	// H
	Hvmin[H_nnz[*J_row]] = -(((v+vl-vmin)*(v+vl-vmin)+2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // y and y (vmin)
	
	Hvmax[H_nnz[*J_row+1]] = -(((vmax-v+vh)*(vmax-v+vh)+2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // z and z (vmax)

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

	  Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*y/pow(sqrtermVmin,3.))*norm;
	  H_nnz[*J_row]++;   // y and v (vmin)
	  
	  Hvmax[H_nnz[*J_row+1]] = -((vmax-v+vh)*z/pow(sqrtermVmax,3.))*norm;
	  H_nnz[*J_row+1]++; // z and v (vmax)
	}

	Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*y/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // y and vl (vmin)
	
	Hvmax[H_nnz[*J_row+1]] = ((vmax-v+vh)*z/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // z and vh (vmax)

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	  // J
	  J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
	  (*J_nnz)++; // dcompVmin/dv
	  
	  J[*J_nnz] = -((1.-(vmax-v+vh)/sqrtermVmax))*norm;
	  (*J_nnz)++; // dcompVmax/dv
	  
	  // H
	  Hvmin[H_nnz[*J_row]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	  H_nnz[*J_row]++;   // v and v (vmin)

	  Hvmax[H_nnz[*J_row+1]] = -((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	  H_nnz[*J_row+1]++; // v and v (vmax)
	  
	  Hvmin[H_nnz[*J_row]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	  H_nnz[*J_row]++;   // v and vl (vmin)
	  
	  Hvmax[H_nnz[*J_row+1]] = ((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	  H_nnz[*J_row+1]++; // v and vh (vmax)
	}

	// J
	J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
	(*J_nnz)++; // dcompVmin/dvl
	
	J[*J_nnz] = (1.-(vmax-v+vh)/sqrtermVmax)*norm;
	(*J_nnz)++; // dcompVmax/dvh
	
	// H
	Hvmin[H_nnz[*J_row]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // vl and vl (vmin)
	
	Hvmax[H_nnz[*J_row+1]] = -((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // vh and vh (vmax)

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // t var

	  // J
	  J[*J_nnz] = -(1.-(bmax-b)/sqrtermBmax)*norm;
	  (*J_nnz)++; // dcompBmax/db

	  J[*J_nnz] = (1.-(b-bmin)/sqrtermBmin)*norm;
	  (*J_nnz)++; // dcompBmin/db

	  // H
	  Hbmax[H_nnz[*J_row+2]] = -((vl*vl + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmax,3.))*norm;
	  H_nnz[*J_row+2]++; // b and b (bmax)

	  Hbmin[H_nnz[*J_row+3]] = -((vh*vh + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmin,3.))*norm;
	  H_nnz[*J_row+3]++; // b and b (bmin)
	  
	  Hbmax[H_nnz[*J_row+2]] = -(vl*(bmax-b)/pow(sqrtermBmax,3.))*norm;
	  H_nnz[*J_row+2]++; // b and vl (bmax)
	  
	  Hbmin[H_nnz[*J_row+3]] = (vh*(b-bmin)/pow(sqrtermBmin,3.))*norm;
	  H_nnz[*J_row+3]++; // b and vh (bmin)
	}

	// J
	J[*J_nnz] = (1.-vl/sqrtermBmax)*norm;
	(*J_nnz)++; // dcompBmax/dvl
	
	J[*J_nnz] = (1.-vh/sqrtermBmin)*norm;
	(*J_nnz)++; // dcompBmin/dvh
	
	// H
	Hbmax[H_nnz[*J_row+2]] = -(((bmax-b)*(bmax-b) + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmax,3.))*norm;
	H_nnz[*J_row+2]++; // vl and vl (bmax)
	
	Hbmin[H_nnz[*J_row+3]] = -(((b-bmin)*(b-bmin) + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmin,3.))*norm;
	H_nnz[*J_row+3]++; // vh and vh (bmin)

	// Count
	(*J_row)++; // compVmin
	(*J_row)++; // compVmax
	(*J_row)++; // compBmax
	(*J_row)++; // compBmin
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_SHUNT_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  int* J_row;
  char* bus_counted;
  int bus_index_t[2];
  REAL lamCompVmin;
  REAL lamCompVmax;
  REAL lamCompBmax;
  REAL lamCompBmin;
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
      
      // Shunts
      //*******

      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {

	lamCompVmin = VEC_get(sf,*J_row);
	(*J_row)++; // compVmin
	lamCompVmax = VEC_get(sf,*J_row);
	(*J_row)++; // compVmax
	lamCompBmax = VEC_get(sf,*J_row);
	(*J_row)++; // compBmax
	lamCompBmin = VEC_get(sf,*J_row);
	(*J_row)++; // compBmin

	if (fabs(lamCompVmin) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompVmin,t);
	if (fabs(lamCompVmax) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompVmax,t);
	if (fabs(lamCompBmax) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompBmax,t);
	if (fabs(lamCompBmin) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompBmin,t);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_SHUNT_free(Constr* c) {
  // Nothing
}
