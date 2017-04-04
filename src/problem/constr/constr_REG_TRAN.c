/** @file constr_REG_TRAN.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_TRAN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_TRAN.h>

Constr* CONSTR_REG_TRAN_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_REG_TRAN_init);
  CONSTR_set_func_count_step(c, &CONSTR_REG_TRAN_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_REG_TRAN_allocate);
  CONSTR_set_func_clear(c, &CONSTR_REG_TRAN_clear);
  CONSTR_set_func_analyze_step(c, &CONSTR_REG_TRAN_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_REG_TRAN_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_REG_TRAN_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_REG_TRAN_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_REG_TRAN_init(Constr* c) {

  // Local variables
  Net* net;
  int num_Jconstr;

  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 4*NET_get_num_tap_changers_v(CONSTR_get_network(c))*NET_get_num_periods(net);;
  CONSTR_set_H_nnz(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_name(c,"voltage regulation by transformers");
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_TRAN_clear(Constr* c) {

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
}

void CONSTR_REG_TRAN_count_step(Constr* c, Branch* br, int tau) {

  // Local variables
  Bus* reg_bus;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row ||
      !J_row || !H_nnz)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);

    // Linear
    //*******
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&     // t var
	BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // A
      (*A_nnz)++; // t
      (*A_nnz)++; // y
      (*A_nnz)++; // z

      (*A_row)++;
    }

    // Nonlinear constraints 1 (vmax,vmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // J
      (*J_nnz)++; // dcompVmin/dy
      (*J_nnz)++; // dcompVmax/dz

      // H
      H_nnz[*J_row]++;   // y and y (vmin)
      H_nnz[*J_row+1]++; // z and z (vmax)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	H_nnz[*J_row]++;   // y and v (vmin)
	H_nnz[*J_row+1]++; // z and v (vmax)
      }
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	H_nnz[*J_row]++;   // y and vl (vmin)
	H_nnz[*J_row+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

      // J
      (*J_nnz)++; // dcompVmin/dv
      (*J_nnz)++; // dcompVmax/dv

      // H
      H_nnz[*J_row]++;   // v and v (vmin)
      H_nnz[*J_row+1]++; // v and v (vmax)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	H_nnz[*J_row]++;   // v and vl (vmin)
	H_nnz[*J_row+1]++; // v and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      (*J_nnz)++; // dcompVmin/dvl
      (*J_nnz)++; // dcompVmax/dvh

      // H
      H_nnz[*J_row]++;   // vl and vl (vmin)
      H_nnz[*J_row+1]++; // vh and vh (vmax)
    }

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var

      // J
      (*J_nnz)++; // dcompTmax/dt
      (*J_nnz)++; // dcompTmin/dt

      // H
      H_nnz[*J_row+2]++; // t and t (tmax)
      H_nnz[*J_row+3]++; // t and t (tmin)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	H_nnz[*J_row+2]++; // t and vl (tmax)
	H_nnz[*J_row+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      (*J_nnz)++; // dcompTmax/dvl
      (*J_nnz)++; // dcompTmin/dvh

      // H
      H_nnz[*J_row+2]++; // vl and vl (tmax)
      H_nnz[*J_row+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*J_row)++; // compVmin
    (*J_row)++; // compVmax
    (*J_row)++; // compTmax
    (*J_row)++; // compTmin
  }
}

void CONSTR_REG_TRAN_allocate(Constr* c) {

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
  int i;

  A_nnz = CONSTR_get_A_nnz(c);
  J_nnz = CONSTR_get_J_nnz(c);
  A_row = CONSTR_get_A_row(c);
  J_row = CONSTR_get_J_row(c);
  H_nnz = CONSTR_get_H_nnz(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // b
  CONSTR_set_b(c,VEC_new(A_row));

  // A
  CONSTR_set_A(c,MAT_new(A_row, // size1 (rows)
			 num_vars,      // size2 (cols)
			 A_nnz));    // nnz

  // f
  CONSTR_set_f(c,VEC_new(J_row));

  // J
  CONSTR_set_J(c,MAT_new(J_row, // size1 (rows)
			 num_vars,      // size2 (cols)
			 J_nnz));    // nnz

  // H
  H_comb_nnz = 0;
  H_array = MAT_array_new(J_row);
  CONSTR_set_H_array(c,H_array,J_row);
  for (i = 0; i < J_row; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,H_nnz[i]);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(H_nnz[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)malloc(H_nnz[i]*sizeof(REAL)));
    H_comb_nnz += H_nnz[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
}

void CONSTR_REG_TRAN_analyze_step(Constr* c, Branch* br, int tau) {

  // Local variables
  Bus* reg_bus;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hvmin;
  Mat* Hvmax;
  Mat* Htmin;
  Mat* Htmax;
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
  int k;
  int m;
  int temp;
  int index_yz_vmin;
  int index_yz_vmax;
  int index_vvio_tmax;
  int index_vvio_tmin;
  int index_v;
  int index_vl;
  int index_vh;
  int index_t;
  int index_y;
  int index_z;
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
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row ||
      !J_row || !H_nnz)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);

    // Hessians (NOTE ORDER!!!)
    Hvmin = MAT_array_get(H_array,*J_row);
    Hvmax = MAT_array_get(H_array,*J_row+1);
    Htmax = MAT_array_get(H_array,*J_row+2);
    Htmin = MAT_array_get(H_array,*J_row+3);

    // Indices
    if (BRANCH_has_pos_ratio_v_sens(br)) {
      index_yz_vmin = BRANCH_get_index_ratio_y(br,tau);
      index_yz_vmax = BRANCH_get_index_ratio_z(br,tau);
      index_vvio_tmax = BUS_get_index_vl(reg_bus,tau);
      index_vvio_tmin = BUS_get_index_vh(reg_bus,tau);
    }
    else {
      index_yz_vmin = BRANCH_get_index_ratio_z(br,tau);
      index_yz_vmax = BRANCH_get_index_ratio_y(br,tau);
      index_vvio_tmax = BUS_get_index_vh(reg_bus,tau);
      index_vvio_tmin = BUS_get_index_vl(reg_bus,tau);
    }
    index_v = BUS_get_index_v_mag(reg_bus,tau);
    index_vl = BUS_get_index_vl(reg_bus,tau);
    index_vh = BUS_get_index_vh(reg_bus,tau);
    index_t = BRANCH_get_index_ratio(br,tau);
    index_y = BRANCH_get_index_ratio_y(br,tau);
    index_z = BRANCH_get_index_ratio_z(br,tau);

    // Linear (t = t_0 + y - z)
    //*************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&     // t var
	BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // b
      VEC_set(b,*A_row,BRANCH_get_ratio(br,tau)); // current ratio value

      // A
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_t);
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++; // t

      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_y);
      MAT_set_d(A,*A_nnz,-1.);
      (*A_nnz)++; // y

      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_z);
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++; // z

      (*A_row)++;
    }

    // Nonlinear constraints 1 (vmin,vmax)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,index_yz_vmin);
      (*J_nnz)++; // dcompVmin/dy

      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,index_yz_vmax);
      (*J_nnz)++; // dcompVmax/dz

      // H
      MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
      MAT_set_j(Hvmin,H_nnz[*J_row],index_yz_vmin);
      H_nnz[*J_row]++;   // y and y (vmin)

      MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
      MAT_set_j(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
      H_nnz[*J_row+1]++; // z and z (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {

	MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_v);
	H_nnz[*J_row]++;   // y and v (vmin)

	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_v);
	H_nnz[*J_row+1]++; // z and v (vmax)
      }

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	MAT_set_i(Hvmin,H_nnz[*J_row],index_yz_vmin);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
	H_nnz[*J_row]++;   // y and vl (vmin)

	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_yz_vmax);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vl);
	H_nnz[*J_row+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

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

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	MAT_set_i(Hvmin,H_nnz[*J_row],index_v);
	MAT_set_j(Hvmin,H_nnz[*J_row],index_vl);
	H_nnz[*J_row]++;   // v and vl (vmin)

	MAT_set_i(Hvmax,H_nnz[*J_row+1],index_v);
	MAT_set_j(Hvmax,H_nnz[*J_row+1],index_vh);
	H_nnz[*J_row+1]++; // v and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

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
    }

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var

      // J
      MAT_set_i(J,*J_nnz,*J_row+2);
      MAT_set_j(J,*J_nnz,index_t);
      (*J_nnz)++; // dcompTmax/dt

      MAT_set_i(J,*J_nnz,*J_row+3);
      MAT_set_j(J,*J_nnz,index_t);
      (*J_nnz)++; // dcompTmin/dt

      // H
      MAT_set_i(Htmax,H_nnz[*J_row+2],index_t);
      MAT_set_j(Htmax,H_nnz[*J_row+2],index_t);
      H_nnz[*J_row+2]++; // t and t (tmax)

      MAT_set_i(Htmin,H_nnz[*J_row+3],index_t);
      MAT_set_j(Htmin,H_nnz[*J_row+3],index_t);
      H_nnz[*J_row+3]++; // t and t (tmin)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	MAT_set_i(Htmax,H_nnz[*J_row+2],index_t);
	MAT_set_j(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
	H_nnz[*J_row+2]++; // t and vl (tmax)

	MAT_set_i(Htmin,H_nnz[*J_row+3],index_t);
	MAT_set_j(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
	H_nnz[*J_row+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      MAT_set_i(J,*J_nnz,*J_row+2);
      MAT_set_j(J,*J_nnz,index_vvio_tmax);
      (*J_nnz)++; // dcompTmax/dvl

      MAT_set_i(J,*J_nnz,*J_row+3);
      MAT_set_j(J,*J_nnz,index_vvio_tmin);
      (*J_nnz)++; // dcompTmin/dvh

      // H
      MAT_set_i(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
      MAT_set_j(Htmax,H_nnz[*J_row+2],index_vvio_tmax);
      H_nnz[*J_row+2]++; // vl and vl (tmax)

      MAT_set_i(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
      MAT_set_j(Htmin,H_nnz[*J_row+3],index_vvio_tmin);
      H_nnz[*J_row+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*J_row)++; // compVmin
    (*J_row)++; // compVmax
    (*J_row)++; // compTmax
    (*J_row)++; // compTmin
  }

  // Done
  if ((tau == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {

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

void CONSTR_REG_TRAN_eval_step(Constr* c, Branch* br, int tau, Vec* values) {

  // Local variables
  Bus* reg_bus;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* Hvmin;
  REAL* Hvmax;
  REAL* Htmin;
  REAL* Htmax;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL v;
  REAL vl;
  REAL vh;
  REAL vmin;
  REAL vmax;
  REAL t;
  REAL tmax;
  REAL tmin;
  REAL y;
  REAL z;
  REAL yz_vmin;
  REAL yz_vmax;
  REAL vvio_tmax;
  REAL vvio_tmin;
  REAL sqrtermVmin;
  REAL sqrtermVmax;
  REAL sqrtermTmax;
  REAL sqrtermTmin;
  REAL norm = CONSTR_REG_TRAN_NORM;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!f || !J || !J_nnz ||
      !J_row || !H_nnz)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);

    // v values
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG))
      v = VEC_get(values,BUS_get_index_v_mag(reg_bus,tau));
    else
      v = BUS_get_v_mag(reg_bus,tau);
    vmax = BUS_get_v_max_reg(reg_bus);
    vmin = BUS_get_v_min_reg(reg_bus);
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
      vl = VEC_get(values,BUS_get_index_vl(reg_bus,tau));
      vh = VEC_get(values,BUS_get_index_vh(reg_bus,tau));
    }
    else {
      vl = 0;
      vh = 0;
    }

    // t values
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
      t = VEC_get(values,BRANCH_get_index_ratio(br,tau));
    }
    else
      t = BRANCH_get_ratio(br,tau);
    tmax = BRANCH_get_ratio_max(br);
    tmin = BRANCH_get_ratio_min(br);
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) {
      y = VEC_get(values,BRANCH_get_index_ratio_y(br,tau));
      z = VEC_get(values,BRANCH_get_index_ratio_z(br,tau));
    }
    else {
      if (t > BRANCH_get_ratio(br,tau)) {
	y = t-BRANCH_get_ratio(br,tau);
	z = 0;
      }
      else {
	y = 0;
	z = BRANCH_get_ratio(br,tau)-t;
      }
    }

    // values that depend on sensitivity
    if (BRANCH_has_pos_ratio_v_sens(br)) {
      yz_vmin = y;
      yz_vmax = z;
      vvio_tmax = vl;
      vvio_tmin = vh;
    }
    else {
      yz_vmin = z;
      yz_vmax = y;
      vvio_tmax = vh;
      vvio_tmin = vl;
    }

    // Terms
    sqrtermVmin = sqrt( (v+vl-vmin)*(v+vl-vmin) + yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermVmax = sqrt( (vmax-v+vh)*(vmax-v+vh) + yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermTmax = sqrt( (tmax-t)*(tmax-t) + vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermTmin = sqrt( (t-tmin)*(t-tmin) + vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM );

    // Hessians (NOTE ORDER!!!)
    Hvmin = MAT_get_data_array(MAT_array_get(H_array,*J_row));
    Hvmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));
    Htmax = MAT_get_data_array(MAT_array_get(H_array,*J_row+2));
    Htmin = MAT_get_data_array(MAT_array_get(H_array,*J_row+3));

    // f
    f[*J_row] = ((v+vl-vmin) + yz_vmin - sqrtermVmin)*norm;   // vmin
    f[*J_row+1] = ((vmax-v+vh) + yz_vmax - sqrtermVmax)*norm; // vmax
    f[*J_row+2] = ((tmax-t) + vvio_tmax - sqrtermTmax)*norm;  // tmax
    f[*J_row+3] = ((t-tmin) + vvio_tmin - sqrtermTmin)*norm;  // tmin

    // Nonlinear constraints 1 (vmin,vmax)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // J
      J[*J_nnz] = (1.-yz_vmin/sqrtermVmin)*norm;
      (*J_nnz)++; // dcompVmin/dy

      J[*J_nnz] = (1.-yz_vmax/sqrtermVmax)*norm;
      (*J_nnz)++; // dcompVmax/dz

      // H
      Hvmin[H_nnz[*J_row]] = -(((v+vl-vmin)*(v+vl-vmin)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // y and y (vmin)

      Hvmax[H_nnz[*J_row+1]] = -(((vmax-v+vh)*(vmax-v+vh)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // z and z (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {

	Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // y and v (vmin)

	Hvmax[H_nnz[*J_row+1]] = -((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // z and v (vmax)
      }

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	Hvmin[H_nnz[*J_row]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // y and vl (vmin)

	Hvmax[H_nnz[*J_row+1]] = ((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

      // J
      J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
      (*J_nnz)++; // dcompVmin/dv

      J[*J_nnz] = -((1.-(vmax-v+vh)/sqrtermVmax))*norm;
      (*J_nnz)++; // dcompVmax/dv

      // H
      Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // v and v (vmin)

      Hvmax[H_nnz[*J_row+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // v and v (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
	H_nnz[*J_row]++;   // v and vl (vmin)

	Hvmax[H_nnz[*J_row+1]] = ((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
	H_nnz[*J_row+1]++; // v and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      J[*J_nnz] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
      (*J_nnz)++; // dcompVmin/dvl

      J[*J_nnz] = (1.-(vmax-v+vh)/sqrtermVmax)*norm;
      (*J_nnz)++; // dcompVmax/dvh

      // H
      Hvmin[H_nnz[*J_row]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      H_nnz[*J_row]++;   // vl and vl (vmin)

      Hvmax[H_nnz[*J_row+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      H_nnz[*J_row+1]++; // vh and vh (vmax)
    }

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var

      // J
      J[*J_nnz] = -(1.-(tmax-t)/sqrtermTmax)*norm;
      (*J_nnz)++; // dcompTmax/dt

      J[*J_nnz] = (1.-(t-tmin)/sqrtermTmin)*norm;
      (*J_nnz)++; // dcompTmin/dt

      // H
      Htmax[H_nnz[*J_row+2]] = -((vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
      H_nnz[*J_row+2]++; // t and t (tmax)

      Htmin[H_nnz[*J_row+3]] = -((vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
      H_nnz[*J_row+3]++; // t and t (tmin)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	Htmax[H_nnz[*J_row+2]] = -(vvio_tmax*(tmax-t)/pow(sqrtermTmax,3.))*norm;
	H_nnz[*J_row+2]++; // t and vl (tmax)

	Htmin[H_nnz[*J_row+3]] = (vvio_tmin*(t-tmin)/pow(sqrtermTmin,3.))*norm;
	H_nnz[*J_row+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      J[*J_nnz] = (1.-vvio_tmax/sqrtermTmax)*norm;
      (*J_nnz)++; // dcompTmax/dvl

      J[*J_nnz] = (1.-vvio_tmin/sqrtermTmin)*norm;
      (*J_nnz)++; // dcompTmin/dvh

      // H
      Htmax[H_nnz[*J_row+2]] = -(((tmax-t)*(tmax-t) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
      H_nnz[*J_row+2]++; // vl and vl (tmax)

      Htmin[H_nnz[*J_row+3]] = -(((t-tmin)*(t-tmin) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
      H_nnz[*J_row+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*J_row)++; // compVmin
    (*J_row)++; // compVmax
    (*J_row)++; // compTmax
    (*J_row)++; // compTmin
  }
}

void CONSTR_REG_TRAN_store_sens_step(Constr* c, Branch* br, int tau, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* reg_bus;
  int* J_row;
  REAL lamCompVmin;
  REAL lamCompVmax;
  REAL lamCompTmax;
  REAL lamCompTmin;

  // Constr data
  J_row = CONSTR_get_J_row_ptr(c);

  // Check pointer
  if (!J_row)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);

    lamCompVmin = VEC_get(sf,*J_row);
    (*J_row)++; // compVmin
    lamCompVmax = VEC_get(sf,*J_row);
    (*J_row)++; // compVmax
    lamCompTmax = VEC_get(sf,*J_row);
    (*J_row)++; // compTmax
    lamCompTmin = VEC_get(sf,*J_row);
    (*J_row)++; // compTmin

    if (fabs(lamCompVmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmin,tau);
    if (fabs(lamCompVmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmax,tau);
    if (fabs(lamCompTmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmax,tau);
    if (fabs(lamCompTmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmin,tau);
  }
}

void CONSTR_REG_TRAN_free(Constr* c) {
  // Nothing
}
