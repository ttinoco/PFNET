/** @file constr_FACTS_EQ.c
 *  @brief This file defines the data structure and routines associated with the constraint of type FACTS_EQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_FACTS_EQ.h>

Constr* CONSTR_FACTS_EQ_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_FACTS_EQ_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_FACTS_EQ_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_FACTS_EQ_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_FACTS_EQ_store_sens_step);
  CONSTR_set_name(c,"FACTS equations");
  return c;
}

void CONSTR_FACTS_EQ_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Bus* bus_k;
  Bus* bus_m;
  Facts* facts;
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
  if (!A_nnz || !J_nnz || !A_row || !J_row || !H_nnz || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {

    // Facts which are not regulating voltage
    if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_Q) && !FACTS_is_regulator(facts)) {

      // Q_sh = 0, because i_max_sh = 0

      // A
      (*A_nnz)++; // Qsh

      // Count row
      (*A_row)++;
    }

    if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_P | FACTS_VAR_Q)) {

      // Linear Power balance equations
      // P_k + P_m - P_dc = 0
      // Q_k + Q_m - Q_sh - Qs = 0

      // A
      (*A_nnz)++; // P_k
      (*A_nnz)++; // P_m
      (*A_nnz)++; // P_dc

      // Count row
      (*A_row)++;

      // A
      (*A_nnz)++; // Q_k
      (*A_nnz)++; // Q_m
      (*A_nnz)++; // Q_sh
      (*A_nnz)++; // Q_s
          
      // Count row
      (*A_row)++;

      if (FACTS_get_P_max_dc(facts) == 0 || FACTS_is_series_link_disabled(facts)) {

        // P_dc = 0

        // A
        (*A_nnz)++; // P_dc

        // Count row
        (*A_row)++;

      }

      if (FACTS_is_series_link_disabled(facts) &&
          FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {

        // linear equations
        // P_m = 0, Q_m = 0, Q_s = 0, vmag_s = 0, vang_s = 0

        (*A_nnz)++; // P_m

        // Count row
        (*A_row)++;

        (*A_nnz)++; // Q_m

        // Count row
        (*A_row)++;

        (*A_nnz)++; // Q_s

        // Count row
        (*A_row)++;

        (*A_nnz)++; // vmag_s

        // Count row
        (*A_row)++;

        (*A_nnz)++; // vang_s

        // Count row
        (*A_row)++;

      }

      if (FACTS_is_in_normal_series_mode(facts)) {

        // Get the buses connected to the FACTS device
        bus_k = FACTS_get_bus_k(facts);
        bus_m = FACTS_get_bus_m(facts);

        // Non linear voltage balance equations
        // -vmag_k*cos(vang_k)+vmag_m*cos(vang_m)-vmag_s*cos(vang_s) = 0 (f1)
        // -vmag_k*sin(vang_k)+vmag_m*sin(vang_m)-vmag_s*sin(vang_s) = 0 (f2)
        // Non linear series converter power balance equations
        // vmag_s*P_m*cos(vang_s)-vmag_s*Q_m*sin(vang_s)-vmag_m*P_dc*cos(vang_m)+vmag_m*Q_s*sin(vang_m) = 0 (f3)
        // vmag_s*P_m*sin(vang_s)+vmag_s*Q_m*cos(vang_s)-vmag_m*P_dc*sin(vang_m)-vmag_m*Q_s*cos(vang_m) = 0 (f4)

        (*J_nnz)++; // df3/P_m
        (*J_nnz)++; // df3/Q_m
        (*J_nnz)++; // df3/P_dc
        (*J_nnz)++; // df3/Q_s
        (*J_nnz)++; // df4/P_m
        (*J_nnz)++; // df4/Q_m
        (*J_nnz)++; // df4/P_dc
        (*J_nnz)++; // df4/Q_s

        if (BUS_has_flags(bus_k, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          (*J_nnz)++; // df1/dvmag_k
          (*J_nnz)++; // df1/vang_k
          (*J_nnz)++; // df2/dvmag_k
          (*J_nnz)++; // df2/vang_k

          // H
          H_nnz[*J_row]++; // df1/(dvang_k and dvmag_k)
          H_nnz[*J_row]++; // df1/(dvang_k and dvang_k)
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvmag_k)
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvang_k)

        }

        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          (*J_nnz)++; // df1/dvmag_m
          (*J_nnz)++; // df1/vang_m
          (*J_nnz)++; // df2/dvmag_m
          (*J_nnz)++; // df2/vang_m
          (*J_nnz)++; // df3/dvmag_m
          (*J_nnz)++; // df3/vang_m
          (*J_nnz)++; // df4/dvmag_m
          (*J_nnz)++; // df4/vang_m

          // H
          H_nnz[*J_row]++; // df1/(dvang_m and dvmag_m)
          H_nnz[*J_row]++; // df1/(dvang_m and dvang_m)
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvmag_m)
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvang_m)
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvmag_m)
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvmag_m)
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvmag_m)
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvang_m)
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvang_m)
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvang_m)
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvmag_m)
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvmag_m)
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvmag_m)
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvang_m)
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvang_m)
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvang_m)

        }

        if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {

          // J
          (*J_nnz)++; // df1/dvmag_s
          (*J_nnz)++; // df1/vang_s
          (*J_nnz)++; // df2/dvmag_s
          (*J_nnz)++; // df2/vang_s
          (*J_nnz)++; // df3/dvmag_s
          (*J_nnz)++; // df3/vang_s
          (*J_nnz)++; // df4/dvmag_s
          (*J_nnz)++; // df4/vang_s

          // H
          H_nnz[*J_row]++; // df1/(dvang_s and dvmag_s)
          H_nnz[*J_row]++; // df1/(dvang_s and dvang_s)
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvmag_s)
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvang_s)
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvmag_s)
          H_nnz[*J_row+2]++; // df3/(dP_m and dvmag_s)
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvmag_s)
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvang_m)
          H_nnz[*J_row+2]++; // df3/(dP_m and dvang_s)
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvang_s)
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvmag_s)
          H_nnz[*J_row+3]++; // df4/(dP_m and dvmag_s)
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvmag_s)
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvang_s)
          H_nnz[*J_row+3]++; // df4/(dP_m and dvang_s)
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvang_s)

        }

        // Count rows
        (*J_row)++; // df1
        (*J_row)++; // df2
        (*J_row)++; // df3
        (*J_row)++; // df4
      }
    }
  }
}

void CONSTR_FACTS_EQ_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Bus* bus_k;
  Bus* bus_m;
  Facts* facts;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* H1;
  Mat* H2;
  Mat* H3;
  Mat* H4;
  int* A_nnz;
  int* J_nnz;
  int* A_row;
  int* J_row;
  int* H_nnz;

  // Constr data
  A = CONSTR_get_A(c);
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !J_row || !H_nnz || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {

    // Facts which are not regulating voltage
    if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_Q) && !FACTS_is_regulator(facts)) {

      // Q_sh = 0, because i_max_sh = 0

      // A
      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_Q_sh(facts, t));
      MAT_set_d(A, *A_nnz, 1.);

      (*A_nnz)++; // Qsh

      // Count row
      (*A_row)++;
    }

    if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_P | FACTS_VAR_Q)) {

      // Linear Power balance equations
      // P_k + P_m - P_dc = 0
      // Q_k + Q_m - Q_sh - Qs = 0

      // A
      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_P_k(facts, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // P_k

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_P_m(facts, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // P_m

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_P_dc(facts, t));
      MAT_set_d(A, *A_nnz, -1.);
      (*A_nnz)++; // P_dc

      // Count row
      (*A_row)++;

      // A
      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_Q_k(facts, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // Q_k

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_Q_m(facts, t));
      MAT_set_d(A, *A_nnz, 1.);
      (*A_nnz)++; // Q_m

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_Q_sh(facts, t));
      MAT_set_d(A, *A_nnz, -1.);
      (*A_nnz)++; // Q_sh

      MAT_set_i(A, *A_nnz, *A_row);
      MAT_set_j(A, *A_nnz, FACTS_get_index_Q_s(facts, t));
      MAT_set_d(A, *A_nnz, -1.);
      (*A_nnz)++; // Q_s

      // Count row
      (*A_row)++;

      if (FACTS_get_P_max_dc(facts) == 0 || FACTS_is_series_link_disabled(facts)) {

        // P_dc = 0

        // A
        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_P_dc(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // P_dc

        // Count row
        (*A_row)++;

      }

      if (FACTS_is_series_link_disabled(facts) &&
          FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {

        // linear equations
        // P_m = 0, Q_m = 0, Q_s = 0, vmag_s = 0, vang_s = 0

        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_P_m(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // P_m

        // Count row
        (*A_row)++;

        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_Q_m(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // Q_m

        // Count row
        (*A_row)++;

        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_Q_s(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // Q_s

        // Count row
        (*A_row)++;

        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_v_mag_s(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // vmag_s

        // Count row
        (*A_row)++;

        MAT_set_i(A, *A_nnz, *A_row);
        MAT_set_j(A, *A_nnz, FACTS_get_index_v_ang_s(facts, t));
        MAT_set_d(A, *A_nnz, 1.);
        (*A_nnz)++; // vang_s

        // Count row
        (*A_row)++;

      }

      if (FACTS_is_in_normal_series_mode(facts)) {

        // Get the buses connected to the FACTS device
        bus_k = FACTS_get_bus_k(facts);
        bus_m = FACTS_get_bus_m(facts);

        // Non linear voltage balance equations
        // -vmag_k*cos(vang_k)+vmag_m*cos(vang_m)-vmag_s*cos(vang_s) = 0 (f1)
        // -vmag_k*sin(vang_k)+vmag_m*sin(vang_m)-vmag_s*sin(vang_s) = 0 (f2)
        // Non linear series converter power balance equations
        // vmag_s*P_m*cos(vang_s)-vmag_s*Q_m*sin(vang_s)-vmag_m*P_dc*cos(vang_m)+vmag_m*Q_s*sin(vang_m) = 0 (f3)
        // vmag_s*P_m*sin(vang_s)+vmag_s*Q_m*cos(vang_s)-vmag_m*P_dc*sin(vang_m)-vmag_m*Q_s*cos(vang_m) = 0 (f4)

        // Hessians
        H1 = MAT_array_get(H_array,*J_row);
        H2 = MAT_array_get(H_array,*J_row+1);
        H3 = MAT_array_get(H_array,*J_row+2);
        H4 = MAT_array_get(H_array,*J_row+3);

        MAT_set_i(J, *J_nnz, *J_row+2);
        MAT_set_j(J, *J_nnz, FACTS_get_index_P_m(facts, t));
        (*J_nnz)++; // df3/P_m

        MAT_set_i(J, *J_nnz, *J_row+2);
        MAT_set_j(J, *J_nnz, FACTS_get_index_Q_m(facts, t));
        (*J_nnz)++; // df3/Q_m

        MAT_set_i(J, *J_nnz, *J_row+2);
        MAT_set_j(J, *J_nnz, FACTS_get_index_P_dc(facts, t));
        (*J_nnz)++; // df3/P_dc

        MAT_set_i(J, *J_nnz, *J_row+2);
        MAT_set_j(J, *J_nnz, FACTS_get_index_Q_s(facts, t));
        (*J_nnz)++; // df3/Q_s

        MAT_set_i(J, *J_nnz, *J_row+3);
        MAT_set_j(J, *J_nnz, FACTS_get_index_P_m(facts, t));
        (*J_nnz)++; // df4/P_m

        MAT_set_i(J, *J_nnz, *J_row+3);
        MAT_set_j(J, *J_nnz, FACTS_get_index_Q_m(facts, t));
        (*J_nnz)++; // df4/Q_m

        MAT_set_i(J, *J_nnz, *J_row+3);
        MAT_set_j(J, *J_nnz, FACTS_get_index_P_dc(facts, t));
        (*J_nnz)++; // df4/P_dc

        MAT_set_i(J, *J_nnz, *J_row+3);
        MAT_set_j(J, *J_nnz, FACTS_get_index_Q_s(facts, t));
        (*J_nnz)++; // df4/Q_s

        if (BUS_has_flags(bus_k, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_k, t));
          (*J_nnz)++; // df1/dvmag_k

          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_k, t));
          (*J_nnz)++; // df1/dvang_k

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_k, t));
          (*J_nnz)++; // df2/dvmag_k

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_k, t));
          (*J_nnz)++; // df2/dvang_k

          // H
          MAT_set_i(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_k, t));
          MAT_set_j(H1, H_nnz[*J_row], BUS_get_index_v_mag(bus_k, t));
          H_nnz[*J_row]++; // df1/(dvang_k and dvmag_k)

          MAT_set_i(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_k, t));
          MAT_set_j(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_k, t));
          H_nnz[*J_row]++; // df1/(dvang_k and dvang_k)

          MAT_set_i(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_k, t));
          MAT_set_j(H2, H_nnz[*J_row+1], BUS_get_index_v_mag(bus_k, t));
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvmag_k)

          MAT_set_i(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_k, t));
          MAT_set_j(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_k, t));
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvang_k)

        }

        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_m, t));
          (*J_nnz)++; // df1/dvmag_m

          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_m, t));
          (*J_nnz)++; // df1/dvang_m

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_m, t));
          (*J_nnz)++; // df2/dvmag_m

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_m, t));
          (*J_nnz)++; // df2/dvang_m

          MAT_set_i(J, *J_nnz, *J_row+2);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_m, t));
          (*J_nnz)++; // df3/dvmag_m

          MAT_set_i(J, *J_nnz, *J_row+2);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_m, t));
          (*J_nnz)++; // df3/dvang_m

          MAT_set_i(J, *J_nnz, *J_row+3);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_mag(bus_m, t));
          (*J_nnz)++; // df4/dvmag_m

          MAT_set_i(J, *J_nnz, *J_row+3);
          MAT_set_j(J, *J_nnz, BUS_get_index_v_ang(bus_m, t));
          (*J_nnz)++; // df4/dvang_m

          // H
          MAT_set_i(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H1, H_nnz[*J_row], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row]++; // df1/(dvang_m and dvmag_m)

          MAT_set_i(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H1, H_nnz[*J_row], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row]++; // df1/(dvang_m and dvang_m)

          MAT_set_i(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H2, H_nnz[*J_row+1], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvmag_m)

          MAT_set_i(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H2, H_nnz[*J_row+1], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvang_m)

          MAT_set_i(H3, H_nnz[*J_row+2], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvmag_m)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_P_dc(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvmag_m)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_Q_s(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvmag_m)

          MAT_set_i(H3, H_nnz[*J_row+2], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvang_m)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_P_dc(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvang_m)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_Q_s(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvang_m)

          MAT_set_i(H4, H_nnz[*J_row+3], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvmag_m)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_P_dc(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvmag_m)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_Q_s(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_mag(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvmag_m)

          MAT_set_i(H4, H_nnz[*J_row+3], BUS_get_index_v_ang(bus_m, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvang_m)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_P_dc(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvang_m)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_Q_s(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], BUS_get_index_v_ang(bus_m, t));
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvang_m)

        }

        if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {

          // J
          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_mag_s(facts, t));
          (*J_nnz)++; // df1/dvmag_s

          MAT_set_i(J, *J_nnz, *J_row);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_ang_s(facts, t));
          (*J_nnz)++; // df1/dvang_s

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_mag_s(facts, t));
          (*J_nnz)++; // df2/dvmag_s

          MAT_set_i(J, *J_nnz, *J_row+1);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_ang_s(facts, t));
          (*J_nnz)++; // df2/dvang_s

          MAT_set_i(J, *J_nnz, *J_row+2);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_mag_s(facts, t));
          (*J_nnz)++; // df3/dvmag_s

          MAT_set_i(J, *J_nnz, *J_row+2);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_ang_s(facts, t));
          (*J_nnz)++; // df3/dvang_s

          MAT_set_i(J, *J_nnz, *J_row+3);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_mag_s(facts, t));
          (*J_nnz)++; // df4/dvmag_s

          MAT_set_i(J, *J_nnz, *J_row+3);
          MAT_set_j(J, *J_nnz, FACTS_get_index_v_ang_s(facts, t));
          (*J_nnz)++; // df4/dvang_s

          // H
          MAT_set_i(H1, H_nnz[*J_row], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H1, H_nnz[*J_row], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row]++; // df1/(dvang_s and dvmag_s)

          MAT_set_i(H1, H_nnz[*J_row], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H1, H_nnz[*J_row], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row]++; // df1/(dvang_s and dvang_s)

          MAT_set_i(H2, H_nnz[*J_row+1], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H2, H_nnz[*J_row+1], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvmag_s)

          MAT_set_i(H2, H_nnz[*J_row+1], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H2, H_nnz[*J_row+1], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvang_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvmag_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_P_m(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dP_m and dvmag_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_Q_m(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvmag_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvang_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_P_m(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dP_m and dvang_s)

          MAT_set_i(H3, H_nnz[*J_row+2], FACTS_get_index_Q_m(facts, t));
          MAT_set_j(H3, H_nnz[*J_row+2], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvang_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvmag_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_P_m(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dP_m and dvmag_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_Q_m(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_mag_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvmag_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_v_ang_s(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvang_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_P_m(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dP_m and dvang_s)

          MAT_set_i(H4, H_nnz[*J_row+3], FACTS_get_index_Q_m(facts, t));
          MAT_set_j(H4, H_nnz[*J_row+3], FACTS_get_index_v_ang_s(facts, t));
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvang_s)

        }

        // Count rows
        (*J_row)++; // df1
        (*J_row)++; // df2
        (*J_row)++; // df3
        (*J_row)++; // df4

      }
    }
  }
}

void CONSTR_FACTS_EQ_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Bus* bus_k;
  Bus* bus_m;
  Facts* facts;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* H1;
  REAL* H2;
  REAL* H3;
  REAL* H4;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL v_mag_k;
  REAL v_ang_k;
  REAL v_mag_m;
  REAL v_ang_m;
  REAL v_mag_s;
  REAL v_ang_s;
  REAL P_m;
  REAL Q_m;
  REAL P_dc;
  REAL Q_s;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  J_row = CONSTR_get_J_row_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!f || !J || !J_nnz || !J_row || !H_nnz || !bus)
    return;

  // FACTS
  for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {

    if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_P | FACTS_VAR_Q)) {

      if (FACTS_is_in_normal_series_mode(facts)) {
            
        // Get the buses connected to the FACTS device
        bus_k = FACTS_get_bus_k(facts);
        bus_m = FACTS_get_bus_m(facts);
            
        // Get variable values
        P_m = VEC_get(values, FACTS_get_index_P_m(facts, t));
        Q_m = VEC_get(values, FACTS_get_index_Q_m(facts, t));
        P_dc = VEC_get(values, FACTS_get_index_P_dc(facts, t));
        Q_s = VEC_get(values, FACTS_get_index_Q_s(facts, t));
            
        if (BUS_has_flags(bus_k, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {
          v_mag_k = VEC_get(values, BUS_get_index_v_mag(bus_k, t));
          v_ang_k = VEC_get(values, BUS_get_index_v_ang(bus_k, t));
        }
        else {
          v_mag_k = BUS_get_v_mag(bus_k, t);
          v_ang_k = BUS_get_v_ang(bus_k, t);
        }
            
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {
          v_mag_m = VEC_get(values, BUS_get_index_v_mag(bus_m, t));
          v_ang_m = VEC_get(values, BUS_get_index_v_ang(bus_m, t));
        }
        else {
          v_mag_m = BUS_get_v_mag(bus_m, t);
          v_ang_m = BUS_get_v_ang(bus_m, t);
        }
            
        if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {
          v_mag_s = VEC_get(values, FACTS_get_index_v_mag_s(facts, t));
          v_ang_s = VEC_get(values, FACTS_get_index_v_ang_s(facts, t));
        }
        else {
          v_mag_s = FACTS_get_v_mag_s(facts, t);
          v_ang_s = FACTS_get_v_ang_s(facts, t);
        }

        // Non linear voltage balance equations
        // -v_mag_k*cos(v_ang_k)+v_mag_m*cos(v_ang_m)-v_mag_s*cos(v_ang_s) = 0 (f1)
        // -v_mag_k*sin(v_ang_k)+v_mag_m*sin(v_ang_m)-v_mag_s*sin(v_ang_s) = 0 (f2)
        // Non linear series converter power balance equations
        // v_mag_s*P_m*cos(v_ang_s)-v_mag_s*Q_m*sin(v_ang_s)-v_mag_m*P_dc*cos(v_ang_m)+v_mag_m*Q_s*sin(v_ang_m) = 0 (f3)
        // v_mag_s*P_m*sin(v_ang_s)+v_mag_s*Q_m*cos(v_ang_s)-v_mag_m*P_dc*sin(v_ang_m)-v_mag_m*Q_s*cos(v_ang_m) = 0 (f4)
            
        f[*J_row] = -v_mag_k*cos(v_ang_k)+v_mag_m*cos(v_ang_m)-v_mag_s*cos(v_ang_s);
        f[*J_row+1] = -v_mag_k*sin(v_ang_k)+v_mag_m*sin(v_ang_m)-v_mag_s*sin(v_ang_s);
        f[*J_row+2] = v_mag_s*P_m*cos(v_ang_s)-v_mag_s*Q_m*sin(v_ang_s)-v_mag_m*P_dc*cos(v_ang_m)+v_mag_m*Q_s*sin(v_ang_m);
        f[*J_row+3] = v_mag_s*P_m*sin(v_ang_s)+v_mag_s*Q_m*cos(v_ang_s)-v_mag_m*P_dc*sin(v_ang_m)-v_mag_m*Q_s*cos(v_ang_m);
            
        // Hessian arrays
        H1 = MAT_get_data_array(MAT_array_get(H_array,*J_row));
        H2 = MAT_get_data_array(MAT_array_get(H_array,*J_row+1));
        H3 = MAT_get_data_array(MAT_array_get(H_array,*J_row+2));
        H4 = MAT_get_data_array(MAT_array_get(H_array,*J_row+3));
            
        J[*J_nnz] = v_mag_s*cos(v_ang_s);
        (*J_nnz)++; // df3/P_m

        J[*J_nnz] = -v_mag_s*sin(v_ang_s);
        (*J_nnz)++; // df3/Q_m

        J[*J_nnz] = -v_mag_m*cos(v_ang_m);
        (*J_nnz)++; // df3/P_dc

        J[*J_nnz] = v_mag_m*sin(v_ang_m);
        (*J_nnz)++; // df3/Q_s

        J[*J_nnz] = v_mag_s*sin(v_ang_s);
        (*J_nnz)++; // df4/P_m

        J[*J_nnz] = v_mag_s*cos(v_ang_s);
        (*J_nnz)++; // df4/Q_m

        J[*J_nnz] = -v_mag_m*sin(v_ang_m);
        (*J_nnz)++; // df4/P_dc

        J[*J_nnz] = -v_mag_m*cos(v_ang_m);
        (*J_nnz)++; // df4/Q_s

        if (BUS_has_flags(bus_k, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          J[*J_nnz] = -cos(v_ang_k);
          (*J_nnz)++; // df1/dvmag_k

          J[*J_nnz] = v_mag_k*sin(v_ang_k);
          (*J_nnz)++; // df1/dvang_k

          J[*J_nnz] = -sin(v_ang_k);
          (*J_nnz)++; // df2/dvmag_k

          J[*J_nnz] = -v_mag_k*cos(v_ang_k);
          (*J_nnz)++; // df2/dvang_k

          // H
          H1[H_nnz[*J_row]] = sin(v_ang_k);
          H_nnz[*J_row]++; // df1/(dvang_k and dvmag_k)

          H1[H_nnz[*J_row]] = v_mag_k*cos(v_ang_k);
          H_nnz[*J_row]++; // df1/(dvang_k and dvang_k)

          H2[H_nnz[*J_row+1]] = -cos(v_ang_k);
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvmag_k)

          H2[H_nnz[*J_row+1]] = v_mag_k*sin(v_ang_k);
          H_nnz[*J_row+1]++; // df2/(dvang_k and dvang_k)

        }

        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG | BUS_VAR_VANG)) {

          // J
          J[*J_nnz] = cos(v_ang_m);
          (*J_nnz)++; // df1/dvmag_m

          J[*J_nnz] = -v_mag_m*sin(v_ang_m);
          (*J_nnz)++; // df1/dvang_m

          J[*J_nnz] = sin(v_ang_m);
          (*J_nnz)++; // df2/dvmag_m

          J[*J_nnz] = v_mag_m*cos(v_ang_m);
          (*J_nnz)++; // df2/dvang_m

          J[*J_nnz] = -P_dc*cos(v_ang_m)+Q_s*sin(v_ang_m);
          (*J_nnz)++; // df3/dvmag_m

          J[*J_nnz] = v_mag_m*P_dc*sin(v_ang_m)+v_mag_m*Q_s*cos(v_ang_m);
          (*J_nnz)++; // df3/dvang_m

          J[*J_nnz] = -P_dc*sin(v_ang_m)-Q_s*cos(v_ang_m);
          (*J_nnz)++; // df4/dvmag_m

          J[*J_nnz] = -v_mag_m*P_dc*cos(v_ang_m)+v_mag_m*Q_s*sin(v_ang_m);
          (*J_nnz)++; // df4/dvang_m

          // H
          H1[H_nnz[*J_row]] = -sin(v_ang_m);
          H_nnz[*J_row]++; // df1/(dvang_m and dvmag_m)

          H1[H_nnz[*J_row]] = -v_mag_m*cos(v_ang_m);
          H_nnz[*J_row]++; // df1/(dvang_m and dvang_m)

          H2[H_nnz[*J_row+1]] = cos(v_ang_m);
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvmag_m)

          H2[H_nnz[*J_row+1]] = -v_mag_m*sin(v_ang_m);
          H_nnz[*J_row+1]++; // df2/(dvang_m and dvang_m)

          H3[H_nnz[*J_row+2]] = P_dc*sin(v_ang_m)+Q_s*cos(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvmag_m)

          H3[H_nnz[*J_row+2]] = -cos(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvmag_m)

          H3[H_nnz[*J_row+2]] = sin(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvmag_m)

          H3[H_nnz[*J_row+2]] = v_mag_m*P_dc*cos(v_ang_m)-v_mag_m*Q_s*sin(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dvang_m and dvang_m)

          H3[H_nnz[*J_row+2]] = v_mag_m*sin(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dP_dc and dvang_m)

          H3[H_nnz[*J_row+2]] = v_mag_m*cos(v_ang_m);
          H_nnz[*J_row+2]++; // df3/(dQ_s and dvang_m)

          H4[H_nnz[*J_row+3]] = -P_dc*cos(v_ang_m)+Q_s*sin(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvmag_m)

          H4[H_nnz[*J_row+3]] = -sin(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvmag_m)

          H4[H_nnz[*J_row+3]] = -cos(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvmag_m)

          H4[H_nnz[*J_row+3]] = v_mag_m*P_dc*sin(v_ang_m)+v_mag_m*Q_s*cos(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dvang_m and dvang_m)

          H4[H_nnz[*J_row+3]] = -v_mag_m*cos(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dP_dc and dvang_m)

          H4[H_nnz[*J_row+3]] = v_mag_m*sin(v_ang_m);
          H_nnz[*J_row+3]++; // df4/(dQ_s and dvang_m)

        }

        if (FACTS_has_flags(facts, FLAG_VARS, FACTS_VAR_VMAG_S | FACTS_VAR_VANG_S)) {

          // J
          J[*J_nnz] = -cos(v_ang_s);
          (*J_nnz)++; // df1/dvmag_s

          J[*J_nnz] = v_mag_s*sin(v_ang_s);
          (*J_nnz)++; // df1/dvang_s

          J[*J_nnz] = -sin(v_ang_s);
          (*J_nnz)++; // df2/dvmag_s

          J[*J_nnz] = -v_mag_s*cos(v_ang_s);
          (*J_nnz)++; // df2/dvang_s

          J[*J_nnz] = P_m*cos(v_ang_s)-Q_m*sin(v_ang_s);
          (*J_nnz)++; // df3/dvmag_s

          J[*J_nnz] = -v_mag_s*P_m*sin(v_ang_s)-v_mag_s*Q_m*cos(v_ang_s);
          (*J_nnz)++; // df3/dvang_s

          J[*J_nnz] = P_m*sin(v_ang_s)+Q_m*cos(v_ang_s);
          (*J_nnz)++; // df4/dvmag_s

          J[*J_nnz] = v_mag_s*P_m*cos(v_ang_s)-v_mag_s*Q_m*sin(v_ang_s);
          (*J_nnz)++; // df4/dvang_s

          // H
          H1[H_nnz[*J_row]] = sin(v_ang_s);
          H_nnz[*J_row]++; // df1/(dvang_s and dvmag_s)

          H1[H_nnz[*J_row]] = v_mag_s*cos(v_ang_s);
          H_nnz[*J_row]++; // df1/(dvang_s and dvang_s)

          H2[H_nnz[*J_row+1]] = -cos(v_ang_s);
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvmag_s)

          H2[H_nnz[*J_row+1]] = v_mag_s*sin(v_ang_s);
          H_nnz[*J_row+1]++; // df2/(dvang_s and dvang_s)

          H3[H_nnz[*J_row+2]] = -P_m*sin(v_ang_s)-Q_m*cos(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvmag_s)

          H3[H_nnz[*J_row+2]] = cos(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dP_m and dvmag_s)

          H3[H_nnz[*J_row+2]] = -sin(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvmag_s)

          H3[H_nnz[*J_row+2]] = -v_mag_s*P_m*cos(v_ang_s)+v_mag_s*Q_m*sin(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dvang_s and dvang_s)

          H3[H_nnz[*J_row+2]] = -v_mag_s*sin(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dP_m and dvang_s)

          H3[H_nnz[*J_row+2]] = -v_mag_s*cos(v_ang_s);
          H_nnz[*J_row+2]++; // df3/(dQ_m and dvang_s)

          H4[H_nnz[*J_row+3]] = P_m*cos(v_ang_s)-Q_m*sin(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvmag_s)

          H4[H_nnz[*J_row+3]] = sin(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dP_m and dvmag_s)

          H4[H_nnz[*J_row+3]] = cos(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvmag_s)

          H4[H_nnz[*J_row+3]] = -v_mag_s*P_m*sin(v_ang_s)-v_mag_s*Q_m*cos(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dvang_s and dvang_s)

          H4[H_nnz[*J_row+3]] = v_mag_s*cos(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dP_m and dvang_s)

          H4[H_nnz[*J_row+3]] = -v_mag_s*sin(v_ang_s);
          H_nnz[*J_row+3]++; // df4/(dQ_m and dvang_s)

        }

        // Count rows
        (*J_row)++; // df1
        (*J_row)++; // df2
        (*J_row)++; // df3
        (*J_row)++; // df4
      }
    }
  }
}

void CONSTR_FACTS_EQ_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
