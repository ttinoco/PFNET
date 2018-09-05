/** @file constr_REG_PF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_PF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_PF.h>

Constr* CONSTR_REG_PF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_REG_PF_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_REG_PF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_REG_PF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_REG_PF_store_sens_step);
  CONSTR_set_name(c,"power factor regulation");
  return c;
}

void CONSTR_REG_PF_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
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
  
  // VSC
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {

    // Candidate
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {
          
      // Linear
      //*******
          
      // A
      (*A_nnz)++; // P
      (*A_nnz)++; // Q
      (*A_nnz)++; // y
      (*A_nnz)++; // z
          
      // Count row
      (*A_row)++;
          
      // Nonlinear
      //**********
          
      // J
      (*J_nnz)++; // dCompY/dQ
      (*J_nnz)++; // dCompY/dy
      (*J_nnz)++; // dCompZ/dQ
      (*J_nnz)++; // dCompZ/dz
          
      // H
      H_nnz[*J_row]++;   // y and y (CompY)
      H_nnz[*J_row]++;   // Q and Q (CompY)
      H_nnz[*J_row]++;   // y and Q (CompY)
         
      H_nnz[*J_row+1]++; // z and z (CompZ)
      H_nnz[*J_row+1]++; // Q and Q (CompZ)
      H_nnz[*J_row+1]++; // z and Q (CompZ)
          
      // Count rows
      (*J_row)++; // dCompY
      (*J_row)++; // dCompZ
          
      // Num extra vars
      CONSTR_set_num_extra_vars(c,*J_row);
    }
  }
}

void CONSTR_REG_PF_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
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
  int index_y;
  int index_z;
  int num_vars;
  REAL gamma;
  REAL factor;

  // Number of vars
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

  // Check pointers
  if (!A_nnz || !J_nnz || !A_row || !H_array || !J_row || !H_nnz || !bus)
    return;

  // VSC
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {

    // Candidate
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {

      // PF
      gamma = CONVVSC_get_target_power_factor(vsc);
      factor = sqrt((1.-gamma*gamma)/(gamma*gamma));
          
      // Hessians
      Hy = MAT_array_get(H_array,*J_row);
      Hz = MAT_array_get(H_array,*J_row+1);
          
      // Indices
      index_y = num_vars+(*J_row);
      index_z = num_vars+(*J_row+1);
          
      // Linear
      // Q - factor*P - factor*y + factor*z = 0 (pf >= 0)
      // Q + factor*P - factor*y + factor*z = 0 (pf <  0)
      //*************************************************
          
      // b
      VEC_set(b,*A_row,0);
          
      // A
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,CONVVSC_get_index_P(vsc,t));
      if (gamma >= 0)
        MAT_set_d(A,*A_nnz,-factor);
      else
        MAT_set_d(A,*A_nnz,factor);
      (*A_nnz)++; // P

      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,CONVVSC_get_index_Q(vsc,t));
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++; // Q
          
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_y);
      MAT_set_d(A,*A_nnz,-factor);
      (*A_nnz)++; // y
          
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,index_z);
      MAT_set_d(A,*A_nnz,factor);
      (*A_nnz)++; // z
          
      // Count row
      (*A_row)++;
          
      // Nonlinear
      // (Q-Qmin) + y - sqrt( (Q-Qmin)^2 + y^2 + 2*eta )
      // (Qmax-Q) + z - sqrt( (Qmax-Q)^2 + z^2 + 2*eta )
      //************************************************
          
      // J
      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,CONVVSC_get_index_Q(vsc,t));
      (*J_nnz)++; // dCompY/dQ

      MAT_set_i(J,*J_nnz,*J_row);
      MAT_set_j(J,*J_nnz,index_y);
      (*J_nnz)++; // dCompY/dy
          
      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,CONVVSC_get_index_Q(vsc,t));
      (*J_nnz)++; // dCompZ/dQ
                    
      MAT_set_i(J,*J_nnz,*J_row+1);
      MAT_set_j(J,*J_nnz,index_z);
      (*J_nnz)++; // dCompZ/dz
          
      // H
      MAT_set_i(Hy,H_nnz[*J_row],index_y);
      MAT_set_j(Hy,H_nnz[*J_row],index_y);
      H_nnz[*J_row]++; // y and y (CompY)

      MAT_set_i(Hy,H_nnz[*J_row],CONVVSC_get_index_Q(vsc,t));
      MAT_set_j(Hy,H_nnz[*J_row],CONVVSC_get_index_Q(vsc,t));
      H_nnz[*J_row]++; // Q and Q (CompY)

      MAT_set_i(Hy,H_nnz[*J_row],index_y);
      MAT_set_j(Hy,H_nnz[*J_row],CONVVSC_get_index_Q(vsc,t));
      H_nnz[*J_row]++; // y and Q (CompY)

      MAT_set_i(Hz,H_nnz[*J_row+1],index_z);
      MAT_set_j(Hz,H_nnz[*J_row+1],index_z);
      H_nnz[*J_row+1]++; // z and z (CompZ)
                                          
      MAT_set_i(Hz,H_nnz[*J_row+1],CONVVSC_get_index_Q(vsc,t));
      MAT_set_j(Hz,H_nnz[*J_row+1],CONVVSC_get_index_Q(vsc,t));
      H_nnz[*J_row+1]++; // Q and Q (CompZ)
          
      MAT_set_i(Hz,H_nnz[*J_row+1],index_z);
      MAT_set_j(Hz,H_nnz[*J_row+1],CONVVSC_get_index_Q(vsc,t));
      H_nnz[*J_row+1]++; // z and Q (CompZ)
          
      // Extra var limits
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row,-CONSTR_REG_PF_MAX_YZ);   // y
      VEC_set(CONSTR_get_l_extra_vars(c),*J_row+1,-CONSTR_REG_PF_MAX_YZ); // z
          
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row,CONSTR_REG_PF_MAX_YZ);   // y
      VEC_set(CONSTR_get_u_extra_vars(c),*J_row+1,CONSTR_REG_PF_MAX_YZ); // z
          
      // Count rows
      (*J_row)++; // dCompY
      (*J_row)++; // dCompZ
    }
  }
}

void CONSTR_REG_PF_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {

  // Local variables
  ConvVSC* vsc;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* Hy;
  REAL* Hz;
  int* J_nnz;
  int* J_row;
  int* H_nnz;
  REAL y;
  REAL z;
  REAL Q;
  REAL Qmin;
  REAL Qmax;
  REAL Qy;
  REAL Qz;
  REAL sqrt_termY;
  REAL sqrt_termZ;
  
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
  
  // VSC
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {

    // Candidate
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {
      
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
          
      // Q values
      Q = VEC_get(values,CONVVSC_get_index_Q(vsc,t)); // p.u.
      Qmax = CONVVSC_get_Q_max(vsc); // p.u.
      Qmin = CONVVSC_get_Q_min(vsc); // p.u.
      Qy = Q-Qmin;
      Qz = Qmax-Q;
          
      // Terms
      sqrt_termY = sqrt( Qy*Qy + y*y + 2*CONSTR_REG_PF_PARAM );
      sqrt_termZ = sqrt( Qz*Qz + z*z + 2*CONSTR_REG_PF_PARAM );
          
      // f
      f[*J_row] = Qy + y - sqrt_termY;   // CompY
      f[*J_row+1] = Qz + z - sqrt_termZ; // CompZ
          
      // J
      J[*J_nnz] = 1. - Qy/sqrt_termY;
      (*J_nnz)++; // dcompY/dQ

      J[*J_nnz] = 1. - y/sqrt_termY;
      (*J_nnz)++; // dCompY/dy
          
      J[*J_nnz] = -1. + Qz/sqrt_termZ;
      (*J_nnz)++; // dcompZ/dQ
                    
      J[*J_nnz] = 1. - z/sqrt_termZ;
      (*J_nnz)++; // dCompZ/dz
          
      // H
      Hy[H_nnz[*J_row]] = -(Qy*Qy+2*CONSTR_REG_PF_PARAM)/pow(sqrt_termY,3.);
      H_nnz[*J_row]++; // y and y (CompY)

      Hy[H_nnz[*J_row]] = -(y*y+2*CONSTR_REG_PF_PARAM)/pow(sqrt_termY,3.);
      H_nnz[*J_row]++; // Q and Q (CompY)
          
      Hy[H_nnz[*J_row]] = Qy*y/pow(sqrt_termY,3.);
      H_nnz[*J_row]++; // y and Q (CompZ)
          
      Hz[H_nnz[*J_row+1]] = -(Qz*Qz+2*CONSTR_REG_PF_PARAM)/pow(sqrt_termZ,3.);
      H_nnz[*J_row+1]++; // z and z (CompZ)
                              
      Hz[H_nnz[*J_row+1]] = -(z*z+2*CONSTR_REG_PF_PARAM)/pow(sqrt_termZ,3.);
      H_nnz[*J_row+1]++; // Q and Q (CompZ)
          
      Hz[H_nnz[*J_row+1]] = -Qz*z/pow(sqrt_termZ,3.);
      H_nnz[*J_row+1]++; // z and Q (CompZ)
          
      // Count rows
      (*J_row)++; // dCompY
      (*J_row)++; // dCompZ
    }
  }
}

void CONSTR_REG_PF_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
