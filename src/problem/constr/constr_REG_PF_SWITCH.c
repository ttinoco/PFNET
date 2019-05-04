/** @file constr_REG_PF_SWITCH.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_PF_SWITCH.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/reg_obj.h>
#include <pfnet/constr_REG_PF_SWITCH.h>

struct Constr_REG_PF_SWITCH_Data {
  
  char* fix_flag;
};

Constr* CONSTR_REG_PF_SWITCH_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_REG_PF_SWITCH_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_REG_PF_SWITCH_allocate);
  CONSTR_set_func_analyze_step(c, &CONSTR_REG_PF_SWITCH_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_REG_PF_SWITCH_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_REG_PF_SWITCH_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_REG_PF_SWITCH_free);
  CONSTR_set_name(c,"switching power factor regulation");
  return c;
}

void CONSTR_REG_PF_SWITCH_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
  int* A_nnz;
  int* A_row;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;
      
  // VSC
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {
    
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {
      (*A_nnz) += 2;
      (*A_row) += 1;
    }
  }
}

void CONSTR_REG_PF_SWITCH_allocate(Constr* c) {

  // Local variables
  Constr_REG_PF_SWITCH_Data* data;
  int num_vars;
  ConvVSC* vsc;
  Net* net;
  int i;
  int t;

  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);
  
  // Data (var-dependent)
  CONSTR_REG_PF_SWITCH_free(c);
  data = (Constr_REG_PF_SWITCH_Data*)malloc(sizeof(Constr_REG_PF_SWITCH_Data));
  ARRAY_zalloc(data->fix_flag,char,num_vars);
  for (i = 0; i < NET_get_num_vsc_convs(net,FALSE); i++) {
    vsc = NET_get_vsc_conv(net,i);
    if (CONVVSC_is_in_service(vsc) &&
        CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {
      for (t = 0; t < NET_get_num_periods(net); t++)
        data->fix_flag[CONVVSC_get_index_Q(vsc,t)] = FALSE; // Q starts free
    }
  }
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_REG_PF_SWITCH_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  ConvVSC* vsc;
  int* A_nnz;
  int* A_row;
  Vec* b;
  Mat* A;
  REAL gamma;
  REAL factor;
  REAL Q;
  REAL Q_min;
  REAL Q_max;
  Constr_REG_PF_SWITCH_Data* data;

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  data = (Constr_REG_PF_SWITCH_Data*)CONSTR_get_data(c);

  // Check pointer
  if (!A_nnz || !A_row || !data || !bus)
    return;

  // VSC
  for (vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {
    
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q)) {
      
      gamma = CONVVSC_get_target_power_factor(vsc);
      factor = sqrt((1.-gamma*gamma)/(gamma*gamma));
      
      // Q free
      //*******
      if (!data->fix_flag[CONVVSC_get_index_Q(vsc,t)]) {
        
        // A
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,CONVVSC_get_index_P(vsc,t));
        if (gamma >= 0)
          MAT_set_d(A,*A_nnz,-factor);
        else
          MAT_set_d(A,*A_nnz,factor);
        (*A_nnz)++;
        
        // A
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,CONVVSC_get_index_Q(vsc,t));
        MAT_set_d(A,*A_nnz,1.);
        (*A_nnz)++;
        
        // b
        VEC_set(b,*A_row,0.);
      }
      
      // Q fixed
      //********
      else {
        
        Q = CONVVSC_get_Q(vsc,t);
        Q_max = CONVVSC_get_Q_max(vsc);
        Q_min = CONVVSC_get_Q_min(vsc);
        
        // A
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,CONVVSC_get_index_P(vsc,t));
        MAT_set_d(A,*A_nnz,0.);
        (*A_nnz)++;
        
        // A
        MAT_set_i(A,*A_nnz,*A_row);
        MAT_set_j(A,*A_nnz,CONVVSC_get_index_Q(vsc,t));
        MAT_set_d(A,*A_nnz,1.);
        (*A_nnz)++;
        
        // b
        if (fabs(Q-Q_min) < fabs(Q-Q_max))
          VEC_set(b,*A_row,Q_min);
        else
          VEC_set(b,*A_row,Q_max);
      }
      
      (*A_row)++;
    }
  }
}

void CONSTR_REG_PF_SWITCH_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_REG_PF_SWITCH_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_REG_PF_SWITCH_free(Constr* c) {

  // Local variables
  Constr_REG_PF_SWITCH_Data* data = (Constr_REG_PF_SWITCH_Data*)CONSTR_get_data(c);

  // Free
  if (data) {
    if (data->fix_flag)
      free(data->fix_flag);
    free(data);
  }

  // Clear
  CONSTR_set_data(c,NULL);
}

char* CONSTR_REG_PF_SWITCH_get_flags(Constr* c) {

  // Local variables
  Constr_REG_PF_SWITCH_Data* data = (Constr_REG_PF_SWITCH_Data*)CONSTR_get_data(c);
  
  // Check
  if (!data)
    return NULL;

  // Return
  return data->fix_flag;
}
