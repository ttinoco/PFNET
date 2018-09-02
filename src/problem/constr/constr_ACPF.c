/** @file constr_ACPF.c
 *  @brief This file defines the data structure and routines associated with the constraint of type ACPF.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_ACPF.h>

Constr* CONSTR_ACPF_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_ACPF_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_ACPF_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_ACPF_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_ACPF_store_sens_step);
  CONSTR_set_name(c,"AC power balance");
  return c;
}

void CONSTR_ACPF_count_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Branch* br;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  ConvVSC* vsc_conv;
  int* J_row;
  int* J_nnz;
  int* H_nnz;
  int* HP_nnz;
  int* HQ_nnz;
  BOOL var_v;
  BOOL var_w;

  // Constr data
  J_row = CONSTR_get_J_row_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!J_row || !J_nnz || !H_nnz)
    return;
  
  // Bus data
  var_w = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG);
  var_v = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG);
  HP_nnz = H_nnz+BUS_get_index_P(bus,t);
  HQ_nnz = H_nnz+BUS_get_index_Q(bus,t);

  if (var_w) {
    
    BUS_set_dPdw_index(bus,*J_nnz, t);
    (*J_nnz)++; // dPk/dwk
    
    BUS_set_dQdw_index(bus,*J_nnz, t);
    (*J_nnz)++; // dQk/dwk
    
    BUS_set_dwdw_index(bus,*HP_nnz, t);
    (*HP_nnz)++;
    (*HQ_nnz)++; // dwkdwk
    
    if (var_v) {
      BUS_set_dwdv_index(bus,*HP_nnz, t);
      (*HP_nnz)++;
      (*HQ_nnz)++; // dwkdvk
    }
  }
  
  if (var_v) {
    
    BUS_set_dPdv_index(bus,*J_nnz, t);
    (*J_nnz)++; // dPk/dvk
    
    BUS_set_dQdv_index(bus,*J_nnz, t);
    (*J_nnz)++; // dQk/dvk
    
    BUS_set_dvdv_index(bus,*HP_nnz, t);
    (*HP_nnz)++;
    (*HQ_nnz)++; // dvkdvk
  }
  
  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var
      
      // J
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var
      
      // J
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
      
      // J
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var
      
      // J
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    //*****************************
    if (var_v) { // vk var
      
      // J
      // dPk/dvk
      // dQk/dvk
      
      // H
      // vk an vk
    }
    
    //*****************************
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
      // J
      (*J_nnz)++; // dQk/db
      
      // H
      if (var_v) {
        (*HP_nnz)++;
        (*HQ_nnz)++; // b an vk
      }
    }
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var
      
      // J
      (*J_nnz)++; // dPk/dPl
    }
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var
      
      // J
      (*J_nnz)++; // dQk/dQl
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    
    //*****************************
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) { // Pc and Pd var
      
      // J
      (*J_nnz)++; // Pc
      (*J_nnz)++; // Pd
    }
  }

  // VSC converters
  for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv)) {
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P)) { // P var
      
      // J
      (*J_nnz)++; // dPk/dP
    }
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q)) { // Q var
      
      // J
      (*J_nnz)++; // dQk/dQ
    }
  }
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    HP_nnz = H_nnz+BUS_get_index_P(BRANCH_get_bus_k(br),t);
    HQ_nnz = H_nnz+BUS_get_index_Q(BRANCH_get_bus_k(br),t);
    BRANCH_power_flow_count(br,
                            J_nnz,
                            HP_nnz,
                            t,
                            TRUE,  // Pkm, Qkm
                            TRUE); // ext_idx
    *HQ_nnz = *HP_nnz;
    
    HP_nnz = H_nnz+BUS_get_index_P(BRANCH_get_bus_m(br),t);
    HQ_nnz = H_nnz+BUS_get_index_Q(BRANCH_get_bus_m(br),t);
    BRANCH_power_flow_count(br,
                            J_nnz,
                            HP_nnz,
                            t,
                            FALSE, // Pmk, Qmk
                            TRUE); // ext_idx
    *HQ_nnz = *HP_nnz;
  }
    
  // Rows
  (*J_row)++;
  (*J_row)++;  
}


void CONSTR_ACPF_analyze_step(Constr* c, Bus* bus, int t) {

  // Local variables
  Branch* br;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  ConvVSC* vsc_conv;
  Mat* J;
  int* J_row;
  int* J_nnz;
  int* H_nnz;
  Mat* H_array;
  Mat* HP;
  Mat* HQ;
  int v_index;
  int w_index;
  int P_index;
  int Q_index;
  int* HP_nnz;
  int* HQ_nnz;
  BOOL var_v;
  BOOL var_w;
  
  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_row = CONSTR_get_J_row_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!J_row || !J_nnz || !H_nnz || !H_array)
    return;

  // Bus data
  w_index = BUS_get_index_v_ang(bus,t);
  v_index = BUS_get_index_v_mag(bus,t);
  var_w = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG);
  var_v = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG);
  P_index = BUS_get_index_P(bus,t);
  Q_index = BUS_get_index_Q(bus,t);
  HP_nnz = H_nnz+P_index;
  HQ_nnz = H_nnz+Q_index;
  HP = MAT_array_get(H_array,P_index);
  HQ = MAT_array_get(H_array,Q_index);

  if (var_w) {	
    
    MAT_set_i(J,BUS_get_dPdw_index(bus,t),P_index);
    MAT_set_j(J,BUS_get_dPdw_index(bus,t),w_index);
    (*J_nnz)++; // dPk/dwk
    
    MAT_set_i(J,BUS_get_dQdw_index(bus,t),Q_index);
    MAT_set_j(J,BUS_get_dQdw_index(bus,t),w_index);
    (*J_nnz)++; // dQk/dwk
    
    MAT_set_i(HP,BUS_get_dwdw_index(bus,t),w_index);
    MAT_set_j(HP,BUS_get_dwdw_index(bus,t),w_index);
    (*HP_nnz)++;
    MAT_set_i(HQ,BUS_get_dwdw_index(bus,t),w_index);
    MAT_set_j(HQ,BUS_get_dwdw_index(bus,t),w_index);
    (*HQ_nnz)++; // dwkdwk
    
    if (var_v) {
      
      MAT_set_i(HP,BUS_get_dwdv_index(bus,t),w_index);
      MAT_set_j(HP,BUS_get_dwdv_index(bus,t),v_index);
      (*HP_nnz)++;
      MAT_set_i(HQ,BUS_get_dwdv_index(bus,t),w_index);
      MAT_set_j(HQ,BUS_get_dwdv_index(bus,t),v_index);
      (*HQ_nnz)++; // dwkdvk
    }
  }
  
  if (var_v) {
    
    MAT_set_i(J,BUS_get_dPdv_index(bus,t),P_index);
    MAT_set_j(J,BUS_get_dPdv_index(bus,t),v_index);
    (*J_nnz)++; // dPk/dvk
    
    MAT_set_i(J,BUS_get_dQdv_index(bus,t),Q_index);
    MAT_set_j(J,BUS_get_dQdv_index(bus,t),v_index);
    (*J_nnz)++; // dQk/dvk
    
    MAT_set_i(HP,BUS_get_dvdv_index(bus,t),v_index);
    MAT_set_j(HP,BUS_get_dvdv_index(bus,t),v_index);
    (*HP_nnz)++;
    MAT_set_i(HQ,BUS_get_dvdv_index(bus,t),v_index);
    MAT_set_j(HQ,BUS_get_dvdv_index(bus,t),v_index);
    (*HQ_nnz)++; // dvkdvk
  }
  
  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var
      
      // J
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,GEN_get_index_P(gen,t));
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var
      
      // J
      MAT_set_i(J,*J_nnz,Q_index);
      MAT_set_j(J,*J_nnz,GEN_get_index_Q(gen,t));
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
      
      // J
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,VARGEN_get_index_P(vargen,t));
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var
      
      // J
      MAT_set_i(J,*J_nnz,Q_index);
      MAT_set_j(J,*J_nnz,VARGEN_get_index_Q(vargen,t));
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    //*****************************
    if (var_v) { // vk var
      
      // J
      // dPk/dvk
      // dQk/dvk
      
      // H
      // vk an vk
    }
    
    //**************************************
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
      // J
      MAT_set_i(J,*J_nnz,Q_index);
      MAT_set_j(J,*J_nnz,SHUNT_get_index_b(shunt,t));
      (*J_nnz)++; // dQk/db
      
      // H
      if (var_v) {
	
        MAT_set_i(HP,*HP_nnz,SHUNT_get_index_b(shunt,t));
        MAT_set_j(HP,*HP_nnz,v_index);
        (*HP_nnz)++;
        MAT_set_i(HQ,*HQ_nnz,SHUNT_get_index_b(shunt,t));
        MAT_set_j(HQ,*HQ_nnz,v_index);
        (*HQ_nnz)++; // b and vk
      }
    }
  }
      
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var
      
      // J
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,LOAD_get_index_P(load,t));
      (*J_nnz)++; // dPk/dPl
    }
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var
      
      // J
      MAT_set_i(J,*J_nnz,Q_index);
      MAT_set_j(J,*J_nnz,LOAD_get_index_Q(load,t));
      (*J_nnz)++; // dQk/dQl
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    
    //*****************************
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {  // Pc and Pd var
      
      // J
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,BAT_get_index_Pc(bat,t));
      (*J_nnz)++; // Pc
      
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,BAT_get_index_Pd(bat,t));
      (*J_nnz)++; // Pd
    }
  }

  // VSC converters
  for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv)) {
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P)) { // P var
      
      // J
      MAT_set_i(J,*J_nnz,P_index);
      MAT_set_j(J,*J_nnz,CONVVSC_get_index_P(vsc_conv,t));
      (*J_nnz)++; // dPk/dP
    }
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q)) { // Q var
      
      // J
      MAT_set_i(J,*J_nnz,Q_index);
      MAT_set_j(J,*J_nnz,CONVVSC_get_index_Q(vsc_conv,t));
      (*J_nnz)++; // dQk/dQ
    }
  }
  
  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    P_index = BUS_get_index_P(BRANCH_get_bus_k(br),t);
    Q_index = BUS_get_index_Q(BRANCH_get_bus_k(br),t);
    BRANCH_power_flow_analyze(br,
                              J_nnz,
                              J,
                              P_index,
                              Q_index,
                              H_nnz+P_index,
                              MAT_array_get(H_array,P_index),
                              MAT_array_get(H_array,Q_index),
                              t,
                              TRUE,  // Pkm, Qkm
                              TRUE); // ext_idx
    *(H_nnz+Q_index) = *(H_nnz+P_index);

    P_index = BUS_get_index_P(BRANCH_get_bus_m(br),t);
    Q_index = BUS_get_index_Q(BRANCH_get_bus_m(br),t);
    BRANCH_power_flow_analyze(br,
                              J_nnz,
                              J,
                              P_index,
                              Q_index,
                              H_nnz+P_index,
                              MAT_array_get(H_array,P_index),
                              MAT_array_get(H_array,Q_index),
                              t,
                              FALSE, // Pmk, Qmk
                              TRUE); // ext_idx
    *(H_nnz+Q_index) = *(H_nnz+P_index);
  }
  
  // Rows
  (*J_row)++;
  (*J_row)++;
}

void CONSTR_ACPF_eval_step(Constr* c, Bus* bus, int t, Vec* values, Vec* values_extra) {

  // Local variables
  Branch* br;
  Gen* gen;
  Vargen* vargen;
  Load* load;
  Bat* bat;
  Shunt* shunt;
  ConvVSC* vsc_conv;
  ConvCSC* csc_conv;
  REAL* f;
  REAL* J;
  int* J_row;
  int* J_nnz;
  int* H_nnz;
  Mat* H_array;
  REAL* HP;
  REAL* HQ;

  BOOL var_v;
  BOOL var_w;
  int P_index;
  int Q_index;
  int* HP_nnz;
  int* HQ_nnz;

  REAL v;
  
  REAL P;
  REAL Q;

  REAL shunt_b;
  REAL shunt_g;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_row = CONSTR_get_J_row_ptr(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  H_nnz = CONSTR_get_H_nnz(c);

  // Check pointers
  if (!f || !J || !J_row || !J_nnz || !H_nnz)
    return;

  // Bus data
  var_w = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG);
  var_v = BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG);
  P_index = BUS_get_index_P(bus,t); // index in f for active power mismatch
  Q_index = BUS_get_index_Q(bus,t); // index in f for reactive power mismatch
  HP = MAT_get_data_array(MAT_array_get(H_array,P_index));
  HQ = MAT_get_data_array(MAT_array_get(H_array,Q_index));
  HP_nnz = H_nnz+P_index;
  HQ_nnz = H_nnz+Q_index;
  if (var_v)
    v = VEC_get(values,BUS_get_index_v_mag(bus,t));
  else
    v = BUS_get_v_mag(bus,t);

  if (var_w) {	
    (*J_nnz)++; // dPk/dwk
    (*J_nnz)++; // dQk/dwk
    (*HP_nnz)++;
    (*HQ_nnz)++; // dwkdwk
    if (var_v) {
      (*HP_nnz)++;
      (*HQ_nnz)++; // dwkdvk
    }
  }
  if (var_v) {
    (*J_nnz)++; // dPk/dvk
    (*J_nnz)++; // dQk/dvk
    (*HP_nnz)++;
    (*HQ_nnz)++; // dvkdvk
  }
  
  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Outage
    if (GEN_is_on_outage(gen))
      continue;
    
    // Var values
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
      P = VEC_get(values,GEN_get_index_P(gen,t)); // p.u.
    else
      P = GEN_get_P(gen,t);                       // p.u.
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
      Q = VEC_get(values,GEN_get_index_Q(gen,t)); // p.u.
    else
      Q = GEN_get_Q(gen,t);                       // p.u.
    
    // f
    f[P_index] += P;
    f[Q_index] += Q;
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // Pg var
      
      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Qg var
      
      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    
    // Var values
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
      P = VEC_get(values,VARGEN_get_index_P(vargen,t)); // p.u.
    else
      P = VARGEN_get_P(vargen,t);                       // p.u.
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q))
      Q = VEC_get(values,VARGEN_get_index_Q(vargen,t)); // p.u.
    else
      Q = VARGEN_get_Q(vargen,t);                       // p.u.
    
    // f
    f[P_index] += P;
    f[Q_index] += Q;
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) { // Pg var
      
      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dPk/dPg
    }
    
    //*****************************
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Qg var
      
      // J
      J[*J_nnz] = 1.;
      (*J_nnz)++; // dQk/dQg
    }
  }
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    // Values
    shunt_g = SHUNT_get_g(shunt);
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC))
      shunt_b = VEC_get(values,SHUNT_get_index_b(shunt,t)); // p.u.
    else
      shunt_b = SHUNT_get_b(shunt,t);
    
    // f
    f[P_index] -= shunt_g*v*v; // p.u.
    f[Q_index] += shunt_b*v*v;  // p.u.
    
    //*****************************
    if (var_v) { // var v
      
      // J
      J[BUS_get_dPdv_index(bus,t)] += -2*shunt_g*v;
      // dPk/dvk
      J[BUS_get_dQdv_index(bus,t)] += 2*shunt_b*v;
      // dQk/dvk
      
      // H
      HP[BUS_get_dvdv_index(bus,t)] += -2*shunt_g;
      HQ[BUS_get_dvdv_index(bus,t)] += 2*shunt_b;
      // vk and vk
    }
    
    //*****************************
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
      // J
      J[*J_nnz] = v*v;
      (*J_nnz)++; // dQk/db
      
      // H
      if (var_v) {
        HP[*HP_nnz] = 0;
        HQ[*HQ_nnz] = 2*v;
        (*HP_nnz)++;
        (*HQ_nnz)++; // b and vk
      }
    }
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    // Var values
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
      P = VEC_get(values,LOAD_get_index_P(load,t)); // p.u.
    else
      P = LOAD_get_P(load,t);                       // p.u.
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q))
      Q = VEC_get(values,LOAD_get_index_Q(load,t)); // p.u.
    else
      Q = LOAD_get_Q(load,t);                       // p.u.
    
    // f
    f[P_index] -= P; // p.u.
    f[Q_index] -= Q; // p.u.
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) { // Pl var
      
      // J
      J[*J_nnz] = -1.; // dPk/dPl
      (*J_nnz)++;
    }
    
    //*****************************
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) { // Ql var
      
      // J
      J[*J_nnz] = -1.; // dQk/dQl
      (*J_nnz)++;
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    
    // var values
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P))
      P = VEC_get(values,BAT_get_index_Pc(bat,t))-VEC_get(values,BAT_get_index_Pd(bat,t)); // p.u.
    else
      P = BAT_get_P(bat,t);                                                                // p.u.
    
    // f
    f[P_index] -= P; // p.u.
    
    //*****************************
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {  // Pc and Pd var
      
      // J
      J[*J_nnz] = -1.; // Pc
      (*J_nnz)++;
      
      J[*J_nnz] = 1.; // Pd
      (*J_nnz)++;
    }
  }

  // CSC converters
  for (csc_conv = BUS_get_csc_conv(bus); csc_conv != NULL; csc_conv = CONVCSC_get_next_ac(csc_conv)) {
    
    // f
    f[P_index] += CONVCSC_get_P(csc_conv,t); // p.u.
    f[Q_index] += CONVCSC_get_Q(csc_conv,t); // p.u.
  }
  
  // VSC converters
  for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv)) {
    
    // Var values
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P))
      P = VEC_get(values,CONVVSC_get_index_P(vsc_conv,t)); // p.u.
    else
      P = CONVVSC_get_P(vsc_conv,t);                       // p.u.
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q))
      Q = VEC_get(values,CONVVSC_get_index_Q(vsc_conv,t)); // p.u.
    else
      Q = CONVVSC_get_Q(vsc_conv,t);                       // p.u.
    
    // f
    f[P_index] += P; // p.u.
    f[Q_index] += Q; // p.u.
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P)) { // P var
      
      // J
      J[*J_nnz] = 1.; // dPk/dP
      (*J_nnz)++;
    }
    
    //*****************************
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q)) { // Q var
      
      // J
      J[*J_nnz] = 1.; // dQk/dQ
      (*J_nnz)++;
    }
  }

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {
    
    P_index = BUS_get_index_P(BRANCH_get_bus_k(br),t);
    Q_index = BUS_get_index_Q(BRANCH_get_bus_k(br),t);
    BRANCH_power_flow_eval(br,
                           f+P_index,
                           f+Q_index,
                           J_nnz,
                           J,
                           H_nnz+P_index,
                           MAT_get_data_array(MAT_array_get(H_array,P_index)),
                           MAT_get_data_array(MAT_array_get(H_array,Q_index)),
                           values,
                           -1.,   // flows leaving bus are negative
                           t,
                           TRUE,  // Pkm, Qkm
                           TRUE); // ext_idx
    *(H_nnz+Q_index) = *(H_nnz+P_index);

    P_index = BUS_get_index_P(BRANCH_get_bus_m(br),t);
    Q_index = BUS_get_index_Q(BRANCH_get_bus_m(br),t);
    BRANCH_power_flow_eval(br,
                           f+P_index,
                           f+Q_index,
                           J_nnz,
                           J,
                           H_nnz+P_index,
                           MAT_get_data_array(MAT_array_get(H_array,P_index)),
                           MAT_get_data_array(MAT_array_get(H_array,Q_index)),
                           values,
                           -1.,    // flows leaving bus are negative
                           t,
                           FALSE,  // Pmk, Qmk
                           TRUE);  // ext_idx
    *(H_nnz+Q_index) = *(H_nnz+P_index);
  }
  
  // Rows
  (*J_row)++;
  (*J_row)++;
}

void CONSTR_ACPF_store_sens_step(Constr* c, Bus* bus, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
  BUS_set_sens_P_balance(bus,VEC_get(sf,BUS_get_index_P(bus,t)),t); // sens of P balance
  BUS_set_sens_Q_balance(bus,VEC_get(sf,BUS_get_index_Q(bus,t)),t); // sens of Q balance
}
