/** @file constr_BOUND.c
 *  @brief This file defines the data structure and routines associated with the constraint of type BOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_BOUND.h>

Constr* CONSTR_BOUND_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_BOUND_count_step);
  CONSTR_set_func_analyze_step(c,&CONSTR_BOUND_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_BOUND_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_BOUND_store_sens_step);
  CONSTR_set_name(c,"variable bounds");
  return c;
}

void CONSTR_BOUND_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  CONSTR_set_G_row(c,NET_get_num_vars(CONSTR_get_network(c)));
  CONSTR_set_G_nnz(c,NET_get_num_vars(CONSTR_get_network(c)));
}

void CONSTR_BOUND_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Branch* br;
  Gen* gen;
  Load* load;
  Vargen* vargen;
  Shunt* shunt;
  Bat* bat;
  ConvVSC* vsc_conv;
  ConvCSC* csc_conv;
  Facts* facts;
  Vec* l;
  Vec* u;
  Mat* G;
  int index;
  int index1;
  int index2;
  int index3;
  int index4;

  // Cosntr data
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G = CONSTR_get_G(c);

  // Check pointer
  if (!G || !u || !l)
    return;
  
  // Voltage magnitude (V_MAG)
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
    index = BUS_get_index_v_mag(bus,t);
    MAT_set_i(G,index,index);
    MAT_set_j(G,index,index);
    MAT_set_d(G,index,1.);
    if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG)) {
      VEC_set(u,index,BUS_get_v_max_norm(bus));
      VEC_set(l,index,BUS_get_v_min_norm(bus));
    }
    else {
      VEC_set(u,index,BUS_INF_V_MAG);
      VEC_set(l,index,0.);
    }
    
    // Row info
    CONSTR_set_G_row_info_string(c,
                                 index,
                                 "bus",
                                 BUS_get_index(bus),
                                 "voltage magnitude",
                                 t);
  }
  
  // Voltage angle (V_ANG)
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
    index = BUS_get_index_v_ang(bus,t);
    MAT_set_i(G,index,index);
    MAT_set_j(G,index,index);
    MAT_set_d(G,index,1.);
    VEC_set(u,index,BUS_INF_V_ANG);
    VEC_set(l,index,-BUS_INF_V_ANG);
    
    // Row info
    CONSTR_set_G_row_info_string(c,
                                 index,
                                 "bus",
                                 BUS_get_index(bus),
                                 "voltage angle",
                                 t);
  }
  
  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Active power (P)
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
      index = GEN_get_index_P(gen,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P)) {
        VEC_set(u,index,GEN_get_P_max(gen));
        VEC_set(l,index,GEN_get_P_min(gen));
      }
      else {
        VEC_set(u,index,GEN_INF_P);
        VEC_set(l,index,-GEN_INF_P);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "generator",
                                   GEN_get_index(gen),
                                   "active power",
                                   t);
    }
    
    // Reactive power (Q)
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
      index = GEN_get_index_Q(gen,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q)) {
        VEC_set(u,index,GEN_get_Q_max(gen));
        VEC_set(l,index,GEN_get_Q_min(gen));
      }
      else {
        VEC_set(u,index,GEN_INF_Q);
        VEC_set(l,index,-GEN_INF_Q);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "generator",
                                   GEN_get_index(gen),
                                   "reactive power",
                                   t);
    }
  }

  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    // Active power (P)
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
      index = LOAD_get_index_P(load,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (LOAD_has_flags(load,FLAG_BOUNDED,LOAD_VAR_P)) {
        VEC_set(u,index,LOAD_get_P_max(load,t));
        VEC_set(l,index,LOAD_get_P_min(load,t));
      }
      else {
        VEC_set(u,index,LOAD_INF_P);
        VEC_set(l,index,-LOAD_INF_P);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "load",
                                   LOAD_get_index(load),
                                   "active power",
                                   t);
    }
    
    // Rective power (Q)
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {
      index = LOAD_get_index_Q(load,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (LOAD_has_flags(load,FLAG_BOUNDED,LOAD_VAR_Q)) {
        VEC_set(u,index,LOAD_get_Q_max(load,t));
        VEC_set(l,index,LOAD_get_Q_min(load,t));
      }
      else {
        VEC_set(u,index,LOAD_INF_Q);
        VEC_set(l,index,-LOAD_INF_Q);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "load",
                                   LOAD_get_index(load),
                                   "reactive power",
                                   t);
    }
  }
  
  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
    
    // Active power (P)
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
      index = VARGEN_get_index_P(vargen,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (VARGEN_has_flags(vargen,FLAG_BOUNDED,VARGEN_VAR_P)) {
        VEC_set(u,index,VARGEN_get_P_ava(vargen,t));
        VEC_set(l,index,VARGEN_get_P_min(vargen));
      }
      else {
        VEC_set(u,index,VARGEN_INF_P);
        VEC_set(l,index,-VARGEN_INF_P);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "variable generator",
                                   VARGEN_get_index(vargen),
                                   "active power",
                                   t);
    }
    
    // Reactive power (Q)
    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
      index = VARGEN_get_index_Q(vargen,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (VARGEN_has_flags(vargen,FLAG_BOUNDED,VARGEN_VAR_Q)) {
        VEC_set(u,index,VARGEN_get_Q_max(vargen));
        VEC_set(l,index,VARGEN_get_Q_min(vargen));
      }
      else {
        VEC_set(u,index,VARGEN_INF_Q);
        VEC_set(l,index,-VARGEN_INF_Q);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "variable generator",
                                   VARGEN_get_index(vargen),
                                   "reactive power",
                                   t);
    }
  }
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    // Susceptance (b)
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
      index = SHUNT_get_index_b(shunt,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC)) {
        VEC_set(u,index,SHUNT_get_b_max(shunt));
        VEC_set(l,index,SHUNT_get_b_min(shunt));
      }
      else {
        VEC_set(u,index,SHUNT_INF_SUSC);
        VEC_set(l,index,-SHUNT_INF_SUSC);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "shunt",
                                   SHUNT_get_index(shunt),
                                   "susceptance",
                                   t);
    }
  }
  
  // Batteries
  for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {
    
    // Charging/discharging power (P)
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
      
      index1 = BAT_get_index_Pc(bat,t);
      index2 = BAT_get_index_Pd(bat,t);
      
      MAT_set_i(G,index1,index1);
      MAT_set_j(G,index1,index1);
      MAT_set_d(G,index1,1.);
      
      MAT_set_i(G,index2,index2);
      MAT_set_j(G,index2,index2);
      MAT_set_d(G,index2,1.);
      
      if (BAT_has_flags(bat,FLAG_BOUNDED,BAT_VAR_P)) {
	
        VEC_set(u,index1,BAT_get_P_max(bat));
        VEC_set(l,index1,0.);
	
        VEC_set(u,index2,-BAT_get_P_min(bat));
        VEC_set(l,index2,0.);
      }
      else {
	
        VEC_set(u,index1,BAT_INF_P);
        VEC_set(l,index1,0.);
	
        VEC_set(u,index2,BAT_INF_P);
        VEC_set(l,index2,0.);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index1,
                                   "battery",
                                   BAT_get_index(bat),
                                   "charging power",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index2,
                                   "battery",
                                   BAT_get_index(bat),
                                   "discharging power",
                                   t);
    }
    
    // Energy level (E)
    if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E)) {
      index = BAT_get_index_E(bat,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (BAT_has_flags(bat,FLAG_BOUNDED,BAT_VAR_E)) {
        VEC_set(u,index,BAT_get_E_max(bat));
        VEC_set(l,index,0.);
      }
      else {
        VEC_set(u,index,BAT_INF_E);
        VEC_set(l,index,0.);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "battery",
                                   BAT_get_index(bat),
                                   "energy level",
                                   t);
    }
  }

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {
    
    // Tap ratio
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
      index = BRANCH_get_index_ratio(br,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO)) {
        VEC_set(u,index,BRANCH_get_ratio_max(br));
        VEC_set(l,index,BRANCH_get_ratio_min(br));
      }
      else {
        VEC_set(u,index,BRANCH_INF_RATIO);
        VEC_set(l,index,0.);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "branch",
                                   BRANCH_get_index(br),
                                   "tap ratio",
                                   t);
    }
    
    // Phase shift
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
      index = BRANCH_get_index_phase(br,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE)) {
        VEC_set(u,index,BRANCH_get_phase_max(br));
        VEC_set(l,index,BRANCH_get_phase_min(br));
      }
      else {
        VEC_set(u,index,2*PI);
        VEC_set(l,index,-2*PI);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "branch",
                                   BRANCH_get_index(br),
                                   "phase shift",
                                   t);
    }
  }

  // DC voltages (v)
  if (BUSDC_has_flags(busdc,FLAG_VARS,BUSDC_VAR_V)) {
    index = BUSDC_get_index_v(busdc,t);
    MAT_set_i(G,index,index);
    MAT_set_j(G,index,index);
    MAT_set_d(G,index,1.);
    VEC_set(u,index,BUSDC_INF_V);
    VEC_set(l,index,-BUSDC_INF_V);

    // Row info
    CONSTR_set_G_row_info_string(c,
                                 index,
                                 "dc bus",
                                 BUSDC_get_index(busdc),
                                 "voltage",
                                 t);
  }
  
  // VSC Converters
  for (vsc_conv = BUS_get_vsc_conv(bus); vsc_conv != NULL; vsc_conv = CONVVSC_get_next_ac(vsc_conv)) {
    
    // Active Power (P)
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_P)) {
      index = CONVVSC_get_index_P(vsc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (CONVVSC_has_flags(vsc_conv,FLAG_BOUNDED,CONVVSC_VAR_P)) {
      	VEC_set(u,index,CONVVSC_get_P_max(vsc_conv));
      	VEC_set(l,index,CONVVSC_get_P_min(vsc_conv));
      }
      else {
      	VEC_set(u,index,CONVVSC_INF_P);
      	VEC_set(l,index,-CONVVSC_INF_P);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "vsc converter",
                                   CONVVSC_get_index(vsc_conv),
                                   "active power",
                                   t);
    }
    
    // Reactive Power (Q)
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_Q)) {
      index = CONVVSC_get_index_Q(vsc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (CONVVSC_has_flags(vsc_conv,FLAG_BOUNDED,CONVVSC_VAR_Q)) {
      	VEC_set(u,index,CONVVSC_get_Q_max(vsc_conv));
      	VEC_set(l,index,CONVVSC_get_Q_min(vsc_conv));
      }
      else {
      	VEC_set(u,index,CONVVSC_INF_Q);
      	VEC_set(l,index,-CONVVSC_INF_Q);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "vsc converter",
                                   CONVVSC_get_index(vsc_conv),
                                   "reactive power",
                                   t);
    }
    
    // DC Power (P_dc)
    if (CONVVSC_has_flags(vsc_conv,FLAG_VARS,CONVVSC_VAR_PDC)) {
      index1 = CONVVSC_get_index_P_dc(vsc_conv,t);
      index2 = CONVVSC_get_index_i_dc(vsc_conv,t);
      
      MAT_set_i(G,index1,index1);
      MAT_set_j(G,index1,index1);
      MAT_set_d(G,index1,1.);
      
      MAT_set_i(G,index2,index2);
      MAT_set_j(G,index2,index2);
      MAT_set_d(G,index2,1.);
      
      VEC_set(u,index1,CONVVSC_INF_PDC);
      VEC_set(l,index1,-CONVVSC_INF_PDC);
      
      VEC_set(u,index2,CONVVSC_INF_PDC);
      VEC_set(l,index2,-CONVVSC_INF_PDC);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index1,
                                   "vsc converter",
                                   CONVVSC_get_index(vsc_conv),
                                   "dc power",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index2,
                                   "vsc converter",
                                   CONVVSC_get_index(vsc_conv),
                                   "dc current",
                                   t);
    }
  }

  // CSC Converters
  for (csc_conv = BUS_get_csc_conv(bus); csc_conv != NULL; csc_conv = CONVCSC_get_next_ac(csc_conv)) {
    
    // Active Power (P)
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_P)) {
      index = CONVCSC_get_index_P(csc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      VEC_set(u,index,CONVCSC_INF_P);
      VEC_set(l,index,-CONVCSC_INF_P);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "active power",
                                   t);
    }
    
    // Reactive Power (Q)
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_Q)) {
      index = CONVCSC_get_index_Q(csc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      VEC_set(u,index,CONVCSC_INF_Q);
      VEC_set(l,index,-CONVCSC_INF_Q);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "reactive power",
                                   t);
    }
    
    // DC Power (P_dc)
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_PDC)) {
      index1 = CONVCSC_get_index_P_dc(csc_conv,t);
      index2 = CONVCSC_get_index_i_dc(csc_conv,t);
      
      MAT_set_i(G,index1,index1);
      MAT_set_j(G,index1,index1);
      MAT_set_d(G,index1,1.);
      
      MAT_set_i(G,index2,index2);
      MAT_set_j(G,index2,index2);
      MAT_set_d(G,index2,1.);
      
      VEC_set(u,index1,CONVCSC_INF_PDC);
      VEC_set(l,index1,-CONVCSC_INF_PDC);
      
      VEC_set(u,index2,CONVCSC_INF_PDC);
      VEC_set(l,index2,-CONVCSC_INF_PDC);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index1,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "dc power",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index2,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "dc current",
                                   t);
    }

    // Angle
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_ANGLE)) {
      index = CONVCSC_get_index_angle(csc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      VEC_set(u,index,CONVCSC_INF_ANGLE);
      VEC_set(l,index,-CONVCSC_INF_ANGLE);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "angle",
                                   t);
    }

    // Ratio
    if (CONVCSC_has_flags(csc_conv,FLAG_VARS,CONVCSC_VAR_RATIO)) {
      index = CONVCSC_get_index_ratio(csc_conv,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      VEC_set(u,index,CONVCSC_INF_RATIO);
      VEC_set(l,index,-CONVCSC_INF_RATIO);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "csc converter",
                                   CONVCSC_get_index(csc_conv),
                                   "tap ratio",
                                   t);
    }
  }
  
  // FACTS
  for (facts = BUS_get_facts_k(bus); facts != NULL; facts = FACTS_get_next_k(facts)) {
    
    // Series Voltage Magnitude (v_mag_s)
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_VMAG_S)) {
      index = FACTS_get_index_v_mag_s(facts,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);
      if (FACTS_has_flags(facts,FLAG_BOUNDED,FACTS_VAR_VMAG_S)) {
      	VEC_set(u,index,FACTS_get_v_max_s(facts));
      	VEC_set(l,index,0.);
      }
      else {
      	VEC_set(u,index,FACTS_INF_VMAG_S);
      	VEC_set(l,index,0.);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "series voltage magnitude",
                                   t);
    }
    
    // Series Voltage Angle (v_ang_s)
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_VANG_S)) {
      index = FACTS_get_index_v_ang_s(facts,t);
      MAT_set_i(G,index,index);
      MAT_set_j(G,index,index);
      MAT_set_d(G,index,1.);

      VEC_set(u,index,FACTS_INF_VANG_S);
      VEC_set(l,index,-FACTS_INF_VANG_S);
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "series voltage angle",
                                   t);
    }
    
    // Active Power (P)
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_P)) {
      index1 = FACTS_get_index_P_k(facts,t);
      index2 = FACTS_get_index_P_m(facts,t);
      index3 = FACTS_get_index_P_dc(facts,t);
      
      MAT_set_i(G,index1,index1);
      MAT_set_j(G,index1,index1);
      MAT_set_d(G,index1,1.);

      MAT_set_i(G,index2,index2);
      MAT_set_j(G,index2,index2);
      MAT_set_d(G,index2,1.);

      MAT_set_i(G,index3,index3);
      MAT_set_j(G,index3,index3);
      MAT_set_d(G,index3,1.);
      
      VEC_set(u,index1,FACTS_INF_P);
      VEC_set(l,index1,-FACTS_INF_P);
      
      VEC_set(u,index2,FACTS_INF_P);
      VEC_set(l,index2,-FACTS_INF_P);

      if (FACTS_has_flags(facts,FLAG_BOUNDED,FACTS_VAR_P)) {
      	VEC_set(u,index3,FACTS_get_P_max_dc(facts));
      	VEC_set(l,index3,-FACTS_get_P_max_dc(facts));
      }
      else {
      	VEC_set(u,index3,FACTS_INF_P);
      	VEC_set(l,index3,-FACTS_INF_P);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index1,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "active power k",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index2,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "active power m",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index3,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "dc power",
                                   t);
    }
    
    // Reactive Powers
    if (FACTS_has_flags(facts,FLAG_VARS,FACTS_VAR_Q)) {
      index1 = FACTS_get_index_Q_k(facts,t);
      index2 = FACTS_get_index_Q_m(facts,t);
      index3 = FACTS_get_index_Q_s(facts,t);
      index4 = FACTS_get_index_Q_sh(facts,t);
      
      MAT_set_i(G,index1,index1);
      MAT_set_j(G,index1,index1);
      MAT_set_d(G,index1,1.);
      
      MAT_set_i(G,index2,index2);
      MAT_set_j(G,index2,index2);
      MAT_set_d(G,index2,1.);
      
      MAT_set_i(G,index3,index3);
      MAT_set_j(G,index3,index3);
      MAT_set_d(G,index3,1.);

      MAT_set_i(G,index4,index4);
      MAT_set_j(G,index4,index4);
      MAT_set_d(G,index4,1.);
      
      VEC_set(u,index1,FACTS_INF_Q);
      VEC_set(l,index1,-FACTS_INF_Q);

      VEC_set(u,index2,FACTS_INF_Q);
      VEC_set(l,index2,-FACTS_INF_Q);
      
      if (FACTS_has_flags(facts,FLAG_BOUNDED,FACTS_VAR_Q)) {
      	VEC_set(u,index3,FACTS_get_Q_max_s(facts));
      	VEC_set(l,index3,FACTS_get_Q_min_s(facts));
      	VEC_set(u,index4,FACTS_get_Q_max_sh(facts));
      	VEC_set(l,index4,FACTS_get_Q_min_sh(facts));
      }
      else {
      	VEC_set(u,index3,FACTS_INF_Q);
      	VEC_set(l,index3,-FACTS_INF_Q);
      	VEC_set(u,index4,FACTS_INF_Q);
      	VEC_set(l,index4,-FACTS_INF_Q);
      }
      
      // Row info
      CONSTR_set_G_row_info_string(c,
                                   index1,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "reactive power k",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index2,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "reactive power m",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index3,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "series reactive power",
                                   t);
      CONSTR_set_G_row_info_string(c,
                                   index4,
                                   "facts",
                                   FACTS_get_index(facts),
                                   "shunt reactive power",
                                   t);
    }
  }

}

void CONSTR_BOUND_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_BOUND_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Branch* br;
  Gen* gen;
  Load* load;
  Shunt* shunt;

  // Check
  if (!bus)
    return;

  // Voltage magnitude (V_MAG)
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
    BUS_set_sens_v_mag_u_bound(bus,VEC_get(sGu,BUS_get_index_v_mag(bus,t)),t);
    BUS_set_sens_v_mag_l_bound(bus,VEC_get(sGl,BUS_get_index_v_mag(bus,t)),t);
  }
  
  // Voltage angle (V_ANG)
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
    BUS_set_sens_v_ang_u_bound(bus,VEC_get(sGu,BUS_get_index_v_ang(bus,t)),t);
    BUS_set_sens_v_ang_l_bound(bus,VEC_get(sGl,BUS_get_index_v_ang(bus,t)),t);
  }
  
  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
    
    // Active power (P)
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
      GEN_set_sens_P_u_bound(gen,VEC_get(sGu,GEN_get_index_P(gen,t)),t);
      GEN_set_sens_P_l_bound(gen,VEC_get(sGl,GEN_get_index_P(gen,t)),t);
    }
    
    // Reactive power (Q)
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
      GEN_set_sens_Q_u_bound(gen,VEC_get(sGu,GEN_get_index_Q(gen,t)),t);
      GEN_set_sens_Q_l_bound(gen,VEC_get(sGl,GEN_get_index_Q(gen,t)),t);
    }
  }
  
  // Loads
  for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
    
    // Active power (P)
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
      LOAD_set_sens_P_u_bound(load,VEC_get(sGu,LOAD_get_index_P(load,t)),t);
      LOAD_set_sens_P_l_bound(load,VEC_get(sGl,LOAD_get_index_P(load,t)),t);
    }
  }
  
  // Variable generators
  
  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
    
    // Susceptance (b)
    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
      SHUNT_set_sens_b_u_bound(shunt,VEC_get(sGu,SHUNT_get_index_b(shunt,t)),t);
      SHUNT_set_sens_b_l_bound(shunt,VEC_get(sGl,SHUNT_get_index_b(shunt,t)),t);
    }
  }
  
  // Batteries

  // Branches
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {
  
    // Branch tap ratio
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
      BRANCH_set_sens_ratio_u_bound(br,VEC_get(sGu,BRANCH_get_index_ratio(br,t)),t);
      BRANCH_set_sens_ratio_l_bound(br,VEC_get(sGl,BRANCH_get_index_ratio(br,t)),t);
    }
    
    // Branch phase shift
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
      BRANCH_set_sens_phase_u_bound(br,VEC_get(sGu,BRANCH_get_index_phase(br,t)),t);
      BRANCH_set_sens_phase_l_bound(br,VEC_get(sGl,BRANCH_get_index_phase(br,t)),t);
    }
  }
}

