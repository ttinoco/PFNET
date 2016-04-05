/** @file constr_LBOUND.c
 *  @brief This file defines the data structure and routines associated with the constraint of type LBOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_LBOUND.h>

void CONSTR_LBOUND_init(Constr* c) {

  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_LBOUND_clear(Constr* c) {
  
  // Counters
  CONSTR_set_Gcounter(c,0);
  CONSTR_set_Gconstr_index(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_LBOUND_count_branch(Constr* c, Branch* br) {
  // Nothing
}

void CONSTR_LBOUND_allocate(Constr *c) {
  
  // Local variables
  int num_vars;
  
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // l u G
  CONSTR_set_l(c,VEC_new(num_vars));
  CONSTR_set_u(c,VEC_new(num_vars));
  CONSTR_set_G(c,MAT_new(num_vars,   // size1 (rows)
			 num_vars,   // size2 (cols)
			 num_vars)); // nnz
}

void CONSTR_LBOUND_analyze_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Load* load;
  Vargen* vargen;
  Shunt* shunt;
  char* bus_counted;
  Vec* l;
  Vec* u;
  Mat* G;
  int i;
  int index;
  int index1;
  int index2;
  
  // Cosntr data
  l = CONSTR_get_l(c);
  u = CONSTR_get_u(c);
  G = CONSTR_get_G(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Tap ratio 
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    index = BRANCH_get_index_ratio(br);
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
  }

  // Tap ratio dev
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) {
    index1 = BRANCH_get_index_ratio_y(br);
    index2 = BRANCH_get_index_ratio_z(br);
    MAT_set_i(G,index1,index1);
    MAT_set_j(G,index1,index1);    
    MAT_set_d(G,index1,1.);
    MAT_set_i(G,index2,index2);
    MAT_set_j(G,index2,index2);    
    MAT_set_d(G,index2,1.);
    VEC_set(u,index1,BRANCH_INF_RATIO);
    VEC_set(u,index2,BRANCH_INF_RATIO);
    VEC_set(l,index1,0.);
    VEC_set(l,index2,0.);
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    index = BRANCH_get_index_phase(br);
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
  }

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	index = BUS_get_index_v_mag(bus);
	MAT_set_i(G,index,index);
	MAT_set_j(G,index,index);    
	MAT_set_d(G,index,1.);
	if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG)) {
	  VEC_set(u,index,BUS_get_v_max(bus));     
	  VEC_set(l,index,BUS_get_v_min(bus));
	}
	else {
	  VEC_set(u,index,BUS_INF_V_MAG);
	  VEC_set(l,index,0.);
	}
      }

      // Voltage magnitude dev (V_DEV)
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	index1 = BUS_get_index_y(bus);
	index2 = BUS_get_index_z(bus);
	MAT_set_i(G,index1,index1);
	MAT_set_j(G,index1,index1);    
	MAT_set_d(G,index1,1.);
	MAT_set_i(G,index2,index2);
	MAT_set_j(G,index2,index2);    
	MAT_set_d(G,index2,1.);
	VEC_set(u,index1,BUS_INF_V_MAG);
	VEC_set(u,index2,BUS_INF_V_MAG);
	VEC_set(l,index1,0.);
	VEC_set(l,index2,0.);
      }

      // Voltage magnitude vio (V_VIO)
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	index1 = BUS_get_index_vl(bus);
	index2 = BUS_get_index_vh(bus);
	MAT_set_i(G,index1,index1);
	MAT_set_j(G,index1,index1);    
	MAT_set_d(G,index1,1.);
	MAT_set_i(G,index2,index2);
	MAT_set_j(G,index2,index2);    
	MAT_set_d(G,index2,1.);
	VEC_set(u,index1,BUS_INF_V_MAG);
	VEC_set(u,index2,BUS_INF_V_MAG);
	VEC_set(l,index1,0.);
	VEC_set(l,index2,0.);
      }
      
      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	index = BUS_get_index_v_ang(bus);
	MAT_set_i(G,index,index);
	MAT_set_j(G,index,index);    
	MAT_set_d(G,index,1.);
	VEC_set(u,index,2.*PI);
	VEC_set(l,index,-2.*PI);
      }
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  index = GEN_get_index_P(gen);
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
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  index = GEN_get_index_Q(gen);
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
	}	
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	
	// Active power (P)
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  index = LOAD_get_index_P(load);
	  MAT_set_i(G,index,index);
	  MAT_set_j(G,index,index);    
	  MAT_set_d(G,index,1.);
	  if (LOAD_has_flags(load,FLAG_BOUNDED,LOAD_VAR_P)) {
	    VEC_set(u,index,LOAD_get_P_max(load));     
	    VEC_set(l,index,LOAD_get_P_min(load));
	  }
	  else {
	    VEC_set(u,index,LOAD_INF_P);
	    VEC_set(l,index,-LOAD_INF_P);
	  }
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  index = VARGEN_get_index_P(vargen);
	  MAT_set_i(G,index,index);
	  MAT_set_j(G,index,index);    
	  MAT_set_d(G,index,1.);
	  if (VARGEN_has_flags(vargen,FLAG_BOUNDED,VARGEN_VAR_P)) {
	    VEC_set(u,index,VARGEN_get_P_max(vargen));     
	    VEC_set(l,index,VARGEN_get_P_min(vargen));
	  }
	  else {
	    VEC_set(u,index,VARGEN_INF_P);
	    VEC_set(l,index,-VARGEN_INF_P);
	  }
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  index = VARGEN_get_index_Q(vargen);
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
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  index = SHUNT_get_index_b(shunt);
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
	}

	// Susceptance dev 
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) {
	  index1 = SHUNT_get_index_y(shunt);
	  index2 = SHUNT_get_index_z(shunt);
	  MAT_set_i(G,index1,index1);
	  MAT_set_j(G,index1,index1);    
	  MAT_set_d(G,index1,1.);
	  MAT_set_i(G,index2,index2);
	  MAT_set_j(G,index2,index2);    
	  MAT_set_d(G,index2,1.);
	  VEC_set(u,index1,SHUNT_INF_SUSC);
	  VEC_set(u,index2,SHUNT_INF_SUSC);
	  VEC_set(l,index1,0.);
	  VEC_set(l,index2,0.);
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)] = TRUE;    
  }
}

void CONSTR_LBOUND_eval_branch(Constr* c, Branch* br, Vec* var_values) {
  // Nothing to do
}

void CONSTR_LBOUND_store_sens_branch(Constr* c, Branch* br, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Load* load;
  char* bus_counted;
  int i;
  
  // Constr data
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointer
  if (!bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)]) {
      
      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	BUS_set_sens_v_ang_u_bound(bus,VEC_get(sGu,BUS_get_index_v_ang(bus)));
	BUS_set_sens_v_ang_l_bound(bus,VEC_get(sGl,BUS_get_index_v_ang(bus)));
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  GEN_set_sens_P_u_bound(gen,VEC_get(sGu,GEN_get_index_P(gen)));
	  GEN_set_sens_P_l_bound(gen,VEC_get(sGl,GEN_get_index_P(gen)));
	}
      } 

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {
	
	// Active power (P)
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  LOAD_set_sens_P_u_bound(load,VEC_get(sGu,LOAD_get_index_P(load)));
	  LOAD_set_sens_P_l_bound(load,VEC_get(sGl,LOAD_get_index_P(load)));
	}
      } 
    }
    
    // Update counted flag
    bus_counted[BUS_get_index(bus)] = TRUE;    
  }
}

void CONSTR_LBOUND_free(Constr* c) {
  // Nothing to do
}
