/** @file heur_REG_PF_SWITCH.c
 *  @brief This file defines the data structure and routines associated with the heuristic of type REG_PF_SWITCH.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/heur_REG_PF_SWITCH.h>
#include <pfnet/constr_REG_PF_SWITCH.h>

Heur* HEUR_REG_PF_SWITCH_new(Net* net) {
  Heur* h = HEUR_new(net);
  HEUR_set_func_apply_step(h,&HEUR_REG_PF_SWITCH_apply_step);
  HEUR_set_name(h,"switching power factor regulation");
  HEUR_init(h);
  return h;
}


void HEUR_REG_PF_SWITCH_apply_step(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Net* net;
  Vec* f;
  ConvVSC* vsc;
  char* fix_flag;
  Constr* cpf;
  Constr* creg;  
  REAL Q;
  REAL Qmax;
  REAL Qmin;
  int i;

  // Heur data
  net = HEUR_get_network(h);
  
  // Power flow constraint
  cpf = NULL;
  for (i = 0; i < cnum; i++) {
    if (strcmp(CONSTR_get_name(cptrs[i]),"AC power balance") == 0) {
      cpf = cptrs[i];
      break;
    }
  }
  if (!cpf) {
    HEUR_set_error(h, "unable to find AC power balance constraint");
    return;
  }

  // Switching power factor regulation constraint
  creg = NULL;
  for (i = 0; i < cnum; i++) {
    if (strcmp(CONSTR_get_name(cptrs[i]),"switching power factor regulation") == 0) {
      creg = cptrs[i];
      break;
    }
  }
  if (!creg) {
    HEUR_set_error(h, "unable to find switching power factor regulation constraints");
    return;
  }

  // Fix flags
  fix_flag = CONSTR_REG_PF_SWITCH_get_flags(creg);
  if (!fix_flag)
    return;
  
  // Constr data
  f = CONSTR_get_f(cpf);

  // VSC
  for(vsc = BUS_get_vsc_conv(bus); vsc != NULL; vsc = CONVVSC_get_next_ac(vsc)) {
        
    if (CONVVSC_is_in_f_ac_mode(vsc) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_P) &&
        CONVVSC_has_flags(vsc,FLAG_VARS,CONVVSC_VAR_Q) &&
        CONVVSC_is_in_service(vsc)) {
          
      Q = VEC_get(var_values,CONVVSC_get_index_Q(vsc,t)); // per unit
      Qmax = CONVVSC_get_Q_max(vsc);                      // per unit
      Qmin = CONVVSC_get_Q_min(vsc);                      // per unit
          
      // CASE: Q free
      if (!fix_flag[CONVVSC_get_index_Q(vsc,t)]) {
          
        // Qmax violation
        if (Q >= Qmax) {
            
          // Switch to Q fixed
          CONVVSC_set_Q(vsc,Qmax,t);
          fix_flag[CONVVSC_get_index_Q(vsc,t)] = TRUE;
          if (HEUR_REG_PF_SWITCH_DEBUG)
            printf("HEUR REG_PF_SWITCH: fix vsc at Qmax\n");
        }
          
        // Q is free: Qmin violation
        else if (Q <= Qmin) {
              
          // Switch to Q fixed
          CONVVSC_set_Q(vsc,Qmin,t);
          fix_flag[CONVVSC_get_index_Q(vsc,t)] = TRUE;
          if (HEUR_REG_PF_SWITCH_DEBUG)
            printf("HEUR REG_PF_SWITCH: fix vsc at Qmin\n");
        } 
      }
        
      // CASE: Q fixed
      else {
                            
        // Q at Qmin - check if free Q moves up
        if (fabs(Q-Qmin) < fabs(Q-Qmax)) {
          Q = Q - VEC_get(f,BUS_get_dQ_index(CONVVSC_get_ac_bus(vsc),t)); // per unit (see constr_ACPF)
          if (Qmin < Q) {
            fix_flag[CONVVSC_get_index_Q(vsc,t)] = FALSE; // Switch to Q free
            if (HEUR_REG_PF_SWITCH_DEBUG)
              printf("HEUR REG_PF_SWITCH: free vsc from Qmin\n");		
          }
        }
            
        // Q at Qmax - check if free Q moves down
        else if (fabs(Q-Qmax) < fabs(Q-Qmin)) { 
          Q = Q - VEC_get(f,BUS_get_dQ_index(CONVVSC_get_ac_bus(vsc),t)); // per unit (see constr_ACPF)
          if (Q < Qmax) {
            fix_flag[CONVVSC_get_index_Q(vsc,t)] = FALSE; // Switch to Q free
            if (HEUR_REG_PF_SWITCH_DEBUG)
              printf("HEUR REG_PF_SWITCH: free vsc from Qmax\n");
          }
        }
      }
    }
  }
  
  // Update
  if (BUS_get_index(bus) == NET_get_num_buses(net,FALSE)-1 && t == NET_get_num_periods(net)-1)
    CONSTR_analyze(creg);
}
