/** @file run_tests.c
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include "unit.h"
#include "test_network.h"
#include "test_graph.h"
#include "test_constraints.h"
#include "test_functions.h"
#include "test_problem.h"

int tests_run = 0;
char* test_case = NULL;

static char * all_tests() {

  // Network
  run_test(test_net_new);
  run_test(test_net_load);
  run_test(test_net_check);
  run_test(test_net_variables);
  run_test(test_net_fixed);
  run_test(test_net_properties);
  run_test(test_net_init_point);

  // Graph
  run_test(test_graph_basic);

  // Constraints
  run_test(test_constr_NBOUND);
  run_test(test_constr_FIX);
  run_test(test_constr_PAR_GEN_P);
  run_test(test_constr_PAR_GEN_Q);
  run_test(test_constr_ACPF);
  run_test(test_constr_REG_GEN);
  run_test(test_constr_REG_TRAN);
  run_test(test_constr_REG_SHUNT);

  // Functions
  run_test(test_func_GEN_COST);

  // Problem
  run_test(test_problem_basic);
  
  return 0;
}

int main(int argc, char **argv) {

  // Local variables
  char* result;

  // Check inputs
  if( argc < 2) {
    printf("usage: run_tests test_case\n");
    return -1;
  }
  
  // Get case
  test_case = argv[1];

  // Run tests
  result = all_tests();

  // Show results
  if (result != 0)
    printf("%s\n", result);
  else
    printf("All tests passed.\n");
  printf("Tests run: %d\n", tests_run);

  // Return outcome
  return result != 0;
}

