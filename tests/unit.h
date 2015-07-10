/** @file unit.h
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#define Assert(message, test) do { if (!(test)) {printf("FAIL\n"); return message;} } while (0)
#define run_test(test) do { char *message = test(); tests_run++; \
                            if (message) return message; } while (0)
extern int tests_run;
extern char* test_case;
