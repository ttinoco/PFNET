#AM_CFLAGS+=-DSRCDIR=$(top_srcdir)

# define program to build with 'check'
check_PROGRAMS = run_tests

# sources for the program(s)
## run_tests
run_tests_headers = unit.h test_network.h test_problem.h test_constraints.h test_functions.h test_graph.h

run_tests_SOURCES = run_tests.c $(run_tests_headers)

run_tests_CFLAGS = -I$(top_srcdir)/include

run_tests_LDADD = $(top_builddir)/src/libpfnet.la

# define tests to run with 'check'
TESTS = run_tests.sh

EXTRA_DIST = run_tests.sh
