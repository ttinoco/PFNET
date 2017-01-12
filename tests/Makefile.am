#AM_CFLAGS+=-DSRCDIR=$(top_srcdir)

#EXTRA_DIST = $(top_srcdir)/data/ieee14.mat

# define program to build with 'check'
check_PROGRAMS = run_tests

# sources for the program(s)
## run_tests
run_tests_headers = unit.h test_network.h test_problem.h test_constraints.h test_graph.h

run_tests_SOURCES = run_tests.c $(run_tests_headers)

run_tests_CFLAGS = -I$(top_builddir)/include
# run_tests_LDFLAGS = -shared  # want both static and shared
run_tests_LDADD = -L$(top_srcdir)/lib $(top_builddir)/src/libpfnet.la

# define tests to run with 'check'
TESTS = run_tests.sh
