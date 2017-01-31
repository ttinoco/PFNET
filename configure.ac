#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([PFNET], [1.2.5], [ttinoco5687@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/net/net.c])
AC_CONFIG_HEADERS([include/pfnet/config.h])
AC_CONFIG_MACRO_DIR([m4])

# Add libtools
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for graphviz
AC_CHECK_LIB(gvc, gvContext)
AC_CHECK_LIB(cgraph, agopen)
AC_CHECK_HEADERS(graphviz/gvc.h)

# Checks for raw parser
AC_CHECK_FILES([$RAW_PARSER/lib/libraw_parser.so],
	[AC_DEFINE([HAVE_RAW_PARSER],[1],[Define to 1 if you have the raw parser library.]),
	 LDFLAGS="$LDFLAGS -lraw_parser"])

# Checks for other header files.
AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt strchr strdup strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
