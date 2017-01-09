#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(PFNET, 1.2.4, ttinoco5687@gmail.com)
AM_INIT_AUTOMAKE([foreign])
LT_INIT()
AC_CONFIG_SRCDIR([tools/mpc2mat.m])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
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

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
