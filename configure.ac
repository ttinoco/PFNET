#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(PFNET, 1.2.5, ttinoco5687@gmail.com)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([tools/mpc2mat.m])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
# add libtools
LT_INIT
# LT_LANG([C++])


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

# Checks for optional libraries/packages.

# Check for Graphviz, optional.
#   option 1 use autoconf to check for specific libs/headers
#     i'm only checking one function in each library
#     this first checks if cgraph is available and if so then also check gvc
#     it will automatically add the libraries if successful, otherwise add no graphviz to cflags
no_graphviz=no
AC_CHECK_LIB([cgraph], [agopen], [], 
  [no_graphviz=yes])
AC_CHECK_LIB([gvc], [gvContext], [], 
  [no_graphviz=yes])
AS_IF([test "x$no_graphviz" = "xyes"], [CFLAGS+=" -DNO_GRAPHVIZ"], [])
AC_CHECK_HEADERS(graphviz/cgraph.h)
AC_CHECK_HEADERS(graphviz/gvc.h)
#   option 2 is to use pkg-config (assumes using pkg manager)
#PKG_PROG_PKG_CONFIG()
#PKG_CHECK_MODULES([GRAPHVIZ], [libcgraph >= 2.36 libgvc], [HAVE_GRAPHVIZ=1], [HAVE_GRAPVHIZ=0])

# Check for RAW_PARSER, optional.
# checks if raw-parser is availble and will automatically add to LIBS, otherwise add no raw parser to cflags
AC_MSG_WARN([raw_parser must be built with '-lpfnet' or be able to locate the unresolved symbols for libpfnet])
AC_CHECK_LIB([raw_parser], [RAW_PARSER_new], [], [CFLAGS+=" -DNO_RAW_PARSER"])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt strchr strdup strstr])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
