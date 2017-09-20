cmake_minimum_required(VERSION 3.5.0)
project(PFNET)
enable_testing()

set(PFNET_VERSION 1.3.1)
set(PFNET_PACKAGE "pfnet")
set(PFNET_NAME "PFNET")
set(PFNET_STRING "${PFNET_NAME} ${PFNET_VERSION}")
set(PFNET_BUGREPORT "ttinoco5687@gmail.com")
set(PFNET_URL "https://github.com/ttinoco/PFNET")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)  # so that cmake builds with rpath already added
set(CMAKE_INSTALL_RPATH "@loader_path")  # for relocatable libraries

option(PFNET_DEBUG "set to ON to enable PFNET debug definition" OFF)
option(PFNET_GRAPHVIZ "set to ON to enable graphviz addon" ON)
option(PFNET_LINE_FLOW "set to ON to enable line flow addon" ON)

set(M_LIB,"")
# REQUIRE MATH LIBRARY (NOT SURE IF OS SPECIFIC)
if (UNIX)
  find_library(M_LIB m)
  set(CMAKE_REQUIRED_INCLUDES math.h)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif (UNIX)
# CHECK_INCLUDE_FILES("math.h" HAVE_MATH_H)

# RAW_PARSER_INSTALL_PREFIX
#set(RAW_PARSER_SOURCE_DIR ""
#    CACHE FILEPATH
#    "path to raw-parser source directory")

# find PFNET source files
file(GLOB_RECURSE pfnet_source src/*.c)

# add shared library target
add_library(pfnet SHARED ${pfnet_source})

# add static library target
add_library(pfnet_static STATIC ${pfnet_source})
set_property(TARGET pfnet_static PROPERTY POSITION_INDEPENDENT_CODE ON)

# pfnet tests
file(GLOB_RECURSE pfnet_test_source tests/*.c)

# test dynamic library
add_executable(pfnet_tests ${pfnet_test_source})
add_test(run_pfnet_tests pfnet_tests ${PFNET_SOURCE_DIR}/data/ieee14.mat)
#target_link_libraries(pfnet_tests pfnet m)
target_link_libraries(pfnet_tests pfnet ${M_LIB})

# test static library
add_executable(pfnet_static_tests ${pfnet_test_source})
add_test(run_pfnet_static_tests pfnet_static_tests ${PFNET_SOURCE_DIR}/data/ieee14.mat)
#target_link_libraries(pfnet_static_tests pfnet_static m)
target_link_libraries(pfnet_static_tests pfnet_static ${M_LIB})

# set the debug flag
if(PFNET_DEBUG)
  add_definitions(-DDEBUG)
endif()

# find graphviz
if(PFNET_GRAPHVIZ)
  find_library(GRAPHVIZ_LIB gvc)
else()
  set(GRAPHVIZ_LIB "")
endif()

if(GRAPHVIZ_LIB)
  message("Graphiz found: " ${GRAPHVIZ_LIB})
  target_link_libraries(pfnet gvc cgraph)
  target_link_libraries(pfnet_tests gvc cgraph)
  target_link_libraries(pfnet_static_tests gvc cgraph)
else()
  message("Graphiz not enabled.")
endif()

# add raw-parser if specified
#if(NOT ${RAW_PARSER_SOURCE_DIR} STREQUAL "")
#  include_directories(${RAW_PARSER_SOURCE_DIR}/include)
#  set(raw_parser_source ${RAW_PARSER_SOURCE_DIR}/src/RAW_net.c
#                        ${RAW_PARSER_SOURCE_DIR}/src/RAW_parser.c)
#  target_sources(pfnet PUBLIC ${raw_parser_source})
#  target_sources(pfnet_static PUBLIC ${raw_parser_source})
#  message("raw-parser enabled.")
#  # install raw-parser headers also
#  install(DIRECTORY ${RAW_PARSER_SOURCE_DIR}/include/
#          DESTINATION include
#          FILES_MATCHING PATTERN "*.h")
#  set(HAVE_RAW_PARSER 1)  # need to set it to 1 so that the dummy parser in pfnet does not conflict
#else()
#  message("raw-parser not enabled.")
#  set(HAVE_RAW_PARSER 0)
#endif()

# add line-flow if specified
# if(PFNET_LINE_FLOW)
#   message("line-flow enabled.")
#   set(HAVE_LINE_FLOW 1 CACHE INTERNAL "")
#   # find_library(LINE_FLOW_LIB ?)  # TODO check with tomas
# else()
#   message("line-flow not enabled.")
#   set(HAVE_LINE_FLOW 0 CACHE INTERNAL "")
# endif()

# add checks needed to produce pfnet_config.h
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckTypeSize)
INCLUDE (CheckFunctionExists)
INCLUDE (CheckLibraryExists)

# Check for include files
# Checks for graphviz
CHECK_LIBRARY_EXISTS(gvc gvContext "" HAVE_LIBGVC)
CHECK_LIBRARY_EXISTS(cgraph agopen "" HAVE_LIBCGRAPH)
CHECK_INCLUDE_FILES("graphviz/gvc.h"  HAVE_GRAPHVIZ_GVC_H)

# check for line-flow
CHECK_LIBRARY_EXISTS(line_flow LINE_FLOW_new "${LINE_FLOW}/lib/" HAVE_LINE_FLOW)

# Checks std header files.
CHECK_INCLUDE_FILES("stddef.h"  HAVE_STDDEF_H)
CHECK_INCLUDE_FILES("stdint.h"  HAVE_STDINT_H)
CHECK_INCLUDE_FILES("stdlib.h"  HAVE_STDLIB_H)

# Checks other header files
CHECK_INCLUDE_FILES("string.h"  HAVE_STRING_H)
CHECK_INCLUDE_FILES("strings.h"  HAVE_STRING_H)
CHECK_INCLUDE_FILES("memory.h"  HAVE_MEMORY_H)
CHECK_INCLUDE_FILES("sys/stat.h"  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES("sys/types.h"  HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES("unistd.h"  HAVE_UNISTD_H)

# Checks for raw parser, not needed since adding files directly to package above
# Option "-Wl,-rpath,$PWD/lib" needed for "make check" without LD_LIBRARY_PATH
CHECK_LIBRARY_EXISTS(raw_parser RAW_PARSER_new "${RAW_PARSER}/lib/" HAVE_RAW_PARSER)

# Checks for typedefs, structures, and compiler characteristics.
CHECK_INCLUDE_FILES("inttypes.h" HAVE_INTTYPES_H)
CHECK_TYPE_SIZE(ptrdiff_t HAVE_PTRDIFF_T)

# Checks for library functions.v
CHECK_FUNCTION_EXISTS(malloc HAVE_MALLOC)
CHECK_FUNCTION_EXISTS(memset HAVE_MEMSET)
CHECK_FUNCTION_EXISTS(pow HAVE_POW)
CHECK_FUNCTION_EXISTS(sqrt HAVE_SQRT)
CHECK_FUNCTION_EXISTS(strchr HAVE_STRCHR)
CHECK_FUNCTION_EXISTS(strdup HAVE_STRDUP)
CHECK_FUNCTION_EXISTS(strstr HAVE_STRSTR)

# create pfnet_config.h file
configure_file(${PFNET_SOURCE_DIR}/include/pfnet/pfnet_config.cmake.in ${PFNET_SOURCE_DIR}/include/pfnet/pfnet_config.h)

# add include directory
include_directories(${PFNET_SOURCE_DIR}/include)

# install library to <prefix>/lib
install(TARGETS pfnet pfnet_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION lib)
# install header files
install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
