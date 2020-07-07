cmake_minimum_required(VERSION 3.5.0)
project(PFNET)
enable_testing()

set(PFNET_VERSION 1.3.6)
set(PFNET_PACKAGE "pfnet")
set(PFNET_NAME "PFNET")
set(PFNET_STRING "${PFNET_NAME} ${PFNET_VERSION}")
set(PFNET_BUGREPORT "awigington@epri.com")
set(PFNET_URL "https://github.com/romcon/PFNET")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)  # so that cmake builds with rpath already added
set(CMAKE_INSTALL_RPATH "@loader_path")  # for relocatable libraries

option(PFNET_DEBUG "set to ON to enable PFNET debug definition" OFF)
s
set(M_LIB,"")
if (UNIX)
  find_library(M_LIB m)
  set(CMAKE_REQUIRED_INCLUDES math.h)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif (UNIX)
s
# Find PFNET source files
file(GLOB_RECURSE pfnet_source src/*.c)

# Add shared library target
add_library(pfnet SHARED ${pfnet_source})

# Add static library target
add_library(pfnet_static STATIC ${pfnet_source})
set_property(TARGET pfnet_static PROPERTY POSITION_INDEPENDENT_CODE ON)

# PFNET tests
file(GLOB_RECURSE pfnet_test_source tests/*.c)

# Test dynamic library
add_executable(pfnet_tests ${pfnet_test_source})
add_test(run_pfnet_tests pfnet_tests ${PFNET_SOURCE_DIR}/data/ieee14.mat)
#target_link_libraries(pfnet_tests pfnet m)
target_link_libraries(pfnet_tests pfnet ${M_LIB})

# Test static library
add_executable(pfnet_static_tests ${pfnet_test_source})
add_test(run_pfnet_static_tests pfnet_static_tests ${PFNET_SOURCE_DIR}/data/ieee14.mat)
#target_link_libraries(pfnet_static_tests pfnet_static m)
target_link_libraries(pfnet_static_tests pfnet_static ${M_LIB})

# Set the debug flag
if(PFNET_DEBUG)
  add_definitions(-DDEBUG)
endif()

# Add raw-parser if available
if(EXISTS $ENV{RAW_PARSER}/src/RAW_parser.c)
  include_directories($ENV{RAW_PARSER}/include)
  set(raw_parser_source $ENV{RAW_PARSER}/src/RAW_net.c
                        $ENV{RAW_PARSER}/src/RAW_parser.c)
  target_sources(pfnet PUBLIC ${raw_parser_source})
  target_sources(pfnet_static PUBLIC ${raw_parser_source})
  message("raw-parser enabled.")
  # install raw-parser headers also
  install(DIRECTORY $ENV{RAW_PARSER}/include/
          DESTINATION include
          FILES_MATCHING PATTERN "*.h")
  set(HAVE_RAW_PARSER 1)  # need to set it to 1 so that the dummy parser in pfnet does not conflict
else()
  message("raw-parser not enabled.")
  set(HAVE_RAW_PARSER 0)
endif()

# Add epc-parser if available
if(EXISTS $ENV{EPC_PARSER}/src/EPC_parser.c)
  include_directories($ENV{EPC_PARSER}/include)
  set(epc_parser_source $ENV{EPC_PARSER}/src/EPC_net.c
                        $ENV{EPC_PARSER}/src/EPC_parser.c)
  target_sources(pfnet PUBLIC ${epc_parser_source})
  target_sources(pfnet_static PUBLIC ${epc_parser_source})
  message("epc-parser enabled.")
  # install epc-parser headers also
  install(DIRECTORY $ENV{EPC_PARSER}/include/
          DESTINATION include
          FILES_MATCHING PATTERN "*.h")
  set(HAVE_EPC_PARSER 1)  # need to set it to 1 so that the dummy parser in pfnet does not conflict
else()
  message("epc-parser not enabled.")
  set(HAVE_EPC_PARSER 0)
endif()

# Add checks needed to produce pfnet_config.h
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckTypeSize)
INCLUDE (CheckFunctionExists)
INCLUDE (CheckLibraryExists)

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
if (UNIX)
  CHECK_INCLUDE_FILES("unistd.h"  HAVE_UNISTD_H)
endif (UNIX)

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

# Create pfnet_config.h file
configure_file(${PFNET_SOURCE_DIR}/include/pfnet/pfnet_config.cmake.in ${PFNET_SOURCE_DIR}/include/pfnet/pfnet_config.h)

# Add include directory
include_directories(${PFNET_SOURCE_DIR}/include)

# Install library to <prefix>/lib
install(TARGETS pfnet pfnet_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION lib)
# Install header files
install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
