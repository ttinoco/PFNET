cmake_minimum_required(VERSION 3.5.0)
project(PFNET)
enable_testing()

# example cmake usage from raw-parser directory
# $ mkdir build
# $ cd build
# $ cmake -DCMAKE_INSTALL_PREFIX=../.. -DRAW_PARSER_SOURCE_DIR=../../raw-parser ..
# $ make
# $ ls # see libpfnet.{so,dylib} depending on configuration

option(PFNET_DEBUG "set to ON to enable PFNET debug definition" OFF)
option(PFNET_GRAPHVIZ "set to OFF to disable search for graphviz" ON)

# RAW_PARSER_INSTALL_PREFIX
set(RAW_PARSER_SOURCE_DIR ""
    CACHE FILEPATH
    "path to raw-parser source directory")

# find PFNET source files
file(GLOB_RECURSE pfnet_source src/*.c)

# add include directory
include_directories(${PFNET_SOURCE_DIR}/include)

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
target_link_libraries(pfnet_tests pfnet m)

# test static library
add_executable(pfnet_static_tests ${pfnet_test_source})
add_test(run_pfnet_static_tests pfnet_static_tests ${PFNET_SOURCE_DIR}/data/ieee14.mat)
target_link_libraries(pfnet_static_tests pfnet_static m)

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
  add_definitions(-DNO_GRAPHVIZ)
endif()

# add raw-parser if specified
if(NOT ${RAW_PARSER_SOURCE_DIR} STREQUAL "")
  include_directories(${RAW_PARSER_SOURCE_DIR}/include)
  set(raw_parser_source ${RAW_PARSER_SOURCE_DIR}/src/RAW_net.c
                        ${RAW_PARSER_SOURCE_DIR}/src/RAW_parser.c)
  target_sources(pfnet PUBLIC ${raw_parser_source})
  target_sources(pfnet_static PUBLIC ${raw_parser_source})
  message("raw-parser enabled.")
  # install raw-parser headers also
  install(DIRECTORY ${RAW_PARSER_SOURCE_DIR}/include/
          DESTINATION include
          FILES_MATCHING PATTERN "*.h")
else()
  message("raw-parser not enabled.")
  add_definitions(-DNO_RAW_PARSER)
endif()

# install library to <prefix>/lib
install(TARGETS pfnet pfnet_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
# install header files
install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
