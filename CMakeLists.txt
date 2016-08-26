cmake_minimum_required(VERSION 3.5.0)
project(PFNET)

# example cmake usage from raw-parser directory
# $ mkdir build
# $ cd build
# $ cmake -DCMAKE_INSTALL_PREFIX=../.. -DRAW_PARSER_INSTALL_PREFIX=../.. ..
# $ make
# $ ls # see libpfnet.{a,so,dylib} depending on configuration

# pfnet debug mode
set(PFNET_DEBUG OFF CACHE BOOL "set to ON to enable PFNET debug definition")

# RAW_PARSER_INSTALL_PREFIX
set(RAW_PARSER_INSTALL_PREFIX ""
    CACHE FILEPATH
    "path to raw-parser install prefix")

if(NOT ${RAW_PARSER_INSTALL_PREFIX} STREQUAL "")
  list(APPEND CMAKE_PREFIX_PATH ${RAW_PARSER_INSTALL_PREFIX})
endif()

# find source files
file(GLOB_RECURSE pfnet_source src/*.c)

# add include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# add library target
add_library(pfnet SHARED ${pfnet_source})

# set the debug flag
if(PFNET_DEBUG)
  add_definitions(-DDEBUG)
endif()

# find graphviz
find_library(GRAPHVIZ_LIB gvc)
if(GRAPHVIZ_LIB)
  message("Graphiz found: " ${GRAPHVIZ_LIB})
  target_link_libraries(pfnet gvc cgraph)
else()
  message("Graphiz not found.")
  add_definitions(-DNO_GRAPHVIZ)
endif()

# find raw_parser
find_library(RAW_PARSER_LIB libraw_parser_static.a)
if(RAW_PARSER_LIB)
  message("raw-parser found: " ${RAW_PARSER_LIB})
  target_link_libraries(pfnet ${RAW_PARSER_LIB})
else()
  message("raw-parser not found.")
  add_definitions(-DNO_RAW_PARSER)
endif()

# install library to <prefix>/lib
install(TARGETS pfnet
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
## copy header files to <prefix>/include
install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
