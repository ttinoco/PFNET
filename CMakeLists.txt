cmake_minimum_required(VERSION 3.5.0)
project(PFNET)

# example cmake usage from raw-parser directory
# $ mkdir build
# $ cd build
# $ cmake -DCMAKE_INSTALL_PREFIX=../.. ..
# $ make
# $ ls # see libpfnet.{a,so,dylib} depending on configuration

# find source files
file(GLOB_RECURSE pfnet_source src/*.c)

# list source files
#message(pfnet_source: "${pfnet_source}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# add library target
add_library(pfnet SHARED ${pfnet_source})

# find graphviz
find_library(GRAPHVIZ gvc)
if ("${GRAPHVIZ}" STREQUAL "GRAPHVIZ-NOTFOUND")
  message("Graphiz not found.")
  add_definitions(-DNO_GRAPHVIZ)
else()
  target_link_libraries(pfnet gvc cgraph)
endif()

# find raw_parser

find_library(RAW_PARSER raw_parser)
if ("${RAW_PARSER}" STREQUAL "RAW_PARSER-NOTFOUND")
 message("raw-parser not found.")
 add_definitions(-DNO_RAW_PARSER)
else()
 message("raw-parser found.")
endif()

# install library to <prefix>/lib
install(TARGETS pfnet
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
## copy header files to <prefix>/include
install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
