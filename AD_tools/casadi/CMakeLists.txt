cmake_minimum_required (VERSION 2.6.0)
project (casadi)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -g -Wall -Wno-sign-compare")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

set(CASADI_INCLUDE_DIR "$ENV{CASADIPATH}")
set(CASADI_LIBRARY_DIR "$ENV{CASADIPATH}/build/lib")
set(CASADI_LIBRARY "casadi")

#add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")

set(MY_LIBRARIES "dl")

message("CASADI_INCLUDE_DIR: ${CASADI_INCLUDE_DIR}")
message("CASADI_LIBRARY_DIR: ${CASADI_LIBRARY_DIR}")
message("CASADI_LIBRARY: ${CASADI_LIBRARY}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CASADI_INCLUDE_DIR})
link_directories(${CASADI_LIBRARY_DIR})

add_executable(test-inequality test-inequality.cpp)
target_link_libraries(test-inequality ${CASADI_LIBRARY} ${MY_LIBRARIES})
