cmake_minimum_required (VERSION 2.6.0)
project (fadbad)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -g -Wall -Wno-sign-compare")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")

#set(MY_LIBRARIES "dl")

set(FADBAD_INCLUDE_DIR "/home/gkahn/source/FADBAD++")

find_package(Armadillo)
find_package(Eigen)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${FADBAD_INCLUDE_DIR} ${Eigen_INCLUDE_DIRS})
link_directories(${ARMADILLO_LIBRARIES})

add_executable(test-forward test-forward.cpp)

add_executable(test-backward test-backward.cpp)

add_executable(test-branching test-branching.cpp)

add_executable(test-with-arma test-with-arma.cpp)
target_link_libraries(test-with-arma ${ARMADILLO_LIBRARIES})

add_executable(test-with-eigen test-with-eigen.cpp)
