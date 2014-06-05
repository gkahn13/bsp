cmake_minimum_required (VERSION 2.6.0)
project (pf)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(PythonLibs 2.7 REQUIRED)

find_package(OpenRAVE REQUIRED)
#find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams python thread filesystem system program_options)
find_package(Boost COMPONENTS iostreams python thread filesystem system program_options)

find_package(Armadillo)
find_package(LAPACK)
find_package(BLAS)

set(CASADI_INCLUDE_DIR "$ENV{CASADIPATH}")
set(CASADI_LIBRARY_DIR "$ENV{CASADIPATH}/build/lib")
set(CASADI_LIBRARY "casadi")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")
add_definitions(${CMAKE_CXX_FLAGS} "-Wno-sign-compare")

#add_definitions(${CMAKE_C_FLAGS} "-std=gnu99")
add_definitions(${CMAKE_C_FLAGS} "-Wno-sign-compare")

add_definitions(${CMAKE_CXX_FLAGS_DEBUG} "-g")
#add_definitions(${CMAKE_CXX_FLAGS_DEBUG} "-Wall")
add_definitions(${CMAKE_C_FLAGS_DEBUG} "-g")
#add_definitions(${CMAKE_C_FLAGS_DEBUG} "-Wall")

set(MY_LIBRARIES "dl;rt;python2.7")

message("CASADI_INCLUDE_DIR: ${CASADI_INCLUDE_DIR}")
message("CASADI_LIBRARY_DIR: ${CASADI_LIBRARY_DIR}")
message("CASADI_LIBRARY: ${CASADI_LIBRARY}")

message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

message("LAPACK_FOUND: ${LAPACK_FOUND}")
message("BLAS FOUND: ${BLAS_FOUND}")

message("BOOST LIBRARIES: ${Boost_LIBRARIES}")
message("PYTHON LIBRARIES: ${PYTHON_LIBRARIES}")
message("PYTHON INCLUDES: ${PYTHON_INCLUDE_DIRS}")
message("PYTHON VERSION: ${PYTHONLIBS_VERSION_STRING}")

include_directories(${OpenRAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS} ${CASADI_INCLUDE_DIR})
if( Boost_INCLUDE_DIRS )
  include_directories(${Boost_INCLUDE_DIRS})
endif()

link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${PYTHON_LIBRARIES} ${CASADI_LIBRARY_DIR}
                 ${ARMADILLO_LIBRARIES})

add_executable(test-eih eih/tests/test.cpp 
               system.cpp casadi-system.cpp
               eih/src/eih_system.cpp eih/src/pr2_sim.cpp eih/src/rave_utils.cpp
               ../util/logging.cpp)
set_target_properties(test-eih PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-eih PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-eih ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES}
                                   ${ARMADILLO_LIBRARIES} ${CASADI_LIBRARY} ${PYTHON_LIBRARIES})

add_executable(eih eih/eih.cpp eih/eihMPC.c
               system.cpp casadi-system.cpp
               eih/src/eih_system.cpp eih/src/pr2_sim.cpp eih/src/rave_utils.cpp
               ../util/logging.cpp)
set_target_properties(eih PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(eih PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(eih ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES}
                          ${ARMADILLO_LIBRARIES} ${CASADI_LIBRARY} ${PYTHON_LIBRARIES} ${MY_LIBRARIES})

add_executable(explore system.cpp casadi-system.cpp
               explore/explore-system.cpp explore/casadi-explore-system.cpp
               explore/explore.cpp explore/exploreMPC.c
               ../util/logging.cpp)
target_link_libraries(explore ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${CASADI_LIBRARY} ${ARMADILLO_LIBRARIES} ${MY_LIBRARIES})

add_executable(boxes system.cpp casadi-system.cpp
               boxes/boxes-system.cpp boxes/casadi-boxes-system.cpp
               boxes/boxes.cpp boxes/boxesMPC.c
               ../util/logging.cpp)
target_link_libraries(boxes ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${CASADI_LIBRARY} ${ARMADILLO_LIBRARIES} ${MY_LIBRARIES})

add_executable(test-entropy tests/test-entropy.cpp ../util/logging.cpp)
target_link_libraries(test-entropy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MY_LIBRARIES})