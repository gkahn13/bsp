cmake_minimum_required (VERSION 2.6.0)
project (pr2_eih)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall -Wno-sign-compare -Wno-enum-compare -Wno-reorder -std=c++11")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(catkin REQUIRED COMPONENTS pr2_utils roscpp roslib pcl_utils actionlib tf tf_conversions)
find_package(OpenRAVE REQUIRED)
find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS system program_options)
find_package(Eigen)

message("catkin_INCLUDE_DIRS: ${catkin_INCLUDE_DIRS}")
message("catkin_LIBRARY_DIRS: ${catkin_LIBRARY_DIRS}")
message("catkin_LIBRARIES: ${catkin_LIBRARIES}")

include_directories(${catkin_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
                    ${OpenRAVE_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS})

link_directories(${catkin_LIBRARY_DIRS} ${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})

add_executable(pr2-eih-sim  pr2_eih_sim.cpp system/pr2_eih_system.cpp sqp/pr2_eih_sqp.cpp mpc/pr2eihMPC.c ../util/logging.cpp)
set_target_properties(pr2-eih-sim PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(pr2-eih-sim PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(pr2-eih-sim ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(bsp  bsp.cpp system/pr2_eih_system.cpp system/pr2_eih_mapping.cpp
                                    sqp/pr2_eih_sqp.cpp mpc/pr2eihMPC.c ../util/logging.cpp)
set_target_properties(bsp PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(bsp PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(bsp ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(nbv  nbv.cpp system/pr2_eih_system.cpp system/pr2_eih_mapping.cpp
                                    sqp/pr2_eih_sqp.cpp mpc/pr2eihMPC.c ../util/logging.cpp)
set_target_properties(nbv PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(nbv PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(nbv ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test-pr2-eih-system  tests/test_pr2_eih_system.cpp system/pr2_eih_system.cpp ../util/logging.cpp)
set_target_properties(test-pr2-eih-system PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-pr2-eih-system PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-pr2-eih-system ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES})

add_executable(test-particles-sd  tests/test_particles_sd.cpp system/pr2_eih_system.cpp ../util/logging.cpp)
set_target_properties(test-particles-sd PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-particles-sd PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-particles-sd ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES})

add_executable(test-entropy  tests/test_entropy.cpp system/pr2_eih_system.cpp ../util/logging.cpp)
set_target_properties(test-entropy PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-entropy PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-entropy ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES})

add_executable(test-cost-teleop  tests/test_cost_teleop.cpp system/pr2_eih_system.cpp ../util/logging.cpp)
set_target_properties(test-cost-teleop PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-cost-teleop PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-cost-teleop ${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES})