cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)

rosbuild_add_library(source_est src/types.cpp src/particles.cpp src/planner.cpp
								src/estimate.cpp src/measurement.cpp)
rosbuild_add_executable(static_source src/static_source_node.cpp)
rosbuild_add_executable(mi_planner src/mi_plan_node.cpp)
rosbuild_add_executable(range_sim src/range_sim_node.cpp)
rosbuild_add_executable(bagread src/bagread.cpp)

target_link_libraries(static_source source_est)
target_link_libraries(mi_planner source_est)
target_link_libraries(range_sim source_est)
target_link_libraries(bagread source_est)

rosbuild_add_gtest(test_info src/test_main.cpp src/information_test.cpp)
rosbuild_add_gtest(test_plan src/test_main.cpp src/planner_test.cpp)
rosbuild_add_gtest(test_particles src/test_main.cpp src/particles_test.cpp)

target_link_libraries(test_info source_est)
target_link_libraries(test_plan source_est)
target_link_libraries(test_particles source_est)

rosbuild_add_executable(static_source_bag src/static_source_bag.cpp)
target_link_libraries(static_source_bag source_est)
