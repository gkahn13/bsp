cmake_minimum_required (VERSION 2.6.0)
project (pr2_bsp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall -Wno-sign-compare -Wno-enum-compare -std=c++11")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(OpenRAVE REQUIRED)
find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams thread filesystem system program_options graph)
find_package(Eigen)
find_package(CUDA)


set(PCL_INCLUDE_DIR "/usr/local/include/pcl-1.7")
set(PCL_LIBRARIES_DIR "/usr/local/lib")
#set(PCL_LIBRARIES "pcl_gpu_kinfu_large_scale;pcl_apps;pcl_common;pcl_filters;pcl_io;pcl_kdtree;pcl_keypoints;pcl_octree;pcl_outofcore;pcl_visualization")

#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_2d")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_apps")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_common")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_cuda_features")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_cuda_sample_consensus")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_cuda_segmentation")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_features")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_filters")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_geometry")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_containers")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_features")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_kinfu_large_scale")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_kinfu")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_octree")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_segmentation")
set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_gpu_utils")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_io_ply")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_io")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_kdtree")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_keypoints")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_ml")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_octree")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_outofcore")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_people")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_recognition")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_registration")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_sample_consensus")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_search")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_segmentation")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_stereo")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_surface")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_tracking")
#set(PCL_LIBRARIES "${PCL_LIBRARIES};pcl_visualization")


set(FIGTREE_INCLUDE_DIR "/home/gkahn/source/figtree-0.9.3/include")
set(FIGTREE_LIBRARIES_DIR "/home/gkahn/source/figtree-0.9.3/lib")
set(FIGTREE_LIBRARIES "figtree;ann_figtree_version")


include_directories(${OpenRAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS}
					${Eigen_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${FIGTREE_INCLUDE_DIR} ${PCL_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})

link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${FIGTREE_LIBRARIES_DIR} ${PCL_LIBRARIES_DIR})

add_executable(pr2  pr2.cpp pr2MPC.c
					geometry/geometry2d.cpp geometry/geometry3d.cpp
					system/pr2-sim.cpp system/pr2-system.cpp system/voxel-grid.cpp
					utils/pr2-utils.cpp utils/rave-utils.cpp
                    ../util/logging.cpp)
set_target_properties(pr2 PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(pr2 PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(pr2 ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES} ${PCL_LIBRARIES})


add_executable(test-pr2-system  tests/test-pr2-system.cpp
								geometry/geometry2d.cpp geometry/geometry3d.cpp
								system/pr2-sim.cpp system/pr2-system.cpp system/voxel-grid.cpp
								utils/pr2-utils.cpp utils/rave-utils.cpp
		                        ../util/logging.cpp)
set_target_properties(test-pr2-system PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-pr2-system PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-pr2-system ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES} ${PCL_LIBRARIES})


add_executable(test-fov tests/test-fov.cpp
						geometry/geometry2d.cpp geometry/geometry3d.cpp
						system/pr2-sim.cpp system/pr2-system.cpp
						utils/pr2-utils.cpp utils/rave-utils.cpp
                        ../util/logging.cpp)
set_target_properties(test-fov PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-fov PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-fov ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES})