cmake_minimum_required (VERSION 2.6.0)
project (pr2_bsp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall -Wno-sign-compare -std=c++11")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(OpenRAVE REQUIRED)
find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams thread filesystem system program_options)
find_package(Eigen)

set(FIGTREE_INCLUDE_DIR "/home/gkahn/source/figtree-0.9.3/include")
set(FIGTREE_LIBRARIES_DIR "/home/gkahn/source/figtree-0.9.3/lib")
set(FIGTREE_LIBRARIES "figtree;ann_figtree_version")

include_directories(${OpenRAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS}
					${Eigen_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${FIGTREE_INCLUDE_DIR})

link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${FIGTREE_LIBRARIES_DIR})

add_executable(pr2  pr2.cpp pr2MPC.c
					geometry/geometry2d.cpp geometry/geometry3d.cpp
					system/pr2-sim.cpp system/pr2-system.cpp
					utils/pr2-utils.cpp utils/rave-utils.cpp
                    ../util/logging.cpp)
set_target_properties(pr2 PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(pr2 PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(pr2 ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES})


add_executable(test-pr2-system  tests/test-pr2-system.cpp
								geometry/geometry2d.cpp geometry/geometry3d.cpp
								system/pr2-sim.cpp system/pr2-system.cpp
								utils/pr2-utils.cpp utils/rave-utils.cpp
		                        ../util/logging.cpp)
set_target_properties(test-pr2-system PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-pr2-system PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-pr2-system ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES})


add_executable(test-fov tests/test-fov.cpp
						geometry/geometry2d.cpp geometry/geometry3d.cpp
						system/pr2-sim.cpp system/pr2-system.cpp
						utils/pr2-utils.cpp utils/rave-utils.cpp
                        ../util/logging.cpp)
set_target_properties(test-fov PROPERTIES COMPILE_FLAGS "${OPENRAVE_CXXFLAGS}")
set_target_properties(test-fov PROPERTIES LINK_FLAGS "${OPENRAVE_LINK_FLAGS}")
target_link_libraries(test-fov ${OpenRAVE_LIBRARIES} ${OpenRAVE_CORE_LIBRARIES} ${Boost_LIBRARIES} ${FIGTREE_LIBRARIES})