cmake_minimum_required (VERSION 2.6.0)
project (planar)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

#---- Now check the compiler
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall -Wno-sign-compare")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

find_package(PythonLibs 2.7 REQUIRED)
find_package(Boost COMPONENTS iostreams python thread filesystem system numpy)
find_package(Eigen)

set(FIGTREE_INCLUDE_DIR "/home/gkahn/source/figtree-0.9.3/include")
set(FIGTREE_LIBRARIES_DIR "/home/gkahn/source/figtree-0.9.3/lib")
set(FIGTREE_LIBRARIES "figtree;ann_figtree_version")

add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")

add_definitions(${CMAKE_CXX_FLAGS_DEBUG} "-g")
add_definitions(${CMAKE_C_FLAGS_DEBUG} "-g")

set(FADBAD_INCLUDE_DIR "/home/gkahn/source/FADBAD++")
set(MY_LIBRARIES "dl;rt;python2.7")

include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS}
					${Eigen_INCLUDE_DIRS} ${FADBAD_INCLUDE_DIR} ${FIGTREE_INCLUDE_DIR})

link_directories(${Boost_LIBRARY_DIRS} ${PYTHON_LIBRARIES} ${FIGTREE_LIBRARIES_DIR})


add_executable(planar planar.cpp planarMPC.c
                      src/planar-system.cpp src/planar-utils.cpp src/geometry2d.cpp src/gmm.cpp
                      fadbad/fadbad-planar-system.cpp fadbad/fadbad-geometry2d.cpp
                      ../util/logging.cpp)
target_link_libraries(planar ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MY_LIBRARIES}) 


add_executable(test-planar-system tests/test-planar-system.cpp
                                  src/planar-system.cpp src/geometry2d.cpp src/planar-utils.cpp
                                  fadbad/fadbad-planar-system.cpp fadbad/fadbad-geometry2d.cpp
                                  ../util/logging.cpp)
target_link_libraries(test-planar-system ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MY_LIBRARIES})

add_executable(test-geometry2d    tests/test-geometry2d.cpp
                                  src/geometry2d.cpp src/planar-utils.cpp
                                  ../util/logging.cpp)
target_link_libraries(test-geometry2d ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MY_LIBRARIES})

add_executable(test-figtree       tests/test-figtree.cpp src/planar-utils.cpp ../util/logging.cpp)
target_link_libraries(test-figtree ${FIGTREE_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MY_LIBRARIES})