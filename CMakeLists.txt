cmake_minimum_required (VERSION 2.8) 
project (ODE) 

#add_subdirectory () 
#add_subdirectory () 
#include_directories () 
#link_directories () 


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
message(STATUS ${CMAKE_MODULE_PATH})
#find_package(Eigen3 COMPONENTS core required)
#if (EIGEN3_FOUND)
#  include_directories("${EIGEN3_INCLUDE_DIRS}")
#  message(STATUS "${EIGEN3_INCLUDE_DIRS}")
#endif (EIGEN3_FOUND)

set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/Eigen/")
include_directories(${EIGEN3_INCLUDE_DIR})

add_library( ode45 ode45.cpp)




#target_link_libraries () 