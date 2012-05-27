cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)
rosbuild_init()

SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -g")

file( GLOB KDL_SRCS src/[^.]*.cpp )
file( GLOB KDL_HPPS src/[^.]*.hpp [^.]*.inl)

file( GLOB UTIL_SRCS src/utilities/[^.]*.cpp src/utilities/[^.]*.cxx)
file( GLOB UTIL_HPPS src/utilities/[^.]*.h src/utilities/[^.]*.hpp)

rosbuild_add_library(orocos-kdl ${KDL_SRCS} ${UTIL_SRCS})

#Creation of python bindings using sip
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})
include(FindSIPandPython)
add_subdirectory(src/bindings/python)

