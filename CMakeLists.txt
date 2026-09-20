cmake_minimum_required(VERSION 3.24)
project(tedder LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(TedderOptions)

add_library(tedder INTERFACE)
add_library(tedder::tedder ALIAS tedder)
target_include_directories(tedder INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(tedder INTERFACE cxx_std_23)

option(TEDDER_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})

if(TEDDER_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
