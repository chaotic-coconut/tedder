include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1
  GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(Catch2)

add_executable(tedder_tests
  test_domain.cpp)

target_link_libraries(tedder_tests PRIVATE tedder::tedder Catch2::Catch2WithMain)

target_compile_options(tedder_tests PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wshadow;-Wconversion;-Wextra-semi>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(tedder_tests)
