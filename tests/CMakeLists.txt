include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1
  GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(Catch2)

add_executable(tedder_tests
  test_domain.cpp)

target_link_libraries(tedder_tests PRIVATE tedder::tedder Catch2::Catch2WithMain)
tedder_configure_compiled_target(tedder_tests)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(tedder_tests DISCOVERY_MODE PRE_TEST)
