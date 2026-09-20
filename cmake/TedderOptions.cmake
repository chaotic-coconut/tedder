include_guard(GLOBAL)

option(TEDDER_ENABLE_CCACHE "Use ccache when it is available" ON)
option(TEDDER_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)

if(TEDDER_ENABLE_CCACHE)
  find_program(TEDDER_CCACHE_PROGRAM ccache)
  if(TEDDER_CCACHE_PROGRAM AND NOT CMAKE_CXX_COMPILER_LAUNCHER)
    execute_process(
      COMMAND "${TEDDER_CCACHE_PROGRAM}" --get-config cache_dir
      OUTPUT_VARIABLE TEDDER_CCACHE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(IS_DIRECTORY "${TEDDER_CCACHE_DIR}")
      set(TEDDER_CCACHE_TEST_DIR "${TEDDER_CCACHE_DIR}")
    else()
      get_filename_component(TEDDER_CCACHE_PARENT
        "${TEDDER_CCACHE_DIR}" DIRECTORY)
      set(TEDDER_CCACHE_TEST_DIR "${TEDDER_CCACHE_PARENT}")
    endif()
    set(TEDDER_CCACHE_TEST_FILE
      "${TEDDER_CCACHE_TEST_DIR}/.tedder-ccache-write-test")
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E touch "${TEDDER_CCACHE_TEST_FILE}"
      RESULT_VARIABLE TEDDER_CCACHE_WRITE_RESULT
      OUTPUT_QUIET ERROR_QUIET)
    if(TEDDER_CCACHE_WRITE_RESULT EQUAL 0)
      file(REMOVE "${TEDDER_CCACHE_TEST_FILE}")
      set(TEDDER_CCACHE_WRITABLE TRUE)
    endif()
    if(TEDDER_CCACHE_WRITABLE)
      set(CMAKE_CXX_COMPILER_LAUNCHER "${TEDDER_CCACHE_PROGRAM}"
          CACHE FILEPATH "Compiler launcher used for tedder builds")
    else()
      message(STATUS "ccache found, but its cache directory is not writable; disabled")
    endif()
  endif()
endif()

function(tedder_configure_compiled_target target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "Cannot configure missing target: ${target}")
  endif()

  target_compile_options("${target}" PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wshadow;-Wconversion;-Wextra-semi>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>)

  if(TEDDER_ENABLE_SANITIZERS)
    target_compile_options("${target}" PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang>:-fsanitize=address,undefined;-fno-omit-frame-pointer;-fno-sanitize-recover=undefined>)
    target_link_options("${target}" PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang>:-fsanitize=address,undefined;-fno-omit-frame-pointer;-fno-sanitize-recover=undefined>)
  endif()
endfunction()
