# The test suite is the project's name with _test appended
set (TEST_NAME ${PROJECT_NAME}_TEST)

# Include googletest
include_directories(../modules/googletest/googletest/include)

# The test suite uses GOOGLE_TEST
find_library(GOOGLE_TEST gtest PATHS ${PROJECT_SOURCE_DIR}/lib)

# Test sources 
set (TEST_SOURCES
  dtest_single_functions.cpp
  dtest_uni_multi_modal_functions.cpp
  dtest_hybrid_functions.cpp
  dtest_composition_functions.cpp

  cec17_test_func.cpp
  test_utils.cpp
  )

# Build the test executable
add_executable(${TEST_NAME} gtest_main.cpp ${TEST_SOURCES})
set(CMAKE_BUILD_TYPE Debug)

# Add the neccessary libraries
target_link_libraries(${TEST_NAME} ${CMAKE_THREAD_LIBS_INIT} ${GOOGLE_TEST} gcov)

# Add flags for debugging and for code coverage tests
set_target_properties(${TEST_NAME} PROPERTIES COMPILE_FLAGS
  "-g -Wall -Wextra -fprofile-arcs -ftest-coverage")
set_target_properties(${TEST_NAME} PROPERTIES LINKER_FLAGS
  "-fprofile-arcs -ftest-coverage -lgcov")

add_custom_target(gtest
  COMMAND ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Running tests.")

# Use CodeCoverage to analyze the results
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
include(CodeCoverage)
set(LCOV_REMOVE_EXTRA "'${PROJECT_SOURCE_DIR}/modules/*'")
setup_target_for_coverage(coverage ${TEST_NAME} coverage)
