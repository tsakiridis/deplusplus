cmake_minimum_required (VERSION 3.3)
set (PROJECT_NAME DEPLUSPLUS)
project (${PROJECT_NAME})

# Set options for the build
option (BUILD_STATIC "Build static version" OFF)
option (BUILD_TESTS  "Build unit tests" OFF)
option (BUILD_DOC    "Build documentation" OFF)

# Export compile commands for YCM
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set binary directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# Set CXX standard and flags
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-Wall -pedantic")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
                            -fno-omit-frame-pointer -DDEBUG")

# Set CXX flags for static build
if (BUILD_STATIC)
  set (CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -static")
  set (CMAKE_EXE_LINKER_FLAGS
     "${CMAKE_EXE_LINKER_FLAGS} \
     -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()

# Add project dependencies
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

# Project Include directories
include_directories (include/
                     modules/ThreadPool/
                    )

# Build examples
add_executable(example_1 examples/1_shade_basic.cpp src/algorithm/shade.cpp)
target_link_libraries(example_1 Threads::Threads)

add_executable(example_2 examples/2_degl_basic.cpp src/algorithm/degl.cpp)
target_link_libraries(example_2 Threads::Threads)

# The flag CEC_MAX_EVALUATIONS sets an upper limit on the number of function
# evaluations for the algorithms
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCEC_MAX_EVALUATIONS")

add_executable(example_3 examples/3_shade_cec.cpp src/algorithm/shade.cpp)
target_link_libraries(example_3 Threads::Threads)

add_executable(example_4 examples/4_shade_cec.cpp src/algorithm/shade.cpp)
target_link_libraries(example_4 Threads::Threads)

# Build tests
if (BUILD_TESTS)
    add_subdirectory(test)
endif()

# Build documentation
if (BUILD_DOC)
    add_subdirectory(documentation)
endif()
