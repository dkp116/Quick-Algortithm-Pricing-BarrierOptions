cmake_minimum_required(VERSION 3.10)

project(Refactor_Build_Test)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Default build type: Debug ---
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# --- Ensure debug symbols and no optimization for Debug ---
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# --- Source files ---
set(SOURCES
    main.cpp
    Barrier.cpp
    BlackScholesDynamics.cpp
    MertonJumpDynamics.cpp
    UniformSampleEstimate.cpp
    TaylorApproximation.cpp
    EstimateGI.cpp
)

# --- Executable ---
add_executable(Main ${SOURCES})

# --- Include directories ---
target_include_directories(Main PRIVATE ${CMAKE_SOURCE_DIR})

# --- Compiler warnings ---
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(Main PRIVATE -Wall -Wextra -Wpedantic)
endif()
