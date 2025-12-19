cmake_minimum_required(VERSION 3.10)

project(Refactor_Build_Test LANGUAGES CXX)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Default build type: Debug ---
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# --- Ensure debug symbols and no optimization for Debug ---
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# --- Source files (now under src/) ---
set(SOURCES
    src/main.cpp
    src/Options/Barrier.cpp
    src/Dynamics/BlackScholesDynamics.cpp
    src/Dynamics/MertonJumpDynamics.cpp
    src/PricingAlgorithm/UniformSampleEstimate.cpp
    src/PricingAlgorithm/TaylorApproximation.cpp
    src/PricingAlgorithm/EstimateGI.cpp
)

# --- Executable ---
add_executable(Main ${SOURCES})

# --- Include directories (point to include/) ---
target_include_directories(Main PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# --- Compiler warnings ---
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(Main PRIVATE -Wall -Wextra -Wpedantic)
endif()
