cmake_minimum_required(VERSION 3.25)

project(Refactor_Build_Test LANGUAGES CXX)

# -------------------
# C++ standard
# -------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -------------------
# Default build type (only if single-config generator)
# -------------------
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# -------------------
# Compiler warnings
# -------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(DEFAULT_WARNING_FLAGS -Wall -Wextra -Wpedantic)
elseif(MSVC)
    set(DEFAULT_WARNING_FLAGS /W4)
endif()

# -------------------
# Core library
# -------------------
set(LIB_SOURCES
    src/Options/Barrier.cpp
    src/Dynamics/BlackScholesDynamics.cpp
    src/Dynamics/MertonJumpDynamics.cpp
    src/PricingAlgorithm/UniformSampleEstimate.cpp
    src/PricingAlgorithm/TaylorApproximation.cpp
    src/PricingAlgorithm/EstimateGI.cpp
    src/PricingAlgorithm/StandardMonteCarlo.cpp
)

add_library(core_lib STATIC ${LIB_SOURCES})

target_include_directories(core_lib
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(core_lib PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

# -------------------
# Precompiled headers (huge compile speed win)
# -------------------
target_precompile_headers(core_lib PRIVATE
    <vector>
    <memory>
    <random>
    <cmath>
    <iostream>
)

# -------------------
# Unity build (massive speedup for template-heavy projects)
# -------------------
set_target_properties(core_lib PROPERTIES UNITY_BUILD ON)

# -------------------
# Main executable
# -------------------
add_executable(Main src/main.cpp)

target_link_libraries(Main PRIVATE core_lib)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(Main PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

# -------------------
# Catch2 (tests)
# -------------------
include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
)

FetchContent_MakeAvailable(Catch2)

# -------------------
# Test sources (explicit listing — no GLOB)
# -------------------
set(TEST_SOURCES
    tests/test_utils/TestDownAndOutPricing.cpp
    tests/test_first.cpp
    # add other test files explicitly here
)

add_executable(tests ${TEST_SOURCES})

target_include_directories(tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tests/test_utils
)

target_link_libraries(tests
    PRIVATE
        core_lib
        Catch2::Catch2WithMain
)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(tests PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

# Unity build for tests too
set_target_properties(tests PROPERTIES UNITY_BUILD ON)

# -------------------
# Enable CTest
# -------------------
enable_testing()
add_test(NAME MonteCarloTests COMMAND tests)