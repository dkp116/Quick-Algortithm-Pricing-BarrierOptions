cmake_minimum_required(VERSION 3.25)

project(barrier_pricing LANGUAGES CXX)

# -----------------------
# C++ Settings
# -----------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -----------------------
# Library
# -----------------------
add_library(core STATIC
    src/Options/Barrier.cpp
    src/Dynamics/BlackScholesDynamics.cpp
    src/Dynamics/MertonJumpDynamics.cpp
    src/PricingAlgorithm/UniformSampleEstimate.cpp
    src/PricingAlgorithm/TaylorApproximation.cpp
    src/PricingAlgorithm/EstimateGI.cpp
    src/PricingAlgorithm/StandardMonteCarlo.cpp
)

target_include_directories(core PUBLIC include)

# Faster builds
set_target_properties(core PROPERTIES UNITY_BUILD ON)

target_precompile_headers(core PRIVATE
    <vector>
    <memory>
    <random>
    <cmath>
)

# -----------------------
# Main executable
# -----------------------
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE core)

# -----------------------
# Tests (optional)
# -----------------------
option(BUILD_TESTS "Build unit tests" ON)

if(BUILD_TESTS)

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
)

FetchContent_MakeAvailable(Catch2)

add_executable(tests
    tests/test_utils/TestDownAndOutPricing.cpp
    tests/intergration_tests.cpp
)

target_link_libraries(tests PRIVATE core Catch2::Catch2WithMain)

enable_testing()
add_test(NAME unit_tests COMMAND tests)

endif()