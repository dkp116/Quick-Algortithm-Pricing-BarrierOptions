cmake_minimum_required(VERSION 3.10)

project(Refactor_Build_Test LANGUAGES CXX)

# -------------------
# C++ standard
# -------------------
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------
# Default build type
# -------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# -------------------
# Compiler warnings
# -------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(DEFAULT_WARNING_FLAGS -Wall -Wextra -Wpedantic)
endif()

# -------------------
# Library: shared code
# -------------------
set(LIB_SOURCES
    src/Options/Barrier.cpp
    src/Dynamics/BlackScholesDynamics.cpp
    src/Dynamics/MertonJumpDynamics.cpp
    src/PricingAlgorithm/UniformSampleEstimate.cpp
    src/PricingAlgorithm/TaylorApproximation.cpp
    src/PricingAlgorithm/EstimateGI.cpp
)

add_library(core_lib STATIC ${LIB_SOURCES})

target_include_directories(core_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(core_lib PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

# -------------------
# Main executable
# -------------------
add_executable(Main src/main.cpp)
target_link_libraries(Main PRIVATE core_lib)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(Main PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

# -------------------
# Testing with Catch2
# -------------------
include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.2
)

FetchContent_MakeAvailable(Catch2)

file(GLOB_RECURSE TEST_SOURCES tests/*.cpp)

add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests PRIVATE core_lib Catch2::Catch2WithMain)

if(DEFAULT_WARNING_FLAGS)
    target_compile_options(tests PRIVATE ${DEFAULT_WARNING_FLAGS})
endif()

enable_testing()
add_test(NAME MonteCarloTests COMMAND $<TARGET_FILE:tests>)
