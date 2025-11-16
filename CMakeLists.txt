cmake_minimum_required(VERSION 3.31)
project(ASICPlacement)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)

# Loadinitial cache(sets defaults before first configure)
# if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/Cache.cmake)
# include(${CMAKE_SOURCE_DIR}/cmake/Cache.cmake)
# endif()

# fetch spdlog
find_package(spdlog NO_MODULE QUIET)

if(NOT spdlog_FOUND)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.14.1
    )
    FetchContent_MakeAvailable(spdlog)
endif()

find_package(Eigen3 REQUIRED NO_MODULE QUIET)

if(NOT Eigen3_FOUND)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

add_executable(QPEngine QPEngine.cpp)
target_link_libraries(QPEngine Eigen3::Eigen spdlog::spdlog)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # add_compile_definitions(DEBUG_PRINT DEBUG_BREAKPOINT)
    target_compile_definitions(QPEngine PRIVATE DEBUG_PRINT DEBUG_BREAKPOINT)

    # target_compile_definitions(QPEngine PRIVATE DEBUG_PRINT)
endif()