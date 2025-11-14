# Initial cache variables (applied before first configure)
# These are defaults that can be overridden by presets or -D flags

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require C++ standard")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable compiler extensions")

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Build options
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
set(ENABLE_SANITIZERS OFF CACHE BOOL "Enable sanitizers in Debug")