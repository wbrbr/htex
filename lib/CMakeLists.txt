cmake_minimum_required(VERSION 3.8.0)
project(Htex)

option(HTEX_BUILD_STATIC_LIBS "Enable building Htex static libraries" ON)
option(HTEX_BUILD_SHARED_LIBS "Enable building Htex shared libraries" ON)
option(PRMAN_15_COMPATIBLE_HTEX "Enable PRMan 15 compatibility" OFF)

if (DEFINED ENV{CXXFLAGS_STD})
    if ($ENV{CXXFLAGS_STD} STREQUAL "c++17")
        set(CMAKE_CXX_STANDARD 17)
    elseif ($ENV{CXXFLAGS_STD} STREQUAL "c++14")
        set(CMAKE_CXX_STANDARD 14)
    elseif ($ENV{CXXFLAGS_STD} STREQUAL "c++11")
        set(CMAKE_CXX_STANDARD 11)
    else ()
        set(CMAKE_CXX_STANDARD 98)
    endif ()
else ()
    set(CMAKE_CXX_STANDARD 98)
endif ()
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_FIND_PACKAGE_RESOLVE_SYMLINKS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

include(GNUInstallDirs)

include(CTest)
enable_testing()

# Setup platform-specific threading flags.
find_package(Threads REQUIRED)

# Use pkg-config to create a PkgConfig::Ptex_ZLIB imported target
find_package(PkgConfig REQUIRED)
pkg_checK_modules(Ptex_ZLIB REQUIRED zlib IMPORTED_TARGET)


if (NOT DEFINED HTEX_SHA)
    # Query git for current commit ID
    execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE HTEX_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

if (NOT DEFINED HTEX_VER)
    # Get the version string from a "version" file or from git
    if (EXISTS "${PROJECT_SOURCE_DIR}/version")
        file(STRINGS "${PROJECT_SOURCE_DIR}/version" HTEX_VER)
    else ()
        execute_process(
            COMMAND git describe --match=v[0-9]*.[0-9]* --first-parent HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE HTEX_VER
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif ()
endif ()
# Provide a fallback version when the Git repository is unavailable.
if (NOT HTEX_VER)
    set(HTEX_VER "v2.4.X")
endif()

# Transform HTEX_VER into HTEX_MAJOR_VERSION and HTEX_MINOR_VERSION
string(REPLACE "v" "" HTEX_VER_STRIPPED ${HTEX_VER})  # strip leading "v"
string(REPLACE "." ";" HTEX_VER_LIST ${HTEX_VER_STRIPPED})
list(LENGTH HTEX_VER_LIST HTEX_VER_LENGTH)
if (${HTEX_VER_LENGTH} LESS 2)
    message(FATAL_ERROR "Could not determine the Ptex library version")
endif ()

# The version variables are used to generate PtexVersion.h
list(GET HTEX_VER_LIST 0 HTEX_MAJOR_VERSION)
list(GET HTEX_VER_LIST 1 HTEX_MINOR_VERSION)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    # Detect the build type from the $FLAVOR environment variable
    # Default to optimized Release builds when unspecified.
    if ("$ENV{FLAVOR}" MATCHES "debug")
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "type of build" FORCE)
    elseif ("$ENV{FLAVOR}" MATCHES "profile")
        set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "type of build" FORCE)
    else ()
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build" FORCE)
    endif ()
endif ()

if (NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif ()

if (MSVC)
    add_definitions(
        # forcing value to bool 'true' or 'false' (performance warning)
        /wd4800
        # conversion from 'double' to 'float', possible loss of data
        /wd4244
        /D_CRT_SECURE_NO_WARNINGS
    )
endif()

# Set the default Ptex library to use for tools and tests.
if (HTEX_BUILD_STATIC_LIBS)
    set(HTEX_LIBRARY Htex_static)
else()
    set(HTEX_LIBRARY Htex_dynamic)
endif()

if (${PRMAN_15_COMPATIBLE_PTEX})
    add_definitions(-DHTEX_NO_LARGE_METADATA_BLOCKS)
endif ()

include_directories(src/htex)

add_subdirectory(src/htex)
#add_subdirectory(src/utils)
#add_subdirectory(src/tests)
#add_subdirectory(src/doc)
add_subdirectory(src/build)
