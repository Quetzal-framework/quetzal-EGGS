cmake_minimum_required(VERSION 3.12)

# project name and language
project(quetzal-EGGS LANGUAGES CXX)

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags, Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C++ flags, Release configuration: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "C++ flags, Release configuration with Debug info: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "C++ flags, minimal Release configuration: ${CMAKE_CXX_FLAGS_MINSIZEREL}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.61 REQUIRED COMPONENTS program_options)

find_package(SQLite3 REQUIRED)
if (SQLite3_FOUND)
    include_directories(${SQLite3_INCLUDE_DIRS})
endif()

find_package(GDAL REQUIRED)
if (GDAL_FOUND)
    include_directories(${GDAL_INCLUDE_DIR})
endif()

add_subdirectory(src)

enable_testing()
add_subdirectory(test)
