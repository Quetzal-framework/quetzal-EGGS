cmake_minimum_required(VERSION 3.12)

find_package(SQLite3 QUIET)

if(SQLite3_FOUND)
  message(STATUS "Found SQLite3: version ${SQLite3_VERSION}")
  add_library(sqlite3_external INTERFACE)  # dummy
else()
  message(STATUS "SQLite3 could not be located, downloading and building version 3.31.1 instead.")

  # Add CMake project file
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp)
  file(WRITE ${CMAKE_BINARY_DIR}/temp/CMakeLists.txt
      "cmake_minimum_required(VERSION 3.12)\n"
      "project(sqlite3)\n"
      "include_directories(${CMAKE_BINARY_DIR}/subprojects/Source/sqlite3_external)\n"
      "add_library(sqlite3 ${CMAKE_BINARY_DIR}/subprojects/Source/sqlite3_external/sqlite3.c)\n"
      "install(TARGETS sqlite3 DESTINATION ${STAGED_INSTALL_PREFIX}/sqlite3/lib)\n"
      "install(FILES sqlite3.h DESTINATION ${STAGED_INSTALL_PREFIX}/sqlite3/include)\n")

  # Download, configure, build and install.
  include(ExternalProject)
  ExternalProject_Add(sqlite3_external
      URL               https://www.sqlite.org/2020/sqlite-autoconf-3310100.tar.gz
      URL_HASH          SHA1=0c30f5b22152a8166aa3bebb0f4bc1f3e9cc508b
      DOWNLOAD_NO_PROGRESS  1
      UPDATE_COMMAND    ${CMAKE_COMMAND} -E copy
                        ${CMAKE_BINARY_DIR}/temp/CMakeLists.txt
                        <SOURCE_DIR>/CMakeLists.txt
      BUILD_COMMAND     ${CMAKE_COMMAND} --build .
      LOG_BUILD         1
      BUILD_IN_SOURCE   1
    )

    set(
      SQLite3_ROOT ${STAGED_INSTALL_PREFIX}/sqlite3
      CACHE PATH "Path to internally built SQLite3 installation root"
      FORCE
      )
    set(
      SQLite3_INCLUDE_DIRS ${STAGED_INSTALL_PREFIX}/sqlite3/include
      CACHE PATH "Path to internally built SQLite3 include directories"
      FORCE
      )
    set(
      SQLite3_LIBRARIES ${STAGED_INSTALL_PREFIX}/sqlite3/lib
      CACHE PATH "Path to internally built SQL3te2 library directories"
      FORCE
      )

endif()
