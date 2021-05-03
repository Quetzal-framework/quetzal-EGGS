cmake_minimum_required(VERSION 3.12)

project(quetzal-EGGS-core)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package(SQLite3 REQUIRED)
IF (SQLite3_FOUND)
    INCLUDE_DIRECTORIES(${SQLite3_INCLUDE_DIRS})
ENDIF()

find_package(GDAL REQUIRED)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

include_directories(headeronly_src)
include_directories(quetzal/include)

add_executable(qegg1 qegg1.cpp)
target_link_libraries(qegg1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

install(
  TARGETS
    qegg1
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
