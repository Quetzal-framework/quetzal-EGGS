include_directories(headeronly_src)
include_directories(quetzal/include)

add_executable(quetzal-EGG1 qegg1.cpp)
target_include_directories(quetzal-EGG1 PUBLIC ${GDAL_INCLUDE_DIR})
target_link_libraries(quetzal-EGG1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES} stdc++fs)

install(
  TARGETS
    quetzal-EGG1
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
