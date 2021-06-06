include_directories(headeronly_src)
include_directories(quetzal/include)

add_executable(qegg1 qegg1.cpp)
target_include_directories(qegg1 PUBLIC ${GDAL_INCLUDE_DIR})
target_link_libraries(qegg1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES} stdc++fs)


install(
  TARGETS
    qegg1
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
