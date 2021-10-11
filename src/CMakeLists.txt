include_directories(headeronly_src)
include_directories(quetzal/include)

add_executable(EGG1 EGG1.cpp)
target_include_directories(EGG1 PUBLIC ${GDAL_INCLUDE_DIR})
target_link_libraries(EGG1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES} stdc++fs)

add_executable(EGG2 EGG2.cpp)
target_include_directories(EGG2 PUBLIC ${GDAL_INCLUDE_DIR})
target_link_libraries(EGG2 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES} stdc++fs)

install(
  TARGETS
    EGG1
    EGG2
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    COMPONENT bin
  )
