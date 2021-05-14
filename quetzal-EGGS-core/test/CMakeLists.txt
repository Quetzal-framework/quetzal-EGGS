cmake_minimum_required(VERSION 3.12)

add_test(
  NAME test_quetzal_egg_1
  COMMAND ${CMAKE_INSTALL_PREFIX}/bin/qegg1 ${CMAKE_CURRENT_SOURCE_DIR}/data/quetzal_egg_1.config ${CMAKE_CURRENT_SOURCE_DIR}/data/suitability.tif)
