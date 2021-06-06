cmake_minimum_required(VERSION 3.12)

add_test(
  NAME test_help_quetzal_egg_1
COMMAND ${PROJECT_BINARY_DIR}/src/qegg1 --help)

add_test(
  NAME test_version_quetzal_egg_1
COMMAND ${PROJECT_BINARY_DIR}/src/qegg1 --version)

add_test(
  NAME test_quetzal_egg_1
COMMAND ${PROJECT_BINARY_DIR}/src/qegg1 ${CMAKE_CURRENT_SOURCE_DIR}/data/quetzal_egg_1.config ${CMAKE_CURRENT_SOURCE_DIR}/data/suitability.tif ${CMAKE_CURRENT_SOURCE_DIR}/data/sample.csv)
