cmake_minimum_required(VERSION 3.12)

add_test(
  NAME test_help_quetzal_EGG1
COMMAND ${PROJECT_BINARY_DIR}/src/EGG1 --help)

add_test(
  NAME test_version_quetzal_EGG1
COMMAND ${PROJECT_BINARY_DIR}/src/EGG1 --version)

add_test(
  NAME test_quetzal_EGG1
COMMAND ${PROJECT_BINARY_DIR}/src/EGG1 ${CMAKE_CURRENT_SOURCE_DIR}/data/quetzal_EGG1.config ${CMAKE_CURRENT_SOURCE_DIR}/data/suitability.tif ${CMAKE_CURRENT_SOURCE_DIR}/data/sample.csv)

add_test(
  NAME test_help_quetzal_EGG2
COMMAND ${PROJECT_BINARY_DIR}/src/EGG2 --help)

add_test(
  NAME test_version_quetzal_EGG2
COMMAND ${PROJECT_BINARY_DIR}/src/EGG1 --version)

add_test(
  NAME test_quetzal_EGG2
COMMAND ${PROJECT_BINARY_DIR}/src/EGG2 ${CMAKE_CURRENT_SOURCE_DIR}/data/quetzal_EGG2.config ${CMAKE_CURRENT_SOURCE_DIR}/data/suitability.tif ${CMAKE_CURRENT_SOURCE_DIR}/data/sample.csv)


install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_PREFIX})
