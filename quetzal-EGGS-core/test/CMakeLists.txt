cmake_minimum_required(VERSION 3.12)

add_test(test_quetzal_egg_1 ${CMAKE_INSTALL_PREFIX}/bin/qegg1 ${CMAKE_CURRENT_SOURCE_DIR}/data/suitability.tif)
