cmake_minimum_required(VERSION 3.1)

project(epsilon-prepare)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(Boost_REALPATH ON)

find_package(Boost 1.59.0 REQUIRED regex)
find_package(GiNaC 1.7.1 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GINAC_INCLUDE_DIRS})

file(GLOB SOURCES src/*.cpp)

add_executable(epsilon-prepare ${SOURCES})
target_link_libraries(epsilon-prepare ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES} ${GINAC_LIBRARIES})

install (TARGETS epsilon-prepare DESTINATION bin)

