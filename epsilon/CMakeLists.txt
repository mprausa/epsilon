cmake_minimum_required(VERSION 3.1)
project(epsilon)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 11)

find_package(libFermat REQUIRED)
find_path(PSTREAMS_INCLUDE_DIR pstream.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${PSTREAMS_INCLUDE_DIR})
include_directories(${LIBFERMAT_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY})

file(GLOB SOURCES src/*.cpp)

set(FUNCTIONS_FER ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/functions_fer.h)
add_custom_command(
	OUTPUT ${FUNCTIONS_FER}
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fer/functions.fer
	COMMAND sed 's|^|\"|\;s|$$|\",|\;1s|^|static vector<string> _functions_fer = {\\n|\;$$s|,$$|\\n}\;|' ${CMAKE_CURRENT_SOURCE_DIR}/fer/functions.fer > ${FUNCTIONS_FER}
	COMMENT "generating ${FUNCTIONS_FER}"
)
add_custom_target(functions_fer DEPENDS ${FUNCTIONS_FER})

add_executable(epsilon ${SOURCES})
add_dependencies(epsilon functions_fer)
target_link_libraries(epsilon ${LIBFERMAT_LIBRARIES})

install (TARGETS epsilon DESTINATION bin)

