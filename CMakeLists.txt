cmake_minimum_required( VERSION 2.6 )

project( remesher )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_SOURCE_DIR}/bin )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/bin )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin )

find_package( utilities )

include_directories( 
 ${UTILITIES_INCLUDE_DIRS}
 ${CMAKE_SOURCE_DIR}/include 
)

set( REMESHER_HEADERS
 include/mesh.h
 include/remesher.h
)

add_executable( remesher source/remesher_main.cpp ${REMESHER_HEADERS} )