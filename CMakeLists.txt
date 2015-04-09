cmake_minimum_required(VERSION 3.1)
project(SeqIg)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-register")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


link_directories("lib/")
include_directories("include/")
set(SOURCE_FILES
    src/AlignAntibody.cpp
    src/AlignAntibody.h
    src/AntibodyJunction.cpp
    src/AntibodyJunction.h
    src/DatabaseHandler.cpp
    src/DatabaseHandler.h
    src/OutputHandler.cpp
    src/OutputHandler.h
    src/PropertiesHandler.cpp
    src/PropertiesHandler.h
    src/SeqIg.cpp
    src/SeqIg.h
    src/StructDefs.h
    src/Utility.h)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS filesystem system)

#if(Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(SeqIg ${SOURCE_FILES})
target_link_libraries(SeqIg ${Boost_LIBRARIES})

