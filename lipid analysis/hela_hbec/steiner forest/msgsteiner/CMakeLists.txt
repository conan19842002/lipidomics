cmake_minimum_required(VERSION 2.8)
set(Boost_INCLUDE_DIR /boost_1_62_0)
set(Boost_LIBRARY_DIR /boost_1_62_0/libs)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_executable(main main.cpp)
target_link_libraries( main ${Boost_LIBRARIES} )