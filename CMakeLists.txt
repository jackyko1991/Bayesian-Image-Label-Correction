cmake_minimum_required(VERSION 2.8)

project(BayesianFilter)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# source files
file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
file(GLOB HDR ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(BayesianFilter ${SRC} ${HDR})

target_link_libraries(BayesianFilter ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} )