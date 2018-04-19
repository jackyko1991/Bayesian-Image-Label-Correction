cmake_minimum_required(VERSION 2.8)

project(BayesianFilter)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(BayesianFilter BayesianFilter.cxx )

target_link_libraries(BayesianFilter ${ITK_LIBRARIES})