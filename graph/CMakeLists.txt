# Add source files for the library
add_library(graph STATIC graph.cpp)

# Set include directories for the target
target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})