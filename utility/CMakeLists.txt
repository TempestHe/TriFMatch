set(UTILS_SRC utils.cpp
QFilter.h
relation/catalog.cpp
primitive/nlf_filter.cpp
primitive/projection.cpp
primitive/scan.cpp
primitive/search.cpp
primitive/semi_join.cpp
han/intersection_algos.hpp
han/intersection_algos.cpp
han/utils/util.cpp
han/utils/util.hpp
graphoperations.cpp
computesetintersection.cpp
)

# Add source files for the library
add_library(utility STATIC ${UTILS_SRC})

# Set include directories for the target
target_include_directories(utility PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(nucleus_decomposition)