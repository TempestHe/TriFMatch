cmake_minimum_required(VERSION 3.5)

project(PPC-index
LANGUAGES CXX
)

# set libtorch
# add_definitions(-D _GLIBCXX_USE_CXX11_ABI=0)
set(CMAKE_CXX_STANDARD 17)
# find_package(Torch REQUIRED)
# message(STATUS "path=${TORCH_LIBRARIES}")
# message(STATUS "include=${TORCH_INCLUDE_DIR}")
# message(STATUS "flags=${TORCH_CXX_FLAGS}")

# set OGDF
# find_package(OGDF REQUIRED)
# if (NOT OGDF_FOUND)
#     message(FATAL_ERROR "OGDF NOT FOUND")
# endif (NOT OGDF_FOUND)

# Set C++ standard to C++14
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -std=c++11 -pthread -lrt -march=native -g -O3")


# Add subdirectories
# add_subdirectory(configuration)


add_subdirectory(graph)
# add_subdirectory(graph_decomposition)
# add_subdirectory(index)
# add_subdirectory(model)
# add_subdirectory(nd)
add_subdirectory(utility)
add_subdirectory(src)
