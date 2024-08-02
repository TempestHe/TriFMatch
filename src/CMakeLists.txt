set(FILE_SRC
        preprocessor.cpp
        query_plan_generator.cpp
        trie_encoder.cpp
        successor_equivalent_cache.cpp
        subgraph_enumeration.cpp
        containment_filter.cpp
        candidates_recorder.cpp
        # subtree_reducer.cpp
)

# Add source files for the executable
# add_executable(debug.o debug.cpp ${FILE_SRC})
# # add_executable(get_matches.o get_match_result.cpp ${FILE_SRC})

# # Link against the libraries
# # target_link_libraries(debug.o configuration)
# target_link_libraries(debug.o graph)
# target_link_libraries(debug.o nd)
# target_link_libraries(debug.o utility)
# # target_link_libraries(debug.o model)
# target_link_libraries(debug.o rt)
# target_link_libraries(debug.o)

# # target_link_libraries(get_matches.o graph)
# # target_link_libraries(get_matches.o nd)
# # target_link_libraries(get_matches.o utility)
# # target_link_libraries(get_matches.o rt)

add_executable(main.o main.cpp ${FILE_SRC})
target_link_libraries(main.o graph)
target_link_libraries(main.o nd)
target_link_libraries(main.o utility)
target_link_libraries(main.o rt)
# target_link_libraries(main.o)

# Set the output directory for built binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
