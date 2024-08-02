#pragma once
#include <iostream>
#include "../graph/graph.h"
#include "../utility/utils.h"
#include "trie_encoder.h"

#define DEBUG 1


struct HashVertexPair{
    size_t operator() (const pair<Vertex, Vertex>& key) const{
        std::hash<int> hasher;
        size_t seed = 0;
        seed ^= hasher(key.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(key.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct HashPointerPair{
    size_t operator() (const pair<void*, void*>& key) const{
        std::hash<void*> hasher;
        size_t seed = 0;
        seed ^= hasher(key.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(key.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

class SubtreeReducer{
public:
    Graph* query_graph_;
    Graph* data_graph_;
    vector<vector<uint32_t>> successors_depth_;
    vector<vector<uint32_t>> predecessor_depth_;
    vector<vector<uint32_t>> none_succesors_depth_;
    vector<bitset<MAX_QUERY_SIZE>> ancestors_depth_;
    vector<bitset<MAX_QUERY_SIZE>> mask_;
    vector<Vertex> order_;
    uint32_t query_vertex_count_;
    TrieEncoder* encoder_;

    vector<vector<vector<vector<Vertex>>>> history_;
    vector<vector<uint32_t>> history_max_search_depth;
    vector<vector<Vertex>> history_root_vertices_;
    vector<vector<bitset<MAX_QUERY_SIZE>>> history_conflict_source_;

    vector<int> successors_with_same_labels_;

    vector<vector<uint32_t>> empty_record_offset_;
    vector<vector<uint32_t>> conflict_record_offset_;
    vector<unordered_map<Vertex, FailingSet*>> failing_set_recorder_;

    vector<pair<vector<uint32_t>, vector<uint32_t>>> unmatched_query_depth_dist_; // first->successors, second->none_successors
#if DEBUG==1
    vector<pair<vector<Vertex>, vector<Vertex>>> same_label_depth_; // first->successors, second->none_successors;
    vector<vector<pair<vector<int>, vector<int>>>> preprocessed_depth_intersection; //(the depth v_1 resides, current_depth)==>(first->matched_depth that do not connected u_c, second->unmatched...) 
#endif

    SubtreeReducer(Graph* query_graph_, Graph* data_graph_, vector<vector<uint32_t>>& successors_depth, vector<vector<uint32_t>>& predecessor_depth, vector<Vertex>& search_order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, TrieEncoder* encoder);

    void clear_record(uint32_t depth);

    bool validate_standard_containment(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent);
    bool validate_standard_containment_test(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent, Vertex& history_root_test);

    bool validate_extended_containment(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs);
    bool validate_extended_containment_debug(uint32_t depth, vector<vector<vector<Vertex>>>& candidates_stack, Vertex root, FailingSet*& result_fs, bool full_descendent, vector<Vertex>& embedding_depth);

    void mark_candidates_to_be_validated(uint32_t depth, Vertex root, FailingSet* fs, vector<vector<vector<Vertex>>>& candidates_stack, bitset<MAX_QUERY_SIZE>& conflict_source, bool full_descendent, uint32_t max_search_depth);
    
    bool validate_vector_equivalent(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2, bool& contain_r2);

    bool validate_vector_containment_with_conflict(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2);

    bool validate_vector_containment_relaxed(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2);
    bool validate_vector_containment_relaxed(const Vertex* vec1, uint32_t vec1_size, const Vertex* vec2, uint32_t vec2_size, Vertex r1, Vertex r2);
    bool validate_vector_containment_relaxed_full(vector<Vertex>& vec1, vector<Vertex>& vec2, Vertex r1, Vertex r2);
};
