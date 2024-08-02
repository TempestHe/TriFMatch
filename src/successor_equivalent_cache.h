#pragma once

#include "../configuration/config.h"
#include "../graph/graph.h"
#include "../utility/utils.h"


struct MatchInfo{
    uint32_t emb_count;
    bitset<MAX_QUERY_SIZE> empty_set_failure;
    vector<Vertex> conflict_vertex;
    vector<Vertex> conflict_depth_outer;
};




struct HashFunc_Vec{
    size_t operator() (const vector<Vertex>& key) const{
        std::hash<int> hasher;
        size_t seed = 0;
        for(auto v : key){
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct EquivalenceZone{
    uint32_t start_depth;
    uint32_t end_depth;
    bitset<MAX_QUERY_SIZE> equivalent_depth;
    vector<uint32_t> equivalent_depth_map;
    vector<uint32_t> equivalent_relative_offset;
    vector<vector<uint32_t>> candidates_to_validate;
    unordered_set<vector<Vertex>, HashFunc_Vec> recorder;
    // unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec> recorder_emb_count;
    // unordered_map<vector<Vertex>, vector<vector<Vertex>>, HashFunc_Vec> recorder_emb_result;
};


class SuccessorEquivalentCache{
public:
    float max_depth_ratio=0.8;
    int max_zone_size=100;
    uint32_t size_limit_=1000000;
    uint64_t emb_limit_;
    int compressed_;
    int compressed_neq_;
    int query_size_;
    
    vector<Vertex> order_;
    vector<Vertex> order_index_;
    Graph* query_graph_;
    Graph* data_graph_;
    vector<vector<Vertex>> successor_neighbors_in_depth_;
    vector<bitset<MAX_QUERY_SIZE>> ancestors_depth_;

    vector<EquivalenceZone*> start_index_;
    vector<EquivalenceZone*> end_index_;

    vector<EquivalenceZone*> equivalent_zones_;

    // vector<unordered_map<vector<Vertex>, MatchInfo, HashFunc_Vec>> recoder_emb_count_;
    // unordered_map<vector<Vertex>, vector<vector<Vertex>>, HashFunc_Vec> recoder_emb_result_;


    SuccessorEquivalentCache(vector<Vertex>& order, vector<Vertex>& order_index, Graph* query_graph, Graph* data_graph, vector<vector<Vertex>>& successor_neighbors_in_depth, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, uint64_t emb_limit);

    bool validate_cache_count(uint32_t cur_depth, Vertex* embedding);
    // bool validate_cache_result(uint32_t cur_depth, Vertex* embedding, uint64_t& result_count, Vertex* visited_query_vertices, vector<vector<Vertex>>& results);

    vector<uint32_t> result_count_recorder_;
    void start_insert_cache(uint32_t cur_depth, uint64_t& result_count);
    bool finish_insert_cache_count(uint32_t cur_depth, Vertex* embedding, uint64_t& result_count, Vertex* visited_query_vertices
#if PRINT_RESULT==1
    , vector<vector<Vertex>>& matches
#endif
    );
    
    
    // void finish_insert_cache_result(uint32_t cur_depth, Vertex* embedding, bitset<MAX_QUERY_SIZE>& conflict_source_depth, bitset<MAX_QUERY_SIZE>& empty_failure_depth, bitset<MAX_QUERY_SIZE>& failing_set_depth, uint64_t& result_count, Vertex* visited_query_vertices, vector<vector<Vertex>>& results);

    void clear(uint32_t cur_depth);

    ~SuccessorEquivalentCache();
};


struct SuccessorEquivalentPair{
    Vertex parent;
    Vertex u;
};