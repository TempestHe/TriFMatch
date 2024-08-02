#pragma once
#include "../configuration/config.h"
#include "../utility/relation/catalog.h"
#include "../graph/graph.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>

class TrieRelation{
public:
    uint32_t size_;
    unordered_map<Vertex, pair<Vertex*, uint32_t>> offset_map; // id -> start_offset, size
    Vertex* children_;
    Vertex src_, dst_;

    TrieRelation(catalog* storage, Vertex src, Vertex dst);

    void remove_candidate(Vertex src_u, Vertex src_v, Graph* data_graph);

    ~TrieRelation();
};

class TrieEncoder{
public:
    catalog* storage_;
    vector<Vertex> order_;
    vector<Vertex> order_index_;
    Graph* query_graph_;
    Graph* data_graph_;
    vector<vector<TrieRelation*>> candidate_edges;
    vector<unordered_map<Vertex, pair<Vertex*, uint32_t>>> candidates_set_depth;
    vector<Vertex> first_successor_depth_;
    TrieEncoder(catalog* storage, vector<Vertex>& order, vector<Vertex>& order_index, Graph* query_graph, Graph* data_graph, vector<vector<uint32_t>>& successor_depth);
    
    void get_candidates(uint32_t u_depth, vector<Vertex>& result);
    bool validate_neighbor_containment(uint32_t cur_depth, Vertex root1, Vertex root2);
    unordered_map<Vertex, pair<Vertex*, uint32_t>>* get_candidates_map(uint32_t u_depth);
    Vertex* get_edge_candidate(uint32_t src_depth, uint32_t dst_depth, Vertex src_v, uint32_t& count);
    void get_joint_candidate(uint32_t suc_depth, uint32_t dst_depth, Vertex src_v, vector<Vertex>& result);

    ~TrieEncoder();
};
