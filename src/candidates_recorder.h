#pragma once
#include "../configuration/config.h"
#include "../graph/graph.h"
#include "trie_encoder.h"
#include "../utility/computesetintersection.h"
#include <vector>


class CandidatesHistoryStack{
public:
    vector<Vertex*> block_content_;
    Vertex* current_block_ptr;
    uint32_t last_block_remaining_size_;
    uint32_t total_block_size_;
    uint32_t block_size_;
    uint32_t cur_depth_;
    vector<Vertex> successors_;

    // elements to be recorded
    vector<Vertex> history_root_vertices_;
    vector<bitset<MAX_QUERY_SIZE>> history_conflict_sources_;
    vector<bitset<MAX_QUERY_SIZE>> history_failing_set_;
    vector<vector<pair<Vertex*, uint32_t>>> history_candidates_;
    // vector<int> validation_marker_;
    vector<uint32_t> validate_exclusion_offsets_;
    vector<uint32_t> validate_completion_offset_;


    CandidatesHistoryStack(Graph* query_graph, Graph* data_graph, vector<Vertex>& successors, uint32_t cur_depth);
    void rebuild(Graph* query_graph, Graph* data_graph, vector<Vertex>& successors, uint32_t cur_depth);

    void clear();

    void append_history(Vertex v, TrieEncoder* encoder, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack);

    CandidatesHistoryStack(){};

    ~CandidatesHistoryStack();
};