#pragma once
#include <queue>
#include "../graph/graph.h"
#include "../configuration/config.h"
#include "preprocessor.h"
// #include "encoder.h"
#include "trie_encoder.h"
#include "query_plan_generator.h"

#include "candidates_recorder.h"
#include "containment_filter.h"

#include "successor_equivalent_cache.h"
#include "subtree_reducer.h"

#if PRINT_MEM_INFO == 1
#include <thread>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <unistd.h>

#define VMRSS_LINE 22
#define PROCESS_ITEM 14

inline int GetCurrentPid();

inline float GetMemoryUsage(int pid);

void thread_get_mem_info(float& peak_memory, bool& stop);
#endif

enum LeafStateType{
    CD_STATE, FD_STATE, FAILINGSET_STATE
};


class SubgraphEnum{
public:
    Graph* data_graph_;
    Graph* query_graph_;
    vector<Vertex> order_;
    vector<Vertex> order_index_;

    uint32_t time_limit_; // seconds

    // enumeration metrics
    uint64_t emb_count_;
    long long state_count_;
    double enumeration_time_;
    double preprocessing_time_;
    double order_adjust_time_;
    double ordering_time_;
    double query_time_;
    float peak_memory_;
    vector<vector<Vertex>> matches_;

    uint64_t intersection_count_ = 0;

    uint32_t compressed_;
    uint32_t compressed_neq_;

    vector<uint64_t> leaf_states_counter_;
    vector<uint64_t> leaf_states_depth_;

    SubgraphEnum(Graph* data_graph);

    void match(Graph* query_graph, string ordering_method, uint64_t count_limit, uint32_t time_limit);

    void match_homo(Graph* query_graph, string ordering_method, uint64_t count_limit, uint32_t time_limit);

private:
    catalog* storage_;
    preprocessor* pp_;
    uint32_t query_vertex_count_, data_vertex_count_;

    bool stop_;
    
    vector<bitset<MAX_QUERY_SIZE>> ancestors_depth_;
    vector<vector<uint32_t>> successor_neighbors_in_depth_;
    vector<vector<uint32_t>> predecessor_neighbors_in_depth_;

    CandidatesHistoryStack* history_candidates_;
    vector<vector<pair<Vertex*, uint32_t>>> candidates_stack_;
    vector<vector<Vertex>> searching_candidates_;
    Vertex* candidates_offset_;
    vector<Vertex> embedding_depth_;
    Vertex* visited_query_depth_;
    bool* find_matches_;

    // to be removed
    // vector<vector<vector<Vertex>>> candidates_stack_;
    // vector<vector<Vertex>> searching_candidates_;
    // vector<FailingSet> search_failing_set_recorder_;
    bitset<MAX_QUERY_SIZE> full_descendent_; // for extended subtree reduction
    vector<vector<bitset<MAX_QUERY_SIZE>>> parent_failing_set_map_;
    // Vertex* candidates_offset_;
    // vector<Vertex> embedding_depth_;
    // Vertex* visited_query_depth_;
    // bool* find_matches_;


    bitset<MAX_QUERY_SIZE> check_neighbor_conflict_;
    unordered_set<Vertex> neighbor_candidates_;
    vector<vector<vector<Vertex>>> unmatched_group;
    vector<vector<vector<Vertex>>> matched_group;
    vector<Vertex> conflict_checking_order;

    vector<vector<vector<Vertex>>> successor_with_same_label_;

    bool check_conflict(uint32_t cur_depth);
    bool find_conflict(uint32_t cur_depth, uint32_t depth);
    vector<unordered_map<Label, vector<uint32_t>>> same_label_vertices_;
    unordered_set<Vertex> used_set_;
    unordered_map<Vertex, Vertex> reverse_match_;
    vector<vector<bitset<MAX_QUERY_SIZE>>> failing_set_stack_;
    vector<uint32_t> conflict_vec_;
    
    void initialization();

    void order_adjustment();

    void debug_print(int depth);
};