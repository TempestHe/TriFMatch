#include "candidates_recorder.h"
#include "../configuration/config.h"
#include "../graph/graph.h"
#include "../utility/utils.h"
#include "trie_encoder.h"
#include "candidates_recorder.h"

#define DEBUG 1

class ContainmentFilter{
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
    vector<vector<uint32_t>> successor_offset_map_;

    vector<int> successors_with_same_labels_;

    vector<vector<uint32_t>> empty_record_offset_;
    vector<vector<uint32_t>> conflict_record_offset_;
    vector<unordered_map<Vertex, FailingSet*>> failing_set_recorder_;

    vector<pair<vector<uint32_t>, vector<uint32_t>>> unmatched_query_depth_dist_; // first->successors, second->none_successors
#if DEBUG==1
    vector<pair<vector<Vertex>, vector<Vertex>>> same_label_depth_; // first->successors, second->none_successors;
    vector<vector<pair<vector<int>, vector<int>>>> preprocessed_depth_intersection; //(the depth v_1 resides, current_depth)==>(first->matched_depth that do not connected u_c, second->unmatched...) 
#endif
    ContainmentFilter(Graph* query_graph_, Graph* data_graph_, vector<vector<uint32_t>>& successors_depth, vector<vector<uint32_t>>& predecessor_depth, vector<Vertex>& search_order, vector<Vertex>& order_index, vector<bitset<MAX_QUERY_SIZE>>& ancestors_depth, TrieEncoder* encoder);

    bool validate_exclusion_containment(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder);
    bool validate_complete_containment(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder, vector<Vertex>& embedding_depth);

    bool validate_vector_containment_strict(Vertex* vec1, int32_t size1, Vertex* vec2, int32_t size2);
    bool validate_vector_containment_relaxed(const Vertex* vec1, int32_t vec1_size, const Vertex* vec2, int32_t vec2_size, Vertex r1, Vertex r2);
    bool validate_vector_containment_relaxed_full(Vertex* vec1, int32_t size1, Vertex* vec2, int32_t size2, Vertex r1, Vertex r2);
    bool validate_vector_element_existence(Vertex* vec1, int32_t size1, Vertex v);
    bool validate_containment_for_homo(uint32_t depth, vector<vector<pair<Vertex*, uint32_t>>>& candidates_stack, Vertex root, CandidatesHistoryStack& recorder);
};